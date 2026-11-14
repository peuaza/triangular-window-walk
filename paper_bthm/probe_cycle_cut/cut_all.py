# -*- coding: utf-8 -*-
"""cut_all -- every cell of the cut-repair measurement, resumable.

usage: cut_all.py [n ...] [--jobs K]

Cells whose cut_<rule>_<n>.json exists are skipped.  One process per cell,
because the work is one numba loop per tape and cells are independent.
"""
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)


def job(args):
    n, rule, nt = args
    import cut_run
    t0 = time.time()
    f, out = cut_run.one(rule, n, nt)
    r = out['rows']
    e0 = [x['end0'] for x in r]
    eC = [x['endC'] for x in r]
    ntr = [x['ntrim'] for x in r]
    return ('%-10s n=%-5d meet0=%3d cyc0=%3d | cut: meet=%3d empty=%3d cap=%3d '
            'trims max=%-6d [%.0fs]'
            % (rule, n, e0.count(0), e0.count(4), eC.count(0), eC.count(1),
               eC.count(3) + eC.count(6), max(ntr), time.time() - t0))


def main():
    argv = sys.argv[1:]
    jobs = 8
    if '--jobs' in argv:
        k = argv.index('--jobs')
        jobs = int(argv[k + 1])
        del argv[k:k + 2]
    ns = [int(a) for a in argv if a.isdigit()] or [512, 1024, 2048]
    import cut_engine                                             # noqa: F401
    from tr_lib import RULENAME, NT
    rules = [r for r in RULENAME if r != 'monotone']
    tasks = []
    for n in ns:
        for rule in rules:
            f = os.path.join(HERE, 'cut_%s_%d.json' % (
                rule.replace('>', 'g').replace('=', 'e').replace('/', 'd'), n))
            if not os.path.exists(f):
                tasks.append((n, rule, NT))
    print('tasks: %d, jobs: %d' % (len(tasks), jobs), flush=True)
    if not tasks:
        return
    # warm up so workers load the numba cache instead of racing to write it
    import numpy as np
    import cut_engine as E
    X = np.zeros((8, 2), np.int8)
    X[0, 0] = X[1, 1] = 1
    C = E.Cyc(X, hlog=20)
    C.run(0, 0, 0, mode=0)
    C.run(0, 0, 0, mode=4)
    print('warmup done', flush=True)
    from concurrent.futures import ProcessPoolExecutor
    t0 = time.time()
    with ProcessPoolExecutor(max_workers=jobs) as ex:
        for k, line in enumerate(ex.map(job, tasks)):
            print('%3d/%3d %s  total %.0fs' % (k + 1, len(tasks), line,
                                               time.time() - t0), flush=True)
    print('done %.0fs' % (time.time() - t0), flush=True)


if __name__ == '__main__':
    main()
