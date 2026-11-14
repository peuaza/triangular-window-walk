# -*- coding: utf-8 -*-
"""par_all -- every cell of the BTHM-C_2 measurement, resumable.

usage: par_all.py [n ...] [--jobs K]

Cells whose par_<rule>_<n>.json exists are skipped.  One process per cell.
"""
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)


def job(args):
    n, rule, nt = args
    import par_run
    t0 = time.time()
    f, out = par_run.one(rule, n, nt)
    r = out['rows']
    e0 = [x['end0'] for x in r]
    eP = [x['endP'] for x in r]
    ntr = [x['ntrim'] for x in r]
    med = sorted(ntr)[len(ntr) // 2] if ntr else -1
    return ('%-10s n=%-5d meet0=%3d cyc0=%3d | par: meet=%3d empty=%3d '
            'cap=%3d full=%3d trims med=%-6d max=%-6d [%.0fs]'
            % (rule, n, e0.count(0), e0.count(4), eP.count(0), eP.count(1),
               eP.count(3) + eP.count(6), eP.count(5), med, max(ntr),
               time.time() - t0))


def main():
    argv = sys.argv[1:]
    jobs = 8
    if '--jobs' in argv:
        k = argv.index('--jobs')
        jobs = int(argv[k + 1])
        del argv[k:k + 2]
    ns = [int(a) for a in argv if a.isdigit()] or [512, 1024, 2048]
    import par_engine                                             # noqa: F401
    from tr_lib import RULENAME, NT
    rules = [r for r in RULENAME if r != 'monotone']
    tasks = []
    for n in ns:
        for rule in rules:
            f = os.path.join(HERE, 'par_%s_%d.json' % (
                rule.replace('>', 'g').replace('=', 'e').replace('/', 'd'), n))
            if not os.path.exists(f):
                tasks.append((n, rule, NT))
    print('tasks: %d, jobs: %d' % (len(tasks), jobs), flush=True)
    if not tasks:
        return
    # warm up so workers load the numba cache instead of racing to write it
    import numpy as np
    import par_engine as E
    X = np.zeros((8, 2), np.int8)
    X[0, 0] = X[1, 1] = 1
    C = E.Cyc(X, hlog=20)
    C.run(0, 0, 0, mode=0)
    C.run(0, 0, 0, mode=5)
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
