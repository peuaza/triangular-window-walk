# -*- coding: utf-8 -*-
"""cyc_all -- every cell, resumable.

usage: cyc_all.py [n ...] [--trim1] [--trim2]
Cells whose JSON exists are skipped, so a crash loses one cell.
"""
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)


def job(args):
    n, rule, nt, trim = args
    import cyc_run
    t0 = time.time()
    f, out = cyc_run.one(rule, n, nt, trim)
    r = out['rows']
    e0 = [x['end0'] for x in r]
    return '%s n=%-5d trim=%d  meet0=%3d cyc=%3d empty0=%3d meetC=%3d [%.0fs]' % (
        rule, n, trim, e0.count(0), e0.count(4), e0.count(1),
        [x['endC'] for x in r].count(0), time.time() - t0)


def main():
    argv = sys.argv[1:]
    trims = [1]
    if '--trim2' in argv:
        trims = [2]
    if '--trim3' in argv:
        trims = [3]
    if '--trim5' in argv:
        trims = [5]
    if '--both' in argv:
        trims = [1, 2]
    ns = [int(a) for a in argv if a.isdigit()] or [512, 1024, 2048]
    import cyc_engine                                              # noqa: F401  (sets sys.path)
    from tr_lib import RULENAME, NT
    tasks = []
    for n in ns:
        for rule in RULENAME:
            if rule == 'monotone':
                continue
            for trim in trims:
                import cyc_run
                f = os.path.join(HERE, 'cyc_%s_%d%s.json' % (
                    cyc_run.safe(rule), n, '' if trim == 1 else '_trim%d' % trim))
                if not os.path.exists(f):
                    tasks.append((n, rule, NT, trim))
    print('tasks: %d' % len(tasks), flush=True)
    if not tasks:
        return
    # warm up: compile the walk once in the parent so the workers load the numba
    # cache instead of racing to write it
    import numpy as np
    import cyc_engine as E
    X = np.zeros((8, 2), np.int8)
    X[0, 0] = X[1, 1] = 1
    C = E.Cyc(X)
    C.run(0, 0, 0, mode=0)
    C.run(0, 0, 0, mode=1, trim=2)
    C.run(0, 0, 2, mode=0)
    print('warmup done', flush=True)
    from concurrent.futures import ProcessPoolExecutor
    t0 = time.time()
    with ProcessPoolExecutor(max_workers=4) as ex:
        for k, line in enumerate(ex.map(job, tasks)):
            print('%3d/%3d %s  total %.0fs' % (k + 1, len(tasks), line,
                                               time.time() - t0), flush=True)
    print('done %.0fs' % (time.time() - t0), flush=True)


if __name__ == '__main__':
    main()
