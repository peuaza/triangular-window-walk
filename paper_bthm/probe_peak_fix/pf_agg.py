# -*- coding: utf-8 -*-
"""pf_agg -- aggregate the pf_repro segments into the EN_LINEAR 1 500-tape grid.

The 1 500 tapes of EN_LINEAR's peak claim are
    main grid      : n = 512, 1024, 2048 (200 each) + 4096 (100)      =  700
    replay 512     : i = 900 000, 200                                  =  200
    replay 1024    : i = 500 000, 200                                  =  200
    replay 2048    : i = 700 000, 200                                  =  200
    replay 2048    : i = 300 000, 200                                  =  200
                                                                    -------
                                                                     = 1500
The 4096 replay (i = 900 000, 100) is outside the 1 500 and is reported apart.
"""
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
SEGS = [
    ('512 main', 512, 0, 200),
    ('1024 main', 1024, 0, 200),
    ('2048 main', 2048, 0, 200),
    ('4096 main', 4096, 0, 100),
    ('512 repl', 512, 900000, 200),
    ('1024 repl', 1024, 500000, 200),
    ('2048 repl A', 2048, 700000, 200),
    ('2048 repl B', 2048, 300000, 200),
    ('4096 repl', 4096, 900000, 100),
]


def fname(n, off):
    return os.path.join(HERE, 'seg_%d_%d.jsonl' % (n, off))


def load():
    data = []
    for name, n, off, k in SEGS:
        p = fname(n, off)
        if not os.path.exists(p):
            print('MISSING', p)
            continue
        rs = [json.loads(l) for l in open(p)]
        data.append((name, n, off, k, rs))
    return data


def main():
    data = load()
    tot = {p: 0 for p in ('base128', 'fix128', 'base256', 'fix256')}
    ntot = 0
    print('%-12s %5s %6s | %-16s %-16s | %-16s %-16s'
          % ('segment', 'n', 'tapes', 'base W2=128', 'fix W2=128',
             'base W2=256', 'fix W2=256'))
    for name, n, off, k, rs in data:
        cells = []
        for p in ('base128', 'fix128', 'base256', 'fix256'):
            ex = sum(r[p]['exact'] for r in rs)
            e = np.mean([r[p]['reads'] for r in rs]) / n
            cells.append('%3d/%3d %.2f' % (ex, len(rs), e))
        print('%-12s %5d %6d | %-16s %-16s | %-16s %-16s'
              % (name, n, len(rs), cells[0], cells[1], cells[2], cells[3]))
        if name != '4096 repl':
            ntot += len(rs)
            for p in tot:
                tot[p] += sum(r[p]['exact'] for r in rs)
    print('\n1 500-tape total: base128 %d/%d  fix128 %d/%d  base256 %d/%d  fix256 %d/%d'
          % (tot['base128'], ntot, tot['fix128'], ntot,
             tot['base256'], ntot, tot['fix256'], ntot))

    print('\nfailures (W2=256, the EN_LINEAR exactness machine):')
    for name, n, off, k, rs in data:
        for r in rs:
            if not r['base256']['exact']:
                b = r['bad'][0] if r['bad'] else None
                print('  %-11s i=%-7d base256 end=%d delta=%.6f  firstbad=%s  nbad=%d'
                      % (name, r['i'], r['base256']['end'],
                         r['base256']['delta'], b, r['nbad']))
                print('       top:', r['top'][:10])
    print('\nfailures at W2=128 (the 2.0 n machine):')
    for name, n, off, k, rs in data:
        for r in rs:
            if not r['base128']['exact']:
                print('  %-11s i=%-7d fd=(k=%d,M0=%d,pred=%d,out=%d) fix128 exact=%d'
                      % (name, r['i'], r['base128']['fd'], r['base128']['fdM0'],
                         r['base128']['fdpred'], r['base128']['fdout'],
                         r['fix128']['exact']))
    print('\nfailures of the fixed machine:')
    for name, n, off, k, rs in data:
        for r in rs:
            if not r['fix256']['exact'] or not r['fix128']['exact']:
                print('  %-11s i=%-7d fix128=%d fix256=%d firstbad_M0=%s'
                      % (name, r['i'], r['fix128']['exact'], r['fix256']['exact'],
                         r['bad'][0][1] if r['bad'] else None))

    # per-size (the step-3 table): 200 tapes at each size
    print('\nper-size, mode of the fixed machine:')
    for n in (512, 1024, 2048, 4096):
        rows = []
        for name, nn, off, k, rs in data:
            if nn == n and name != '4096 repl':
                rows += rs
        if not rows:
            continue
        for W in (128, 256):
            p = 'fix%d' % W
            ex = sum(r[p]['exact'] for r in rows)
            e = np.mean([r[p]['reads'] for r in rows]) / n
            print('  n=%4d tapes=%3d  fix W2=%3d exact %3d/%3d  e/n %.3f'
                  % (n, len(rows), W, ex, len(rows), e))


if __name__ == '__main__':
    main()
