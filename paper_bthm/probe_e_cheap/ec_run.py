"""Re-measure e with the cheap exact machine.

`probe_en_linear` §8: dropping the re-span from the pcr walk gives a machine that is exact on every
tape (meeting or halting) whose cost is n + steps + cuts cells.  The record's e is the cost of the
re-spanning walk, so it is systematically high.  This script measures both on the same tapes.

e_old  mean of the record machine's reads over the meeting tapes (the record's number)
e_new  n + mean(steps) + mean(cuts) over the same tapes

Usage: python ec_run.py [n ...]
"""
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
PAPER = os.path.dirname(HERE)
TAN = os.path.join(PAPER, 'probe_triple_nat')
for d in (TAN, os.path.join(PAPER, 'probe_e_uniform')):
    if d not in sys.path:
        sys.path.insert(0, d)

import tr_lib as T                                              # noqa: E402
import eu_lib as L                                              # noqa: E402

RULES = ('dense2', 'peak', 'trip', 'oddrun', 'canonical', 'pal', 'hashmod', 'MAXP', 'LAS',
         'KCHG', 'BAL>m/2', 'BAL>m/4', 'NBLK=1', 'NBLK>2', 'NBLK>3', 'AMAX', 'AMIN')
# EN_LINEAR §8: n + steps + cuts is a legal cost only for a sliding statistic, i.e. one whose state
# is O(log n) bits and is carried across the cut.  For a window extremum the statistic needs a scan
# per step, so the row would be a step count, not a cost.
LEGAL = {'dense2', 'peak', 'trip', 'oddrun', 'canonical', 'pal', 'hashmod', 'LAS', 'KCHG',
         'BAL>m/2', 'BAL>m/4', 'NBLK=1', 'NBLK>2', 'NBLK>3'}
NT = 200
CAP = 6000


def cell(rule, n):
    old, new = [], []
    drawn = 0
    for i in range(CAP):
        if len(new) >= NT:
            break
        drawn = i + 1
        a, s, X = L.make_tape(n, 5000 + n + 100003 * i)
        r = T.run(rule, a, s)
        d = r.get('delta')
        if d is None or d != d:
            continue
        st = r.get('steps')
        cu = r.get('cuts')
        if st is None:
            continue
        cu = 0 if cu is None else int(cu)
        old.append(float(r['reads']))
        new.append(float(n + int(st) + cu))
    if not new:
        return None
    return dict(rule=rule, n=n, nt=len(new), drawn=drawn, legal=rule in LEGAL,
                e_old=float(np.mean(old)), e_new=float(np.mean(new)),
                e_old_n=float(np.mean(old)) / n, e_new_n=float(np.mean(new)) / n)


def main(ns):
    out = []
    for n in ns:
        for rule in RULES:
            c = cell(rule, n)
            if c is None:
                print(f'{rule:<10} n={n:>5}  no machine', flush=True)
                continue
            out.append(c)
            print(f"{rule:<10} n={n:>5} nt={c['nt']:>4} e_old/n={c['e_old_n']:>8.2f}"
                  f" e_new/n={c['e_new_n']:>8.2f} ratio={c['e_old']/c['e_new']:>5.2f}", flush=True)
        with open(os.path.join(HERE, 'ec_cells.json'), 'w') as f:
            json.dump(out, f)
    print('written', len(out), 'cells')


if __name__ == '__main__':
    ns = [int(x) for x in sys.argv[1:]] or [512, 1024, 2048]
    main(ns)
