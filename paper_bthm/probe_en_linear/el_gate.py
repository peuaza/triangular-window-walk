# -*- coding: utf-8 -*-
"""el_gate -- the trace must reproduce pcr_lib.Runner bit for bit.

For every rule and every tape at n = 48, 64, 128 (the record's tapes i = 0..11)
compare delta, end, steps, cuts, inc reads, msum.  Also check the trace's own
bookkeeping: sum(eplen) == steps, nep == cuts + 1, and the last epoch's outcome
agrees with the walk's end code.

Usage: python el_gate.py
Writes el_gate.json
"""
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import el_lib as E
import pcr_lib as P

NS = (48, 64, 128)
NT = 12
OUT = os.path.join(HERE, 'el_gate.json')


def ref_run(rule, X):
    if rule == 'BAL>m/4':
        import tr_bal4
        r = tr_bal4.Runner(X).run(thr=4)
        return dict(delta=r['delta'], end=int(r['end']), steps=int(r['steps']),
                    msum=int(r['scan']), inc=int(r['reads']), cuts=int(r['cuts']))
    if rule == 'NBLK>3':
        # pcr_lib has no NBLK>3 row; hmach3 is an independent implementation of
        # the same walk for that rule (its own gate is hgate.py).  Only the walk
        # output is compared -- steps/cuts mean different things there.
        import hmach3
        X = np.asarray(X)
        s = (X[:, 0] - X[:, 1]).astype(np.int64)
        a = (X[:, 0] + X[:, 1]).astype(np.int64)
        r = hmach3.H3([int(v) for v in s], [int(v) for v in a], 'nblk_gt3',
                      B=64, log=False, deep=False).run()
        dv = r['delta']
        ec = {'meet': 0, 'narrow': 1, 'empty': 1, 'cycle': 2, 'repeat': 2,
              'cap': 3}.get(r['end'], 9)
        return dict(delta=(-1.0 if dv is None else float(dv)), end=int(ec))
    r = P.Runner(X).run(rule)
    return dict(delta=r[0], end=int(r[5]), steps=int(r[1]), msum=int(r[2]),
                inc=int(r[3]), cuts=int(r[4]))


def main():
    bad = []
    rows = []
    for rule in E.RULES:
        for n in NS:
            diffs = 0
            book = 0
            for i in range(NT):
                X = E.tape(n, i)
                a_ref = ref_run(rule, X)
                M = E.Mach(X, n)
                t = M.trace(rule)
                a = dict(delta=t['delta'], end=t['end'], steps=t['steps'],
                         msum=t['msum'], inc=t['inc'], cuts=t['cuts'])
                for k in a:
                    if k not in a_ref:
                        continue
                    if k == 'delta':
                        x, y = a[k], a_ref[k]
                        if not ((x == y) or (abs(x - y) < 1e-12)):
                            diffs += 1
                            bad.append((rule, n, i, k, x, y))
                    elif a[k] != a_ref[k]:
                        diffs += 1
                        bad.append((rule, n, i, k, a[k], a_ref[k]))
                # bookkeeping: the epoch lengths must add up to the step count
                # (a cycle-detecting step is counted in `steps` but belongs to
                # no epoch, so a cycle tape may be one short)
                sl = int(t['eplen'].sum())
                if not (sl == t['steps'] or (t['end'] == 2 and sl == t['steps'] - 1)):
                    book += 1
                if t['nep'] != t['cuts'] + 1:
                    book += 1
                last = int(t['epout'][-1])
                if t['end'] == 0 and last != 0:
                    book += 1
                if t['end'] == 1 and last not in (1, 2, 3, 4):
                    book += 1
                if t['end'] == 2 and last != 5:
                    book += 1
            rows.append(dict(rule=rule, n=n, nt=NT, diffs=diffs, book=book))
            print('%-9s n=%4d  diffs=%d book=%d' % (rule, n, diffs, book))
    ok = (not bad) and all(r['book'] == 0 for r in rows)
    d = dict(ns=NS, nt=NT, ok=bool(ok), rows=rows,
             fails=[list(b) for b in bad[:20]])
    with open(OUT, 'w') as fh:
        json.dump(d, fh, indent=1)
    print('GATE', 'PASS' if ok else 'FAIL', '->', OUT)
    return 0 if ok else 1


if __name__ == '__main__':
    sys.exit(main())
