# -*- coding: utf-8 -*-
"""el_gate2 -- gates for the jump machine itself.

A. With W2 >= n the machine never jumps: it must reproduce the trace's
   (delta, end) on every rule.  That validates the endgame's control flow.
B. The lazy range-max device must produce the same answer as the scan device
   (same b1 at every epoch), on MAXP/AMAX/AMIN.
C. The sliding device must be exact-in-principle: at W2 = n it is the literal
   walk, so the same check as A covers it.

Usage: python el_gate2.py     -> el_gate2.json
"""
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import el_cell as C
import el_lib as E

NS = (96, 160)
NT = 8
OUT = os.path.join(HERE, 'el_gate2.json')


def main():
    fails = []
    rows = []
    for rule in E.RULES:
        mode = ('sliding' if E.EVAL[rule] == 'sliding'
                else ('lazy' if E.EVAL[rule] == 'rangemax' else 'scan'))
        for n in NS:
            bad = 0
            for i in range(NT):
                X = E.tape(n, i)
                M = E.Mach(X, n)
                t = M.trace(rule)
                d, en = C.jm_walk(rule, M.S, n, n, mode)[:2]
                ok = (int(en) == int(t['end'])) and (
                    (d is None and t['delta'] < 0) or
                    (d is not None and t['delta'] >= 0
                     and abs(d - t['delta']) < 1e-12))
                if not ok:
                    bad += 1
                    fails.append([rule, n, i, 'W2=n', d, en, t['delta'], t['end']])
            rows.append(dict(rule=rule, n=n, gate='W2=n', bad=bad))
            print('%-9s n=%4d  W2=n mismatches %d/%d' % (rule, n, bad, NT))
    # B: lazy == scan
    for rule in E.RANGE_MAX:
        for n in NS:
            bad = 0
            for i in range(NT):
                X = E.tape(n, i)
                M = E.Mach(X, n)
                a = C.jm_walk(rule, M.S, n, 32, 'lazy')
                b = C.jm_walk(rule, M.S, n, 32, 'scan')
                if (a[0] != b[0]) or (a[1] != b[1]) or (a[3] != b[3]) or (a[4] != b[4]):
                    bad += 1
                    fails.append([rule, n, i, 'lazy-scan', a[:5], b[:5]])
            rows.append(dict(rule=rule, n=n, gate='lazy==scan', bad=bad))
            print('%-9s n=%4d  lazy != scan on %d/%d' % (rule, n, bad, NT))
    ok = not fails
    with open(OUT, 'w') as fh:
        json.dump(dict(ok=bool(ok), rows=rows, fails=fails[:20]), fh, indent=1)
    print('GATE2', 'PASS' if ok else 'FAIL')
    return 0 if ok else 1


if __name__ == '__main__':
    sys.exit(main())
