# -*- coding: utf-8 -*-
"""el_litwalk -- a plain Python four-move walk, for arbitration only.

python el_litwalk.py RULE N TAPE
Prints the literal walk's (delta, end, steps, cuts) and the step-by-step
(ell, r, b1, b2) list around the first divergence from el_lib's trace.
"""
import sys

import numpy as np

sys.path.insert(0, '.')
import el_lib as E
import el_cell as C
import pcr_lib as P

TAB = P.TABLE_NP
rule = sys.argv[1]
n = int(sys.argv[2])
i = int(sys.argv[3])

X = E.tape(n, i)
S = P.St(X)
a = S.a


def lit():
    L, R = 0, n - 1
    ell, r = L, R
    seen = set()
    steps = cuts = 0
    trace = []
    while True:
        if ell >= r:
            return ell / (n - 1.0), 0, steps, cuts, trace
        if R - L + 1 < 3:
            return -1.0, 1, steps, cuts, trace
        steps += 1
        if (ell, r) in seen:
            return -1.0, 2, steps, cuts, trace
        seen.add((ell, r))
        b1 = C.b1_at(S, rule, ell, r)
        b2 = C.b2_at(S, ell, r)
        trace.append((ell, r, b1, b2))
        act = TAB[2 * b2 + b1]
        dy = 1 if act == P.LPLUS else (-1 if act == P.LMINUS else 0)
        dx = 1 if act == P.RPLUS else (-1 if act == P.RMINUS else 0)
        if ell + dy < L or r + dx > R:
            cuts += 1
            if ell == L and r == R:
                L += 1
                R -= 1
            elif ell + dy < L:
                L += 1
            else:
                R -= 1
            if L > R:
                return -1.0, 1, steps, cuts, trace
            ell, r = L, R
            seen.clear()
            continue
        if ell + dy >= r + dx:
            return (ell + dy) / (n - 1.0), 0, steps, cuts, trace
        ell += dy
        r += dx


d, e, st, cu, tr = lit()
print('literal  delta %.6f end %d steps %d cuts %d' % (d, e, st, cu))
M = E.Mach(X, n)
t = M.trace(rule)
print('trace    delta %.6f end %d steps %d cuts %d' % (t['delta'], t['end'],
                                                       t['steps'], t['cuts']))
if rule == 'BAL>m/4':
    import tr_bal4
    r2 = tr_bal4.Runner(X).run(thr=4)
    print('tr_bal4  delta %.6f end %d steps %d cuts %d'
          % (r2['delta'], r2['end'], r2['steps'], r2['cuts']))
print('literal trace (ell,r,b1,b2) head:', tr[:12])
