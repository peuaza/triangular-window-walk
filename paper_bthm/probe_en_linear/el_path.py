# -*- coding: utf-8 -*-
"""el_path -- compare the machine's region path with the walk's, epoch by epoch.

Usage: python el_path.py RULE N TAPE_INDEX W2
"""
import sys

import numpy as np

sys.path.insert(0, '.')
import el_cell as C
import el_lib as E

rule = sys.argv[1]
n = int(sys.argv[2])
i = int(sys.argv[3])
W2 = int(sys.argv[4])

X = E.tape(n, i)
M = E.Mach(X, n)
t = M.trace(rule)
# truth path
L, R = 0, n - 1
truth = []
for k in range(t['nep']):
    truth.append((L, R, int(t['epM'][k]), int(t['epb2'][k]), int(t['epb1'][k]),
                  int(t['epout'][k]), int(t['eplen'][k])))
    o = int(t['epout'][k])
    if o == 1:
        L += 1
    elif o == 2:
        R -= 1
    elif o == 3:
        L += 1
        R -= 1
print('truth: end %d delta %.6f nep %d steps %d' % (t['end'], t['delta'],
                                                    t['nep'], t['steps']))
# machine path
S = M.S
L, R = 0, n - 1
k = 0
while True:
    Mw = R - L + 1
    if Mw < 3:
        print('machine end 1 at region [%d,%d]' % (L, R))
        break
    if Mw <= W2:
        print('  -- endgame at [%d,%d] (epoch %d)' % (L, R, k))
        tt = truth[k] if k < len(truth) else None
        if tt:
            print('     truth at this epoch: [%d,%d] M %d b2 %d b1 %d out %d'
                  % tt[:6])
        d, en = C.jm_walk(rule, S, n, W2, 'sliding')[:2]
        print('     machine result delta %s end %d' % (d, en))
        break
    b1 = C.b1_at(S, rule, L, R)
    b2 = C.b2_at(S, L, R)
    p = C.pred_of(rule, b2, b1)
    tt = truth[k] if k < len(truth) else None
    flag = ''
    if tt and tt[5] != p:
        flag = '   <-- MISMATCH (truth %d)' % tt[5]
    print('  ep %4d region [%4d,%4d] M %4d b2 %d b1 %d pred %d%s'
          % (k, L, R, Mw, b2, b1, p, flag))
    k += 1
    if p == 3:
        L += 1
        R -= 1
    elif p == 2:
        R -= 1
    else:
        L += 1
    if k > n + 4:
        print('runaway')
        break
