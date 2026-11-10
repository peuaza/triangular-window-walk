# -*- coding: utf-8 -*-
"""el_lazydbg -- find the first epoch where the lazy b1 disagrees with b1_at."""
import sys

import numpy as np

sys.path.insert(0, '.')
import el_cell as C
import el_lib as E
import pcr_lib as P

rule = sys.argv[1]
n = int(sys.argv[2])
i = int(sys.argv[3])
X = E.tape(n, i)
S = P.St(X)
a = S.a
Ps = S.Ps

L, R = 0, n - 1
ismax = rule in ('MAXP', 'AMAX')
qhi = (R + 1) if rule == 'MAXP' else R
s = 0
best = 0
p = L
for j in range(L, qhi + 1):
    if j > L:
        s += int(a[j - 1])
    if (s > best) if ismax else (s < best):
        best = s
        p = j
pv = best
Ssum = sum(int(v) for v in a)
k = 0
while R - L + 1 > 32:
    M = R - L + 1
    q = (R + 1) if rule == 'MAXP' else R
    if p < L or p > q:
        s = 0
        best = 0
        p = L
        for j in range(L, q + 1):
            if j > L:
                s += int(a[j - 1])
            if (s > best) if ismax else (s < best):
                best = s
                p = j
        pv = best
    if rule == 'MAXP':
        b1 = 1 if pv * pv > M else 0
    elif rule == 'AMAX':
        b1 = 1 if Ssum > pv else 0
    else:
        b1 = 1 if Ssum < pv else 0
    bt = C.b1_at(S, rule, L, R)
    if b1 != bt:
        print('MISMATCH ep %d region [%d,%d] M %d lazy %d true %d  p %d pv %d '
              'Ssum %d q %d' % (k, L, R, M, b1, bt, p, pv, Ssum, q))
        # what is the true max?
        vals = [int(Ps[j]) - int(Ps[L]) for j in range(min(L, q), max(q, L) + 1)]
        print('   Ps offsets over the query range:', vals[:20], '...' if len(vals) > 20 else '')
        break
    pr = C.pred_of(rule, C.b2_at(S, L, R), b1)
    k += 1
    if pr == 3:
        dw = int(a[L])
        Ssum -= dw + int(a[R])
        pv -= dw
        L += 1
        R -= 1
    elif pr == 2:
        Ssum -= int(a[R])
        R -= 1
    else:
        dw = int(a[L])
        Ssum -= dw
        pv -= dw
        L += 1
else:
    print('no mismatch in', k, 'epochs')
