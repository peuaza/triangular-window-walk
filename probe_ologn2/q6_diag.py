# -*- coding: utf-8 -*-
"""q6: diagnose where the four-quarters sufficient condition breaks.

For every M and every quarter-monotone tape (exhaustive over 3^M, M <= 11),
find the LARGEST m with b2 = 0 and report min over tapes of m*/M, plus the
histogram of the failing m at the boundary.
"""
import itertools
from collections import Counter
import numpy as np


def quarters(M):
    q, r = divmod(M, 4)
    lens = [q + (1 if i < r else 0) for i in range(4)]
    b = [0]
    for x in lens:
        b.append(b[-1] + x)
    return lens, b


def b2(a, L, m):
    k = max(1, m // 2)
    s = int(a[L:L + m].sum())
    sk = int(a[L:L + k].sum())
    return int(s * k > sk * m)


for M in range(5, 12):
    lens, b = quarters(M)
    fails = Counter()
    worst = 1.0
    tested = 0
    for tup in itertools.product((0, 1, 2), repeat=M):
        a = np.array(tup, dtype=np.int64)
        qm = [a[b[i]:b[i + 1]].mean() for i in range(4)]
        if not (qm[0] < qm[1] < qm[2] < qm[3]):
            continue
        tested += 1
        z = [m for m in range(2, M + 1) if b2(a, 0, m) == 0]
        mstar = max(z) if z else None
        if mstar is not None:
            fails[mstar] += 1
            worst = min(worst, mstar / M)
    print('M=%2d lens=%s  monotone tapes=%6d  largest-concavity histogram %s   min m*/M=%.2f'
          % (M, lens, tested, dict(sorted(fails.items())), worst))
