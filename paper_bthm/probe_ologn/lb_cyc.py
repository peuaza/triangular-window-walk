# -*- coding: utf-8 -*-
"""Why does the walk return 'no halt'?  And how long are the cycles?

The write-up wants to know whether the O(log n)-bit simulation story has to
carry cycle detection (which would cost more than O(log n) bits) or whether the
non-halting cases are short cycles that O(1) extra words can spot.

Classification per tape: meet / region<3 / cycle (period recorded) / cap.
Same control flow as lb_cert2.walk.  Usage: python lb_cyc.py 6 8 10 11 12
"""
import sys
from collections import Counter

import numpy as np

sys.path.insert(0, '.')
from p0_tapes import enum_tapes          # noqa: E402


def walk_class(a, n):
    L, R = 0, n - 1
    ell, r = 0, n - 1
    seen = {}
    steps = 0
    while True:
        if ell >= r:
            return 'meet', 0, steps
        if R - L + 1 < 3:
            return 'short', 0, steps
        steps += 1
        if steps > 20 * n * n:
            return 'cap', 0, steps
        key = (ell, r)
        if key in seen:
            return 'cycle', steps - seen[key], steps
        seen[key] = steps
        z = a[ell:r + 1]
        m = r - ell + 1
        k = max(1, m // 2)
        b2 = int(z.sum()) * k > int(z[:k].sum()) * m
        b = (z > 1).astype(np.int8)
        if b.any():
            d = np.diff(np.concatenate((np.zeros(1, np.int8), b, np.zeros(1, np.int8))))
            b1 = int((np.flatnonzero(d == -1) - np.flatnonzero(d == 1)).max()) * 4 > m
        else:
            b1 = False
        if b2 and b1:
            dy, dx = 1, 0
        elif b2 and not b1:
            dy, dx = 0, -1
        elif b1 and not b2:
            dy, dx = 0, 1
        else:
            dy, dx = -1, 0
        if ell + dy < L or r + dx > R:
            if ell == L and r == R:
                L, R = L + 1, R - 1
            elif ell + dy < L:
                L += 1
            else:
                R -= 1
            ell, r = L, R
            seen.clear()
            continue
        if ell + dy >= r + dx:
            return 'meet', 0, steps
        if r + dx < L or ell + dy > R:
            L, R = L + 1, R - 1
            ell, r = L, R
            seen.clear()
            continue
        ell, r = ell + dy, r + dx


if __name__ == '__main__':
    for n in [int(x) for x in sys.argv[1:]] or [8, 10]:
        cls = Counter()
        per = Counter()
        maxper = 0
        for a in enum_tapes(n):
            a = np.asarray(a, dtype=np.int64)
            c, p, s = walk_class(a, n)
            cls[c] += 1
            if c == 'cycle':
                per[p] += 1
                maxper = max(maxper, p)
        tot = sum(cls.values())
        print('n=%2d tapes=%6d  meet %.1f%%  short %.1f%%  cycle %.1f%%  cap %.1f%%'
              % (n, tot, 100.0 * cls['meet'] / tot, 100.0 * cls['short'] / tot,
                 100.0 * cls['cycle'] / tot, 100.0 * cls['cap'] / tot))
        print('       cycle periods: %s   max period %d'
              % (dict(sorted(per.items())), maxper))
        sys.stdout.flush()
