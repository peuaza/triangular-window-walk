# -*- coding: utf-8 -*-
"""Independent check of the load-bearing certificate fact, plus a direct
measurement of the one-way communication complexity across the middle cut.

(1) K_global = n-1 ?   Monotonicity: a superset of a certificate is a
    certificate, so "some S of size k works" is monotone in k.  Hence it is
    enough to test every S with |S| = n-2: if none of them certifies the whole
    ensemble then K_global = n-1 (and n-1 always works: the ensemble's count
    constraint forces the missing cell).
(2) one-way CC across the middle cut = ceil(log2 #distinct rows) of the delta
    matrix (rows = left halves, columns = right halves).
"""
import itertools
import sys
from collections import defaultdict

import numpy as np

sys.path.insert(0, '.')
from p0_tapes import enum_tapes          # noqa: E402


def walk(a, n, rule='longrun'):
    """exact copy of _arr_rules.walk semantics (four-move table, cut-on-exit)"""
    L, R = 0, n - 1
    ell, r = 0, n - 1
    seen = set()
    steps = 0
    while True:
        if ell >= r:
            return ell / (n - 1), steps
        if R - L + 1 < 3:
            return None, steps
        steps += 1
        if steps > 20 * n * n:
            return None, steps
        if (ell, r) in seen:
            return None, steps
        seen.add((ell, r))
        z = a[ell:r + 1]
        m = r - ell + 1
        k = max(1, m // 2)
        b2 = int(z.sum()) * k > int(z[:k].sum()) * m
        if rule == 'longrun':
            b = (z > 1).astype(np.int8)
            if b.any():
                d = np.diff(np.concatenate((np.zeros(1, np.int8), b, np.zeros(1, np.int8))))
                b1 = int((np.flatnonzero(d == -1) - np.flatnonzero(d == 1)).max()) * 4 > m
            else:
                b1 = False
        else:
            raise ValueError(rule)
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
            return (ell + dy) / (n - 1), steps
        if r + dx < L or ell + dy > R:
            L, R = L + 1, R - 1
            ell, r = L, R
            seen.clear()
            continue
        ell, r = ell + dy, r + dx


def deltas(n):
    out = {}
    for a in enum_tapes(n):
        t = tuple(int(v) for v in a)
        d, _ = walk(np.array(t, dtype=np.int64), n)
        out[t] = None if d is None else round(d, 12)
    return out


def key(a, S):
    return tuple(v for i, v in enumerate(a) if i not in S)


def test_n_minus_2(tab, n):
    """is there an S of size n-2 that certifies every tape?"""
    hits = []
    for S in itertools.combinations(range(n), n - 2):
        g = defaultdict(set)
        for a, d in tab.items():
            g[key(a, S)].add(d)
        if all(len(v) == 1 for v in g.values()):
            hits.append(S)
    return hits


def owcc(tab, n, k):
    rows = defaultdict(dict)
    for a, d in tab.items():
        rows[a[:k]][a[k:]] = d
    sigs = defaultdict(int)
    for l, col in rows.items():
        sigs[tuple(sorted(col.items()))] += 1
    return len(sigs), len(rows), 3 ** k


if __name__ == '__main__':
    for n in [int(x) for x in sys.argv[1:]] or [6, 8, 10]:
        tab = deltas(n)
        hits = test_n_minus_2(tab, n)
        nd = sum(1 for v in tab.values() if v is not None)
        print('n=%2d tapes=%6d defined=%d (%.1f%%)  certifying S of size n-2: %d  => K_global = %s'
              % (n, len(tab), nd, 100.0 * nd / len(tab), len(hits),
                 (n - 2) if hits else (n - 1)))
        h = n // 2
        nr, nrows, maxrows = owcc(tab, n, h)
        print('      one-way cut at %d: %d distinct rows / %d left halves (3^%d=%d) -> CC <= %.2f bits;'
              '  n/2*log2(3)=%.2f;  bound T>=CC/s with s=log2 n: %.2f'
              % (h, nr, nrows, h, maxrows, np.log2(nr), h * np.log2(3),
                 np.log2(nr) / np.log2(n)))
        sys.stdout.flush()
