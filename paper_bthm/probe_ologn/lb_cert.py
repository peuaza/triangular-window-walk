# -*- coding: utf-8 -*-
"""Independent check of the two load-bearing certificate facts, plus a direct
measurement of the one-way communication complexity across the middle cut.

(1) K_global = n-1 ?   By monotonicity (a superset of a certificate is a
    certificate) it suffices to test every S with |S| = n-2 : if none of them
    certifies the whole ensemble, no smaller S does either.
(2) one-way CC across the middle cut = log2(#distinct rows of the
    communication matrix), rows = left halves, columns = right halves,
    entries = delta (a sentinel for "no output").
"""
import itertools
import sys
from collections import defaultdict

import numpy as np

sys.path.insert(0, '.')
from p0_tapes import enum_tapes          # noqa: E402
from lb_steps import walk                # noqa: E402  (same control flow as _arr_rules)


def deltas(n, rule='longrun'):
    out = {}
    for a in enum_tapes(n):
        t = tuple(int(v) for v in a)
        d, _ = walk(np.array(t, dtype=np.int64), n, rule)
        out[t] = None if d is None else round(d, 12)
    return out


def key(a, S):
    return tuple(v for i, v in enumerate(a) if i not in S)


def kglobal(tab, n):
    """smallest k for which some S of size k is a certificate for every tape"""
    for k in range(1, n - 1):          # k <= n-2 by the monotonicity shortcut
        for S in itertools.combinations(range(n), k):
            g = defaultdict(set)
            for a, d in tab.items():
                g[key(a, S)].add(d)
            if all(len(v) == 1 for v in g.values()):
                return k, S
    return n - 1, None                  # n-1 always works (counts force the last cell)


def owcc(tab, n, k):
    """# distinct rows of the delta matrix, left = [0,k), right = [k,n)"""
    rows = defaultdict(dict)
    for a, d in tab.items():
        rows[a[:k]][a[k:]] = d
    sigs = defaultdict(int)
    for l, col in rows.items():
        sig = tuple(sorted(col.items()))
        sigs[sig] += 1
    return len(sigs), len(rows), 3 ** k


if __name__ == '__main__':
    for n in [int(x) for x in sys.argv[1:]] or [6, 8, 10]:
        tab = deltas(n)
        k, S = kglobal(tab, n)
        print('n=%2d  tapes=%6d  K_global=%d   (n-1=%d)  witness S=%s'
              % (n, len(tab), k, n - 1, S))
        h = n // 2
        nr, nrows, maxrows = owcc(tab, n, h)
        print('        one-way cut at %d: distinct rows=%d of %d possible left halves'
              ' (<=3^%d=%d)  =>  CC <= %.2f bits'
              % (h, nr, nrows, h, maxrows, np.log2(max(nr, 1))))
        # how much of the ensemble is defined at all
        nd = sum(1 for v in tab.values() if v is not None)
        print('        defined (walk halts): %d/%d = %.1f%%' % (nd, len(tab), 100.0 * nd / len(tab)))
        sys.stdout.flush()
