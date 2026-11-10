# -*- coding: utf-8 -*-
"""Certificate test, corrected: the projection must KEEP the read set S.

S is a certificate for tape a if every admissible tape agreeing with a on S has
the same delta.  Monotonicity (a superset of a certificate is a certificate)
means it suffices to test all S with |S| = n-2 in order to decide whether
K_global = n-1.  Two variants: 'delta' counts the 'no output' tapes as a value,
'defined' restricts to tapes where the walk halts.
"""
import itertools
import sys
from collections import defaultdict

import numpy as np

sys.path.insert(0, '.')
from p0_tapes import enum_tapes          # noqa: E402
from lb_cert2 import walk                # noqa: E402


def deltas(n):
    out = {}
    for a in enum_tapes(n):
        t = tuple(int(v) for v in a)
        d, _ = walk(np.array(t, dtype=np.int64), n)
        out[t] = None if d is None else round(d, 12)
    return out


def proj(a, S):
    return tuple(a[i] for i in S)


if __name__ == '__main__':
    for n in [int(x) for x in sys.argv[1:]] or [6, 8, 10]:
        tab = deltas(n)
        defd = {a: d for a, d in tab.items() if d is not None}
        ok_full, ok_res, gs = 0, 0, []
        for S in itertools.combinations(range(n), n - 2):
            g = defaultdict(set)
            for a, d in tab.items():
                g[proj(a, S)].add(d)
            gs.extend(len(v) for v in g.values())
            if all(len(v) == 1 for v in g.values()):
                ok_full += 1
            gr = defaultdict(set)
            for a, d in defd.items():
                gr[proj(a, S)].add(d)
            if all(len(v) == 1 for v in gr.values()):
                ok_res += 1
        tot = len(list(itertools.combinations(range(n), n - 2)))
        # sanity: n-1 cells must always suffice (counts force the missing cell)
        bad = 0
        for S in itertools.combinations(range(n), n - 1):
            g = defaultdict(set)
            for a, d in tab.items():
                g[proj(a, S)].add(d)
            bad += sum(1 for v in g.values() if len(v) > 1)
        print('n=%2d tapes=%6d defined=%.0f%% | S of size n-2: %d tested, %d certify all,'
              ' %d certify the defined part | mean #distinct deltas per projection %.2f'
              ' max %d | (n-1)-cell classes with >1 delta: %d  => K_global = %s'
              % (n, len(tab), 100.0 * len(defd) / len(tab), tot, ok_full, ok_res,
                 np.mean(gs), max(gs), bad, (n - 2) if ok_full else (n - 1)))
        sys.stdout.flush()
