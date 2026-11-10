# -*- coding: utf-8 -*-
"""Does the certificate statement survive if we ignore the 'no output' tapes?

delta is a PARTIAL function: the walk aborts (cycle, or region < 3) on a large
share of the ensemble.  lb_cert2 treated 'no output' as a value, so a set S can
fail to be a certificate merely because one completion is undefined.  Here we
re-test with the *defined* tapes only: S is a certificate in the restricted
sense if all defined tapes agreeing with a on S have the same delta.  Also
reports how often |group| = 1 (i.e. the counts alone pin the tape).
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


def key(a, S):
    return tuple(v for i, v in enumerate(a) if i not in S)


if __name__ == '__main__':
    for n in [int(x) for x in sys.argv[1:]] or [6, 8, 10]:
        tab = deltas(n)
        defd = {a: d for a, d in tab.items() if d is not None}
        print('n=%2d tapes=%6d defined=%d (%.1f%%)'
              % (n, len(tab), len(defd), 100.0 * len(defd) / len(tab)))
        ok_full, ok_res, sizes = 0, 0, []
        for S in itertools.combinations(range(n), n - 2):
            g = defaultdict(set)
            for a, d in tab.items():
                g[key(a, S)].add(d)
            if all(len(v) == 1 for v in g.values()):
                ok_full += 1
            gr = defaultdict(set)
            for a, d in defd.items():
                gr[key(a, S)].add(d)
            sizes.extend(len(v) for v in gr.values())
            if all(len(v) == 1 for v in gr.values()):
                ok_res += 1
        tot = len(list(itertools.combinations(range(n), n - 2)))
        print('      S of size n-2 (%d tested): %d certify the full ensemble,'
              ' %d certify the defined sub-ensemble; group sizes mean %.2f max %d'
              % (tot, ok_full, ok_res, np.mean(sizes), max(sizes)))
        sys.stdout.flush()
