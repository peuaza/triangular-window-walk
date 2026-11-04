# -*- coding: utf-8 -*-
"""Conjecture R at n = 14: rank of a maximal count-balanced block of the delta
matrix at the middle cut, computed by enumerating only the block.

A block is (left halves with profile p) x (right halves with profile q) where
p+q has #0-#2 = eps.  Every pair in it is a tape of the ensemble, so the block
sits inside the promise and CC_2way >= log2 rank(block).

Usage: python lb_rank14.py [n]
"""
import itertools
import sys
import time
from collections import defaultdict

import numpy as np

sys.path.insert(0, '.')
from lb_cert2 import walk          # noqa: E402  (same walk as p0/cert_n*.npy)


def perms(profile):
    """distinct words of length k over {0,1,2} with the given counts"""
    word = [0] * profile[0] + [1] * profile[1] + [2] * profile[2]
    return sorted(set(itertools.permutations(word)))


def blocks(n):
    k = n // 2
    eps = 0 if n % 2 == 0 else 1
    profs = []
    for c0 in range(k + 1):
        for c1 in range(k - c0 + 1):
            profs.append((c0, c1, k - c0 - c1))

    def multi(p):
        r = 1
        for c in p:
            r *= len(set(itertools.permutations(range(c)))) if False else 1
        from math import factorial
        return factorial(k) // (factorial(p[0]) * factorial(p[1]) * factorial(p[2]))

    out = []
    for p in profs:
        for q in profs:
            if p[0] + q[0] - p[2] - q[2] != eps:
                continue
            out.append((p, q, multi(p) * multi(q)))
    return k, sorted(out, key=lambda t: -t[2])


if __name__ == '__main__':
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 14
    k, bl = blocks(n)
    print('n=%d k=%d  candidate blocks by size: %s' % (n, k, [(p, q, s) for p, q, s in bl[:4]]))
    for (p, q, sz) in bl[:3]:
        t0 = time.time()
        L = perms(p)
        R = perms(q)
        d = np.full((len(L), len(R)), -1, dtype=np.int16)
        for i, lv in enumerate(L):
            for j, rv in enumerate(R):
                a = np.array(lv + rv, dtype=np.int64)
                v, _ = walk(a, n)
                d[i, j] = -1 if v is None else int(round(v * (n - 1)))
        r = int(np.linalg.matrix_rank(d.astype(np.float64)))
        md = min(d.shape)
        print('  block %s x %s  shape %s  rank %d  rank/min-dim %.3f  log2 rank %.2f bits  (%.1fs)'
              % (p, q, d.shape, r, r / md, np.log2(r) if r else 0, time.time() - t0))
        sys.stdout.flush()
