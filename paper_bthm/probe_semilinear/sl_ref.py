# -*- coding: utf-8 -*-
"""sl_ref.py -- the literal reference engine for the threshold-run rule.

A line-by-line copy of rs_common.ref_walk with b1 replaced by the window scan.
Used only for validation (n <= 12).
"""
import numpy as np


def lrun2(z):
    best = cur = 0
    for v in z:
        if v == 2:
            cur += 1
            if cur > best:
                best = cur
        else:
            cur = 0
    return best


def ref_rule(a, n, t):
    a = np.asarray(a, dtype=np.int64)
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
        if steps > 40 * n * n:
            return None, steps
        if (ell, r) in seen:
            return None, steps
        seen.add((ell, r))
        z = a[ell:r + 1]
        m = len(z)
        k = max(1, m // 2)
        b2 = int(z.sum()) * k > int(z[:k].sum()) * m
        b1 = lrun2(z) > t
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
