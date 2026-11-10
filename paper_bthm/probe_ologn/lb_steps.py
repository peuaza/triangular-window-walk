# -*- coding: utf-8 -*-
"""How long are the four rules' trajectories?  steps/n^2 must stay bounded away
from 0 for the O(n polylog n) question to be non-trivial (O(1) reads/step would
then already give O(n^2) -- and if steps were O(n polylog) the claim would be
REFUTED outright)."""
import sys
import numpy as np

RULES = ['longrun', 'runs', 'pattern', 'pal']
PAT = np.array([2, 1, 1])


def longest_run(z):
    b = (z > 1).astype(np.int8)
    if not b.any():
        return 0
    d = np.diff(np.concatenate((np.zeros(1, np.int8), b, np.zeros(1, np.int8))))
    st = np.flatnonzero(d == 1)
    en = np.flatnonzero(d == -1)
    return int((en - st).max())


def bit_b1(z, rule):
    m = len(z)
    if rule == 'longrun':
        return longest_run(z) * 4 > m
    if rule == 'runs':
        return int((np.diff(z) != 0).sum()) * 2 > m
    if rule == 'pattern':
        if m < 3:
            return False
        return bool((np.lib.stride_tricks.sliding_window_view(z, 3) == PAT).all(1).any())
    if rule == 'pal':
        return bool((z == z[::-1]).all())
    raise ValueError


def walk(a, n, rule):
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
        if steps > 4 * n * n:
            return None, steps
        if (ell, r) in seen:
            return None, steps
        seen.add((ell, r))
        z = a[ell:r + 1]
        m = r - ell + 1
        k = max(1, m // 2)
        b2 = int(z.sum()) * k > int(z[:k].sum()) * m
        b1 = bit_b1(z, rule)
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


for N in [int(x) for x in sys.argv[1:]] or [64, 128, 256]:
    NW = max(40, 4000 // N)
    H = N // 2
    print('n=%d tapes=%d' % (N, NW))
    for rule in RULES:
        rng = np.random.default_rng(2026)
        st, ok = [], 0
        for _ in range(NW):
            B = np.zeros((N, 2), dtype=np.int16)
            for c in range(2):
                B[rng.choice(N, H, replace=False), c] = 1
            a = B.sum(1).astype(np.int64)
            d, s = walk(a, N, rule)
            st.append(s)
            ok += (d is not None)
        st = np.array(st, float)
        print('  %-8s meet=%-4d steps=%8.0f  steps/n^2=%.4f  steps/(n logn)=%.2f  max/n^2=%.4f'
              % (rule, ok, st.mean(), st.mean() / N**2,
                 st.mean() / (N * np.log(N)), st.max() / N**2))
    sys.stdout.flush()
