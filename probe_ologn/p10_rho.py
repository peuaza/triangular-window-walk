# -*- coding: utf-8 -*-
"""p10: (a) K_global for small n, (b) rho = global longest 2-run vs window size at
b1-true, (c) step counts of the full walk and of the b1=0 walk."""
import sys, itertools
import numpy as np
sys.path.insert(0, '..')
from p0_tapes import enum_tapes
from _arr_longrun_fast import walk_naive, bit_b1_ref, longest_run_vec
from p2_b0 import walk_b0


def global_min_set(n):
    """smallest k such that some k-subset is pure for all tapes (brute force)."""
    A, D = [], []
    for a in enum_tapes(n):
        d, s, _ = walk_naive(a, n)
        A.append(a); D.append(-1 if d is None else int(round(d * (n - 1))))
    A = np.array(A); D = np.array(D); N = len(A)
    for k in range(1, n + 1):
        for S in itertools.combinations(range(n), k):
            key = A[:, S].astype(np.int64) @ (3 ** np.arange(k))
            ok = True
            for v in np.unique(key):
                if len(np.unique(D[key == v])) > 1:
                    ok = False; break
            if ok:
                return k
    return None


def rho_stats(n, nw, seed):
    rng = np.random.default_rng(seed)
    h = n // 2 if n % 2 == 0 else (n + 1) // 2
    rows = []
    for _ in range(nw):
        B = np.zeros((n, 2), dtype=np.int16)
        for c in range(2):
            B[rng.choice(n, h, replace=False), c] = 1
        a = B.sum(1).astype(np.int64)
        rho = longest_run_vec(a)
        # full walk, tracking m at b1 true
        L, R = 0, n - 1; ell, r = 0, n - 1; seen = set(); steps = 0
        tru = 0; maxm = 0; msum = 0
        while True:
            if ell >= r or R - L + 1 < 3 or steps > 20 * n * n or (ell, r) in seen:
                break
            steps += 1
            seen.add((ell, r))
            m = r - ell + 1
            b2 = int(np.sum(a[ell:r + 1]) * max(1, m // 2) > np.sum(a[ell:ell + max(1, m // 2)]) * m)
            b1 = bit_b1_ref(a[ell:r + 1])
            if b1:
                tru += 1; maxm = max(maxm, m); msum += m
            if b2 and b1: dy, dx = 1, 0
            elif b2: dy, dx = 0, -1
            elif b1: dy, dx = 0, 1
            else: dy, dx = -1, 0
            if ell + dy < L or r + dx > R:
                if ell == L and r == R: L, R = L + 1, R - 1
                elif ell + dy < L: L += 1
                else: R -= 1
                ell, r = L, R; seen.clear(); continue
            if ell + dy >= r + dx: break
            if r + dx < L or ell + dy > R:
                L, R = L + 1, R - 1; ell, r = L, R; seen.clear(); continue
            ell, r = ell + dy, r + dx
        steps0 = walk_b0(a, n)[1]
        rows.append((rho, tru, maxm, steps, steps0))
    R = np.array(rows, dtype=float)
    print('n=%d tapes=%d  rho: mean %.2f median %d max %d | b1true steps/tape %.2f | max m at b1true %.1f'
          % (n, nw, R[:, 0].mean(), np.median(R[:, 0]), R[:, 0].max(), R[:, 1].mean(), R[:, 2].max()))
    print('    violation of m < 4*rho: %d ;  steps/n^2 full %.4f  b0 %.4f'
          % (int((R[:, 2] >= 4 * R[:, 0]).sum()), R[:, 3].mean() / n ** 2, R[:, 4].mean() / n ** 2))


if __name__ == '__main__':
    for n in (6, 8):
        print('K_global(n=%d) =' % n, global_min_set(n))
    for n in (64, 128, 256, 512):
        rho_stats(n, 40, 2026)
