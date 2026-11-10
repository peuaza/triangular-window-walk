# -*- coding: utf-8 -*-
"""p8: prefix-average form of b2 and the structure of the concavity set.

Lemma (verified here): with A_j = mean(a[L..L+j-1]), m = r-L+1,
    b2([L,r]) = 0  <=>  A_floor(m/2) >= A_m .
So the b1=0 walk scans m downwards and cuts left at the largest m <= M with
A_{floor(m/2)} >= A_m.
"""
import sys
import numpy as np
sys.path.insert(0, '..')
from _arr_longrun_fast import bit_b2_ref


def mktapes(n, nw, seed):
    rng = np.random.default_rng(seed)
    h = n // 2 if n % 2 == 0 else (n + 1) // 2
    for _ in range(nw):
        B = np.zeros((n, 2), dtype=np.int16)
        for c in range(2):
            B[rng.choice(n, h, replace=False), c] = 1
        yield B.sum(1).astype(np.int64)


def check_lemma():
    rng = np.random.default_rng(11)
    bad = 0
    tot = 0
    for _ in range(20000):
        n = int(rng.integers(2, 200))
        a = rng.integers(0, 3, size=n)
        L = int(rng.integers(0, n - 1))
        r = int(rng.integers(L + 1, n))
        m = r - L + 1
        k = max(1, m // 2)
        A = np.cumsum(a[L:r + 1]) / np.arange(1, m + 1)
        pred = int(A[k - 1] >= A[m - 1])
        tot += 1
        if pred != (1 - bit_b2_ref(a[L:r + 1])):
            bad += 1
    print('lemma check: %d random windows, mismatches %d' % (tot, bad))


def structure(n, nw, seed):
    tot_L = 0
    sizes = []
    nsign = []
    maxm = []
    for a in mktapes(n, nw, seed):
        for L in range(0, n - 3):
            M = n - L
            A = np.cumsum(a[L:]) / np.arange(1, M + 1)
            mm = np.arange(2, M + 1)
            B = A[mm - 1] - A[mm // 2 - 1]
            Z = mm[B <= 0]
            tot_L += 1
            sizes.append(len(Z))
            maxm.append((int(Z.max()) if len(Z) else -1) / M)
            # sign changes of B along increasing m
            nsign.append(int(np.sum(np.diff(np.sign(B)) != 0)))
    sizes = np.array(sizes)
    nsign = np.array(nsign)
    maxm = np.array(maxm, dtype=float)
    print('n=%d (L,M) pairs=%d:  |Z| mean %.1f median %.1f q90 %.1f max %d ; frac Z empty %.3f'
          % (n, tot_L, sizes.mean(), np.median(sizes), np.quantile(sizes, .9), sizes.max(),
             np.mean(sizes == 0)))
    print('    sign changes of B_m along m: mean %.2f median %d ; frac <=1 change %.3f'
          % (nsign.mean(), np.median(nsign), np.mean(nsign <= 1)))
    print('    max Z / M : mean %.2f  frac(maxZ<M) %.3f' % (np.mean(maxm), np.mean(maxm < 1)))


if __name__ == '__main__':
    check_lemma()
    for n in (128, 256, 512):
        structure(n, 3, 2026)
