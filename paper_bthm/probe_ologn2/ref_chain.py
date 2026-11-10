# -*- coding: utf-8 -*-
"""Referee check 2: exhaustive verification of the section-8 chain on small tapes.

For every tape of the ensemble (uniform on #2-#0 = 2*ceil(n/2)-n) at n=8..12:
  (G) L_i = i, R_i = n-1-delta_i, delta_{i+1}-delta_i in {0,1}, delta_i <= i
  (F) cost_i >= cost(i,n-1) - delta_i   at every realised epoch i
  (S) T >= sum_{i<N} cost(i,n-1) - N(N-1)/2   whenever the first N epochs are realised
Also: exact covariance identity of Lemma 2(i) at M=8 over the ensemble.
"""
import itertools, sys
import numpy as np


def prefix(a):
    P = np.zeros(len(a) + 1, dtype=np.int64)
    np.cumsum(a, out=P[1:])
    return P


def epoch(P, L, R):
    M = R - L + 1
    if M < 2:
        return None, 0
    m = np.arange(2, M + 1, dtype=np.int64)
    k = m // 2
    Sm = P[L + m] - P[L]
    Sk = P[L + k] - P[L]
    conc = Sk * m >= Sm * k
    idx = np.flatnonzero(conc)
    if idx.size == 0:
        return None, M - 1
    mstar = int(m[idx[-1]])
    return mstar, M - mstar + 1


def trace(a):
    n = len(a)
    P = prefix(a)
    L, R = 0, n - 1
    out = []
    while True:
        if L >= R or R - L + 1 < 3:
            break
        mstar, cost = epoch(P, L, R)
        out.append((L, R, R - L + 1, mstar, cost))
        if mstar is None:
            break
        if mstar == R - L + 1:
            L, R = L + 1, R - 1
        else:
            L += 1
    return out


def tapes_even(n):
    """all a in {0,1,2}^n with #2 = #0 (n even), all counts."""
    for k in range(n // 2 + 1):
        for pos2 in itertools.combinations(range(n), k):
            rest = [j for j in range(n) if j not in pos2]
            for pos0 in itertools.combinations(rest, k):
                a = np.ones(n, dtype=np.int64)
                a[list(pos2)] = 2
                a[list(pos0)] = 0
                yield a


badG = badF = badS = 0
ntapes = 0
maxviol = 0
for n in (8, 10, 12):
    cnt = 0
    for a in tapes_even(n):
        cnt += 1
        tr = trace(a)
        delta = 0
        T = sum(e[4] for e in tr)
        for i, (L, R, M, ms, cost) in enumerate(tr):
            if L != i:
                badG += 1
                print("GEOM L!=i", n, a)
            if R != n - 1 - delta:
                badG += 1
                print("GEOM R", n, a)
            if delta > i:
                badG += 1
                print("GEOM delta>i", n, a)
            csuf = epoch(prefix(a), i, n - 1)[1]
            if cost < csuf - delta:
                badF += 1
                print("FRONTIER FAIL n=%d i=%d a=%s cost=%d suf=%d delta=%d" %
                      (n, i, a, cost, csuf, delta))
            # advance delta for next epoch: double cut iff mstar == M
            if ms == M:
                delta += 1
        # sum chain
        P = prefix(a)
        real = len(tr)
        for N in range(1, min(real, n // 2) + 1):
            LB = sum(epoch(P, i, n - 1)[1] for i in range(N)) - N * (N - 1) // 2
            if T < LB:
                badS += 1
                print("SUM FAIL n=%d N=%d a=%s T=%d LB=%d" % (n, N, a, T, LB))
    print("n=%d  tapes=%d  geometry_bad=%d  frontier_bad=%d  sum_bad=%d" %
          (n, cnt, badG, badF, badS))
    ntapes += cnt

print("\ntotal tapes", ntapes, "G/F/S bad:", badG, badF, badS)

# ---------- exact covariance identity at M=8 over the ensemble ----------
print("\nLemma 2(i) check, M=8, exhaustive over #2=#0:")
rows = []
for a in tapes_even(8):
    rows.append(a - 1.0)
A = np.array(rows)
s2 = (A ** 2).mean(1).mean()
def D(a, th):
    M = len(a)
    i = int(round(th * M))
    j = int(round(th * M / 2))
    X = np.concatenate([[0], np.cumsum(a - 1.0)])
    return (X[i] - 2 * X[j]) / np.sqrt(s2 * M)
thetas = [0.25, 0.5, 0.75, 1.0]
Dm = np.array([[D(a, t) for t in thetas] for a in A])
C = Dm.T @ Dm / len(A)
pred = np.zeros((4, 4))
for p, t in enumerate(thetas):
    for q, u in enumerate(thetas):
        lo, hi = min(t, u), max(t, u)
        pred[p, q] = (2 * lo - hi if hi < 2 * lo else 0.0) * 8 / 7
print(" max |empirical - (2min-max)_+ * M/(M-1)| =", np.abs(C - pred).max())
print(" empirical diag:", np.round(np.diag(C), 4), " predicted:", np.round(np.diag(pred), 4))
