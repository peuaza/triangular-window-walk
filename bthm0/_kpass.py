# -*- coding: utf-8 -*-
"""k 遍特征（k=1,2,3，全部 O(n) 读）能解释多少？

特征族（贪心加入，CV）：
  1 遍：ac1, mean, var, runs, centroid, first1
  2 遍：vol_ac = corr((a_i-m)^2, (a_{i+1}-m)^2)
        sign_ac = corr(sign(a_i-m), sign(a_{i+1}-m))
        cross_above = P(a_i>m 且 a_{i+1}>m) - P(a_i>m)P(a_{i+1}>m)
  3 遍：trend_ac = corr(r_i, r_{i+1})，r 为去线性趋势后的残差
        run_above = 在均值线上方的最长游程 / n
"""
import sys

import numpy as np
from sklearn.ensemble import HistGradientBoostingRegressor
from sklearn.model_selection import KFold

N = int(sys.argv[1]) if len(sys.argv) > 1 else 256
NW = int(sys.argv[2]) if len(sys.argv) > 2 else 4000
h = N // 2


def walk(a, n):
    L, R = 0, n - 1
    ell, r = 0, n - 1
    steps = 0
    while True:
        if ell >= r:
            return ell / (n - 1)
        if R - L + 1 < 3:
            return None
        steps += 1
        if steps > 20 * n:
            return None
        m = r - ell + 1
        k = max(1, m // 2)
        z = a[ell:r + 1]
        b2 = z.mean() > z[:k].mean()
        if m == 1:
            b1 = a[ell] > 0.5
        elif m == 2:
            b1 = a[r] > a[ell]
        else:
            u, v = z[:-1], z[1:]
            b1 = (np.corrcoef(u, v)[0, 1] > 0) if (u.std() > 1e-12 and v.std() > 1e-12) else False
        o = ("B" if b1 else "D") if b2 else ("C" if b1 else "A")
        ny, nx = {"A": (ell - 1, r), "B": (ell + 1, r), "C": (ell, r + 1), "D": (ell, r - 1)}[o]
        if ny < L or nx > R:
            if ny < L:
                R = r - 1
            else:
                L = ell + 1
            ell, r = L, R
            continue
        if ny >= nx:
            return ny / (n - 1)
        if nx < L or ny > R:
            L, R = L + 1, R - 1
            ell, r = L, R
            continue
        ell, r = ny, nx


def corr(u, v):
    return float(np.corrcoef(u, v)[0, 1]) if (u.std() > 1e-12 and v.std() > 1e-12) else 0.0


def feats(a):
    n = len(a)
    f = {}
    # ---- 1 遍 ----
    f['ac1'] = corr(a[:-1], a[1:])
    m = float(a.mean())
    f['mean'] = m
    f['var'] = float(a.var())
    f['runs'] = float((np.abs(np.diff(a)) > 1e-9).sum()) / n
    idx = np.arange(n, dtype=float)
    f['centroid'] = float((a * idx).sum() / (a.sum() + 1e-12)) / n
    ones = a > 0.9
    f['first1'] = float(np.argmax(ones)) / n if ones.any() else 1.0
    # ---- 2 遍（用到 m）----
    d = a - m
    f['vol_ac'] = corr(d[:-1] ** 2, d[1:] ** 2)
    sg = np.sign(d)
    f['sign_ac'] = corr(sg[:-1], sg[1:])
    ab = (a > m).astype(float)
    f['cross_above'] = float((ab[:-1] * ab[1:]).mean() - ab.mean() ** 2)
    # ---- 3 遍（去趋势后再相关）----
    t = np.arange(n, dtype=float)
    slope = float(np.polyfit(t, a, 1)[0]) if n > 2 else 0.0
    rr = a - m - slope * (t - t.mean())
    f['trend_ac'] = corr(rr[:-1], rr[1:])
    longest = cur = 0
    for v in ab:
        cur = cur + 1 if v > 0 else 0
        longest = max(longest, cur)
    f['run_above'] = longest / n
    return f


rng = np.random.default_rng(2026)
X, F = [], []
for _ in range(NW):
    B = np.zeros((N, 2), dtype=np.int8)
    for c in range(2):
        B[rng.choice(N, h, replace=False), c] = 1
    a = B.mean(1)
    d = walk(a, N)
    if d is None:
        continue
    X.append(d); F.append(feats(a))

X = np.array(X)
keys = list(F[0].keys())
M = np.array([[r[k] for k in keys] for r in F])
print('n=%d ｜ 相遇 %d ｜ 特征 %d 个' % (N, len(X), len(keys)))

kf = KFold(5, shuffle=True, random_state=0)
folds = list(kf.split(X))


def cv(cols):
    sc = []
    for tr, te in folds:
        mo = HistGradientBoostingRegressor(max_iter=250, learning_rate=.06, max_depth=5, random_state=0)
        mo.fit(M[np.ix_(tr, cols)], X[tr]); p = mo.predict(M[np.ix_(te, cols)])
        sc.append(1 - ((X[te] - p) ** 2).sum() / ((X[te] - X[tr].mean()) ** 2).sum())
    return float(np.mean(sc)), float(np.std(sc, ddof=1) / np.sqrt(5))


one = [i for i, k in enumerate(keys) if k in ('ac1', 'mean', 'var', 'runs', 'centroid', 'first1')]
two = one + [i for i, k in enumerate(keys) if k in ('vol_ac', 'sign_ac', 'cross_above')]
three = list(range(len(keys)))
for name, cols in (('1 pass (6 features)', one), ('2 passes (+3)', two), ('3 passes (+2)', three)):
    mu, se = cv(cols)
    print('%-20s %8.4f  ± %.4f' % (name, mu, se))
