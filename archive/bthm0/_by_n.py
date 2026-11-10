# -*- coding: utf-8 -*-
"""常数随 n 的依赖：同一套估计量，n = 128, 256, 512, 1024

对每个 n：生成 NW 条纸带，跑规范规则（BTHM-0，cap=20n），
再对每次相遇计算 ac1 与 20 个特征，做 5 折 CV：
  affine(ac1) / cubic(ac1) / GB(ac1) / GB(20 features)
"""
import sys

import numpy as np
from sklearn.ensemble import HistGradientBoostingRegressor
from sklearn.model_selection import KFold


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


def feats(a):
    n = len(a)
    f = {}
    u, v = a[:-1], a[1:]
    f['ac1'] = float(np.corrcoef(u, v)[0, 1]) if (u.std() > 1e-12 and v.std() > 1e-12) else 0.0
    f['mean'] = float(a.mean()); f['var'] = float(a.var())
    f['runs'] = float((np.abs(np.diff(a)) > 1e-9).sum()) / n
    idx = np.arange(n, dtype=float)
    f['centroid'] = float((a * idx).sum() / (a.sum() + 1e-12)) / n
    ones = a > 0.9
    f['first1'] = float(np.argmax(ones)) / n if ones.any() else 1.0
    for lag in (2, 3):
        uu, vv = a[:-lag], a[lag:]
        f['ac%d' % lag] = float(np.corrcoef(uu, vv)[0, 1]) if (uu.std() > 1e-12 and vv.std() > 1e-12) else 0.0
    for sc in (2, 4, 8, 16, 32, 64):
        bl = a[: (n // sc) * sc].reshape(-1, sc).mean(1)
        f['blk%d_mean' % sc] = float(bl.mean()); f['blk%d_sd' % sc] = float(bl.std())
    return f


print('%-6s %8s %8s %9s %9s %9s %9s' % ('n', '相遇', '率', 'affine', 'cubic', 'GB-ac1', 'GB-20'))
for n in (128, 256, 512, 1024):
    rng = np.random.default_rng(2026 + n)
    h = n // 2
    NW = 2000 if n <= 512 else 1200
    X, F = [], []
    for _ in range(NW):
        B = np.zeros((n, 2), dtype=np.int8)
        for c in range(2):
            B[rng.choice(n, h, replace=False), c] = 1
        a = B.mean(1)
        d = walk(a, n)
        if d is None:
            continue
        X.append(d); F.append(feats(a))
    X = np.array(X)
    keys = list(F[0].keys())
    M = np.array([[r[k] for k in keys] for r in F])
    i0 = keys.index('ac1')
    kf = KFold(5, shuffle=True, random_state=0)
    folds = list(kf.split(X))

    def cv(cols, deg=None, gb=False):
        sc = []
        for tr, te in folds:
            if gb:
                m = HistGradientBoostingRegressor(max_iter=250, learning_rate=.06, max_depth=5, random_state=0)
                m.fit(M[np.ix_(tr, cols)], X[tr]); p = m.predict(M[np.ix_(te, cols)])
            else:
                c = np.polyfit(M[tr, cols[0]], X[tr], deg); p = np.polyval(c, M[te, cols[0]])
            sc.append(1 - ((X[te] - p) ** 2).sum() / ((X[te] - X[tr].mean()) ** 2).sum())
        return float(np.mean(sc))

    a1 = cv([i0], deg=1)
    a3 = cv([i0], deg=3)
    g1 = cv([i0], gb=True)
    g20 = cv(list(range(len(keys))), gb=True)
    print('%-6d %8d %8.3f %9.4f %9.4f %9.4f %9.4f' % (n, len(X), len(X) / NW, a1, a3, g1, g20))
