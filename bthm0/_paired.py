# -*- coding: utf-8 -*-
"""配对差：同一折上 R²(ac1) − R²(20特征)，给出均值、SE、分位

若差的 CI 覆盖 0 且宽度小 ⇒ 「加特征无用」是等价性陈述，而非「未拒绝」。
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
i0 = keys.index('ac1')
print('n=%d ｜ 相遇 %d' % (N, len(X)))

kf = KFold(5, shuffle=True, random_state=0)
diffs = []
for tr, te in kf.split(X):
    scores = []
    for cols in ([i0], list(range(len(keys)))):
        m = HistGradientBoostingRegressor(max_iter=250, learning_rate=.06, max_depth=5, random_state=0)
        m.fit(M[np.ix_(tr, cols)], X[tr])
        p = m.predict(M[np.ix_(te, cols)])
        scores.append(1 - ((X[te] - p) ** 2).sum() / ((X[te] - X[tr].mean()) ** 2).sum())
    diffs.append(scores[0] - scores[1])   # ac1 − 20特征

d = np.array(diffs)
print('逐折差（ac1 减 20 特征）：%s' % ' '.join('%.4f' % v for v in d))
print('均值 %.4f ｜ SE %.4f ｜ 95%% CI [%.4f, %.4f]' % (d.mean(), d.std(ddof=1) / np.sqrt(5),
                                                    d.mean() - 2.776 * d.std(ddof=1) / np.sqrt(5),
                                                    d.mean() + 2.776 * d.std(ddof=1) / np.sqrt(5)))
