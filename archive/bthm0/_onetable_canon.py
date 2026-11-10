# -*- coding: utf-8 -*-
"""规范规则版（出界双切）—— 由 _onetable.py 改写

一次运行，给出所有口径：样本内 R²、CV R²（1/3 阶）、GB 单特征、贪心嵌套

n=256, N 条纸带，固定步数上限 20n。输出一张可对账的表。
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
    seen = set()
    while True:
        if ell >= r:
            return ell / (n - 1)
        if R - L + 1 < 3:
            return None
        steps += 1
        if steps > 20 * n:
            return None
        if (ell, r) in seen:
            return None          # 成环 = 无输出
        seen.add((ell, r))
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
            if ell == L and r == R:          # 起点即出界 => 双切
                L, R = L + 1, R - 1
            elif ny < L:
                R = r - 1                    # 出底边 => 切列
            else:
                L = ell + 1                  # 出右边 => 切行
            ell, r = L, R
            seen.clear()
            continue
        if ny >= nx:
            return ny / (n - 1)
        if nx < L or ny > R:
            L, R = L + 1, R - 1
            ell, r = L, R
            seen.clear()
            continue
        ell, r = ny, nx


def feats(a):
    n = len(a)
    f = {}
    u, v = a[:-1], a[1:]
    f['ac1'] = float(np.corrcoef(u, v)[0, 1]) if (u.std() > 1e-12 and v.std() > 1e-12) else 0.0
    f['mean'] = float(a.mean())
    f['var'] = float(a.var())
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
        f['blk%d_mean' % sc] = float(bl.mean())
        f['blk%d_sd' % sc] = float(bl.std())
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
    X.append(d)
    F.append(feats(a))

X = np.array(X)
keys = list(F[0].keys())
M = np.array([[r[k] for k in keys] for r in F])
print('n=%d ｜ 纸带 %d ｜ 相遇 %d ｜ 率 %.3f' % (N, NW, len(X), len(X) / NW))

kf = KFold(5, shuffle=True, random_state=0)
folds = list(kf.split(X))


def fit_eval(cols, deg=None, gb=False):
    """返回 (样本内 R², CV R²)"""
    if gb:
        m = HistGradientBoostingRegressor(max_iter=250, learning_rate=.06, max_depth=5, random_state=0)
        m.fit(M[:, cols], X)
        ins = 1 - ((X - m.predict(M[:, cols])) ** 2).sum() / ((X - X.mean()) ** 2).sum()
    else:
        c = np.polyfit(M[:, cols[0]], X, deg)
        ins = 1 - ((X - np.polyval(c, M[:, cols[0]])) ** 2).sum() / ((X - X.mean()) ** 2).sum()
    cv = []
    for tr, te in folds:
        if gb:
            m = HistGradientBoostingRegressor(max_iter=250, learning_rate=.06, max_depth=5, random_state=0)
            m.fit(M[np.ix_(tr, cols)], X[tr]); p = m.predict(M[np.ix_(te, cols)])
        else:
            c = np.polyfit(M[tr, cols[0]], X[tr], deg); p = np.polyval(c, M[te, cols[0]])
        cv.append(1 - ((X[te] - p) ** 2).sum() / ((X[te] - X[tr].mean()) ** 2).sum())
    return float(ins), float(np.mean(cv)), float(np.std(cv, ddof=1) / np.sqrt(5))


print('%-34s %10s %10s %8s' % ('估计量', '样本内', 'CV', 'SE'))
i0 = keys.index('ac1')
for deg in (1, 3):
    a, b, c = fit_eval([i0], deg=deg)
    print('%-34s %10.4f %10.4f %8.4f' % ('ac1, 多项式 %d 阶' % deg, a, b, c))
a, b, c = fit_eval([i0], gb=True)
print('%-34s %10.4f %10.4f %8.4f' % ('ac1, 梯度提升', a, b, c))
a, b, c = fit_eval(list(range(len(keys))), gb=True)
print('%-34s %10.4f %10.4f %8.4f' % ('全部 20 特征, 梯度提升', a, b, c))
