# -*- coding: utf-8 -*-
"""§3 主表的逐估计量版本：affine / cubic / boosting / 21 特征

用法: python _table3.py <n> <N> [seed]
输出一行：n N seed meetings rate r2_affine r2_cubic r2_gb r2_gb21 r2_shuffle
并把（特征, delta）写入 _data/table3_n<n>_N<N>_s<seed>.jsonl
"""
import sys

import numpy as np
from sklearn.ensemble import HistGradientBoostingRegressor
from sklearn.linear_model import LinearRegression
from sklearn.model_selection import KFold
from sklearn.preprocessing import PolynomialFeatures

N = int(sys.argv[1]) if len(sys.argv) > 1 else 256
NW = int(sys.argv[2]) if len(sys.argv) > 2 else 5000
SEED = int(sys.argv[3]) if len(sys.argv) > 3 else 2026
h = N // 2


def walk(a, n):
    L, R = 0, n - 1
    ell, r = 0, n - 1
    seen = set()
    steps = 0
    while True:
        if ell >= r:
            return ell / (n - 1)
        if R - L + 1 < 3:
            return None
        steps += 1
        if steps > 40 * n:
            return None
        if (ell, r) in seen:
            return None
        seen.add((ell, r))
        m = r - ell + 1
        k = max(1, m // 2)
        z = a[ell:r + 1]
        b2 = z.mean() > z[:k].mean()
        if m == 1:
            b1 = z[0] > 0.5
        elif m == 2:
            b1 = z[1] > z[0]
        else:
            u, v = z[:-1], z[1:]
            b1 = (np.corrcoef(u, v)[0, 1] > 0) if (u.std() > 1e-12 and v.std() > 1e-12) else False
        if b2 and b1:
            ny, nx = ell + 1, r
        elif b2 and not b1:
            ny, nx = ell, r - 1
        elif b1 and not b2:
            ny, nx = ell, r + 1
        else:
            ny, nx = ell - 1, r
        if ny < L or nx > R:
            if ell == L and r == R:
                L, R = L + 1, R - 1
            elif ny < L:
                R = r - 1
            else:
                L = ell + 1
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
    f['ac2'] = float(np.corrcoef(a[:-2], a[2:])[0, 1]) if (a[:-2].std() > 1e-12 and a[2:].std() > 1e-12) else 0.0
    f['mean'] = float(a.mean())
    f['sd'] = float(a.std())
    f['runs'] = float(np.sum(np.diff(a) != 0))
    f['inc'] = float((np.diff(a) > 0).mean())
    f['first1'] = float(np.argmax(a > 0.5)) if (a > 0.5).any() else float(n)
    f['last1'] = float(n - 1 - np.argmax((a > 0.5)[::-1])) if (a > 0.5).any() else 0.0
    for sc in (4, 8, 16, 32):
        step = max(1, n // sc)
        bl = np.array([a[i:i + step].mean() for i in range(0, n, step)])
        f['blk%d_mean' % sc] = float(bl.mean())
        f['blk%d_sd' % sc] = float(bl.std())
        f['blk%d_trend' % sc] = float(np.polyfit(np.arange(len(bl)), bl, 1)[0]) if len(bl) > 1 else 0.0
    return f


rng = np.random.default_rng(SEED)
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

if len(X) < 50:
    print('n=%d N=%d seed=%d 相遇 %d 样本不足' % (N, NW, SEED, len(X)))
    raise SystemExit

X = np.array(X)
keys = list(F[0].keys())
M = np.array([[r[k] for k in keys] for r in F])
ac1 = M[:, keys.index('ac1')].reshape(-1, 1)
kf = list(KFold(5, shuffle=True, random_state=SEED).split(M))


def score(pred_fn, y):
    sc = []
    for tr, te in kf:
        p = pred_fn(tr, te, y)
        sc.append(1 - ((y[te] - p) ** 2).sum() / ((y[te] - y[te].mean()) ** 2).sum())
    return float(np.mean(sc))


def lin(deg):
    def fn(tr, te, y):
        if deg == 1:
            g = LinearRegression().fit(ac1[tr], y[tr])
            return g.predict(ac1[te])
        P = PolynomialFeatures(deg, include_bias=False)
        g = LinearRegression().fit(P.fit_transform(ac1[tr]), y[tr])
        return g.predict(P.transform(ac1[te]))
    return fn


def gb(cols):
    def fn(tr, te, y):
        g = HistGradientBoostingRegressor(max_iter=80, learning_rate=.05, max_depth=3,
                                          min_samples_leaf=20, l2_regularization=1.0, random_state=0)
        g.fit(M[tr][:, cols], y[tr])
        return g.predict(M[te][:, cols])
    return fn


allc = list(range(len(keys)))
print('n=%d N=%d seed=%d meetings=%d rate=%.3f affine=%.3f cubic=%.3f gb=%.3f gb21=%.3f shuffle=%.3f'
      % (N, NW, SEED, len(X), len(X) / NW,
         score(lin(1), X), score(lin(3), X), score(gb([keys.index('ac1')]), X), score(gb(allc), X),
         score(gb(allc), np.random.default_rng(SEED + 1).permutation(X))))

import os, json
os.makedirs('_data', exist_ok=True)
with open('_data/table3_n%d_N%d_s%d.jsonl' % (N, NW, SEED), 'w', encoding='utf-8') as fh:
    for i in range(len(X)):
        rec = {k: float(F[i][k]) for k in keys}
        rec['delta'] = float(X[i])
        fh.write(json.dumps(rec) + '\n')
