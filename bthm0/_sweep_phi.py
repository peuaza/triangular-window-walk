# -*- coding: utf-8 -*-
"""策略空间扫描：phi（两个比特的定义）× 动作集

用法: python _sweep_phi.py <b2v> <b1v> <act> <n> <N> [seed]

b2v (趋势比特):
  0  mean(w) > mean(前半)          [规范]
  1  mean(w) > mean(后半)
  2  mean(w) > median(w)
  3  mean(w) > (min+max)/2
  4  median(前半) > median(后半)

b1v (纹理比特):
  0  ac1 > 0                        [规范]
  1  相邻递增比例 > 1/2
  2  ac2 > 0（滞后二阶）
  3  ac1 > 0 且 相邻递增比例 > 1/2
  4  median(w) > 1/2

act (动作集):
  0  四方向各一格                    [规范]
  1  四方向各两格
  2  滑动 (y+1, x-1)
  3  扩张方向一格的变体：一致时也走对角方向
"""
import sys

import numpy as np
from sklearn.ensemble import HistGradientBoostingRegressor
from sklearn.model_selection import KFold

B2V = int(sys.argv[1]) if len(sys.argv) > 1 else 0
B1V = int(sys.argv[2]) if len(sys.argv) > 2 else 0
ACT = int(sys.argv[3]) if len(sys.argv) > 3 else 0
N = int(sys.argv[4]) if len(sys.argv) > 4 else 256
NW = int(sys.argv[5]) if len(sys.argv) > 5 else 2000
SEED = int(sys.argv[6]) if len(sys.argv) > 6 else 2026
h = N // 2


def bits(z, a, ell, r):
    m = len(z)
    k = max(1, m // 2)
    if B2V == 0:
        b2 = z.mean() > z[:k].mean()
    elif B2V == 1:
        b2 = z.mean() > z[-k:].mean()
    elif B2V == 2:
        b2 = z.mean() > np.median(z)
    elif B2V == 3:
        b2 = z.mean() > 0.5 * (z.min() + z.max())
    else:
        b2 = np.median(z[:k]) > np.median(z[-k:])
    if B1V == 0:
        if m < 3:
            b1 = (z[-1] > z[0]) if m == 2 else (z[0] > 0.5)
        else:
            u, v = z[:-1], z[1:]
            b1 = (np.corrcoef(u, v)[0, 1] > 0) if (u.std() > 1e-12 and v.std() > 1e-12) else False
    elif B1V == 1:
        b1 = (np.diff(z) > 0).mean() > 0.5 if m >= 2 else (z[0] > 0.5)
    elif B1V == 2:
        if m < 4:
            b1 = (z[-1] > z[0]) if m == 2 else (z[0] > 0.5)
        else:
            u, v = z[:-2], z[2:]
            b1 = (np.corrcoef(u, v)[0, 1] > 0) if (u.std() > 1e-12 and v.std() > 1e-12) else False
    elif B1V == 3:
        if m < 3:
            b1 = (z[-1] > z[0]) if m == 2 else (z[0] > 0.5)
        else:
            u, v = z[:-1], z[1:]
            ac = (np.corrcoef(u, v)[0, 1] > 0) if (u.std() > 1e-12 and v.std() > 1e-12) else False
            b1 = ac and (np.diff(z) > 0).mean() > 0.5
    else:
        b1 = np.median(z) > 0.5
    return bool(b2), bool(b1)


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
        z = a[ell:r + 1]
        b2, b1 = bits(z, a, ell, r)
        if ACT == 2:
            ny, nx = ell + 1, r - 1
        else:
            st = 2 if ACT == 1 else 1
            ny, nx = 0, 0
            if b2 and b1:
                ny, nx = ell + st, r
            elif b2 and not b1:
                ny, nx = ell, r - st
            elif b1 and not b2:
                ny, nx = ell, r + st
            else:
                ny, nx = ell - st, r
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
    f['mean'] = float(a.mean())
    f['sd'] = float(a.std())
    f['runs'] = float(np.sum(np.diff(a) != 0))
    f['first1'] = float(np.argmax(a > 0.5)) if (a > 0.5).any() else float(n)
    for sc in (4, 8, 16, 32):
        bl = np.array([a[i:i + max(1, n // sc)].mean() for i in range(0, n, max(1, n // sc))])
        f['blk%d_mean' % sc] = float(bl.mean())
        f['blk%d_sd' % sc] = float(bl.std())
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

if len(X) < 30:
    print('phi=%d%d act=%d n=%d ｜ 纸带 %d ｜ 相遇 %d ｜ 率 %.3f ｜ 样本不足'
          % (B2V, B1V, ACT, N, NW, len(X), len(X) / NW))
    raise SystemExit

X = np.array(X)
keys = list(F[0].keys())
M = np.array([[r[k] for k in keys] for r in F])
kf = KFold(5, shuffle=True, random_state=SEED)
sc = []
for tr, te in kf.split(M):
    g = HistGradientBoostingRegressor(max_iter=80, learning_rate=.05, max_depth=3, min_samples_leaf=20, l2_regularization=1.0, random_state=0)
    g.fit(M[tr], X[tr])
    p = g.predict(M[te])
    sc.append(1 - ((X[te] - p) ** 2).sum() / ((X[te] - X[te].mean()) ** 2).sum())
print('phi=%d%d act=%d n=%d ｜ 纸带 %d ｜ 相遇 %d ｜ 率 %.3f ｜ CV R2 %.3f ± %.3f'
      % (B2V, B1V, ACT, N, NW, len(X), len(X) / NW, float(np.mean(sc)), float(np.std(sc) / np.sqrt(5))))

# 落盘：每一条纸带的特征与 delta，供之后任何预测器/统计复用
import os, json
os.makedirs('_data', exist_ok=True)
tagf = '_data/phi%d%d_act%d_n%d_N%d_s%d.jsonl' % (B2V, B1V, ACT, N, NW, SEED)
with open(tagf, 'w', encoding='utf-8') as fh:
    for i in range(len(X)):
        rec = dict(F[i])
        rec['delta'] = float(X[i])
        fh.write(json.dumps(rec) + '\n')
