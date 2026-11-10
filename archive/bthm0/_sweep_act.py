# -*- coding: utf-8 -*-
"""动作空间扫描：切割与停机作为显式动作，以及依赖位置的 phi

用法: python _sweep_act.py <posbit> <cutλ> <haltλ> <n> <N> [seed]

posbit（位置比特，替换 b1）:
  0  ac1 > 0（基准，不看位置）
  1  [ℓ-L > R-r]   左侧剩余比右侧多
  2  [中心 > 1/2]
  3  [m > n/2]

cutλ（主动切割阈值）:
  0    不主动切（只在出界时切）
  λ>0  当 m > λ·n 时，同时在两侧各切一格（作为一步动作，不移动）

haltλ（主动停机阈值）:
  0    不停机
  λ>0  当 m < λ·n 时停机、无输出

b2 固定为规范定义（均值 vs 前半）；动作固定为四方向各一格。
"""
import sys

import numpy as np
from sklearn.ensemble import HistGradientBoostingRegressor
from sklearn.model_selection import KFold

POSBIT = int(sys.argv[1]) if len(sys.argv) > 1 else 0
CUTL = float(sys.argv[2]) if len(sys.argv) > 2 else 0.0
HALTL = float(sys.argv[3]) if len(sys.argv) > 3 else 0.0
N = int(sys.argv[4]) if len(sys.argv) > 4 else 256
NW = int(sys.argv[5]) if len(sys.argv) > 5 else 2000
SEED = int(sys.argv[6]) if len(sys.argv) > 6 else 2026
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
        m = r - ell + 1
        if HALTL > 0 and m < HALTL * (R - L + 1):
            return None
        steps += 1
        if steps > 40 * n:
            return None
        if (ell, r) in seen:
            return None
        seen.add((ell, r))
        # 主动切割：窗口占区域的比例过高时先切（相对尺度）
        if CUTL > 0 and (R - L + 1) < n and m < CUTL * (R - L + 1):
            L, R = L + 1, R - 1
            ell, r = L, R
            seen.clear()
            continue
        z = a[ell:r + 1]
        k = max(1, m // 2)
        b2 = z.mean() > z[:k].mean()
        if POSBIT == 0:
            if m < 3:
                b1 = (z[-1] > z[0]) if m == 2 else (z[0] > 0.5)
            else:
                u, v = z[:-1], z[1:]
                b1 = (np.corrcoef(u, v)[0, 1] > 0) if (u.std() > 1e-12 and v.std() > 1e-12) else False
        elif POSBIT == 1:
            b1 = (ell - L) > (R - r)
        elif POSBIT == 2:
            b1 = (L + R) / 2 / n > 0.5
        else:
            b1 = m > n / 2
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

tag = 'pos=%d cut=%.2f halt=%.2f n=%d' % (POSBIT, CUTL, HALTL, N)
if len(X) < 30:
    print('%s ｜ 纸带 %d ｜ 相遇 %d ｜ 率 %.3f ｜ 样本不足' % (tag, NW, len(X), len(X) / NW))
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
print('%s ｜ 纸带 %d ｜ 相遇 %d ｜ 率 %.3f ｜ CV R2 %.3f ± %.3f'
      % (tag, NW, len(X), len(X) / NW, float(np.mean(sc)), float(np.std(sc) / np.sqrt(5))))

import os, json
os.makedirs('_data', exist_ok=True)
tagf = '_data/act_p%d_c%.2f_h%.2f_n%d_N%d_s%d.jsonl' % (POSBIT, CUTL, HALTL, N, NW, SEED)
with open(tagf, 'w', encoding='utf-8') as fh:
    for i in range(len(X)):
        rec = dict(F[i])
        rec['delta'] = float(X[i])
        fh.write(json.dumps(rec) + '\n')
