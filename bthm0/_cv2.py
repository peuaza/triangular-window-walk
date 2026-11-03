# -*- coding: utf-8 -*-
"""标准测量：正则化提升树 + 打乱标签对照

用法: python _cv2.py <rule> <n> <N> [seed]
  rule: base | jump2 | slide | centre:0.6
输出一行：rule n tapes meetings rate R2 r2_shuffled sd
"""
import sys

import numpy as np
from sklearn.ensemble import HistGradientBoostingRegressor
from sklearn.model_selection import KFold

RULE = sys.argv[1] if len(sys.argv) > 1 else 'base'
N = int(sys.argv[2]) if len(sys.argv) > 2 else 256
NW = int(sys.argv[3]) if len(sys.argv) > 3 else 4000
SEED = int(sys.argv[4]) if len(sys.argv) > 4 else 2026
h = N // 2

KIND, _, ARG = RULE.partition(':')
THR = float(ARG) if ARG else 0.0


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
        if KIND == 'pos':
            b2 = (L + R) / 2 / n > THR
        if m == 1:
            b1 = z[0] > 0.5
        elif m == 2:
            b1 = z[1] > z[0]
        else:
            u, v = z[:-1], z[1:]
            b1 = (np.corrcoef(u, v)[0, 1] > 0) if (u.std() > 1e-12 and v.std() > 1e-12) else False
        st = 2 if KIND == 'jump2' else 1
        if KIND == 'slide':
            ny, nx = ell + 1, r - 1
        elif b2 and b1:
            ny, nx = ell + st, r
        elif b2 and not b1:
            ny, nx = ell, r - st
        elif b1 and not b2:
            ny, nx = ell, r + st
        else:
            ny, nx = ell - st, r
        if ny < L or nx > R:
            if KIND == 'centre' and abs((L + R) / 2 / n - THR) < 0.10:
                L, R = L + 1, R - 1
                ell, r = L, R
                seen.clear()
                continue
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
    print('%s n=%d N=%d seed=%d ｜ 相遇 %d ｜ 样本不足' % (RULE, N, NW, SEED, len(X)))
    raise SystemExit

X = np.array(X)
M = np.array([[r[k] for k in F[0]] for r in F])
kf = list(KFold(5, shuffle=True, random_state=SEED).split(M))


def cv(y):
    sc = []
    for tr, te in kf:
        g = HistGradientBoostingRegressor(max_iter=80, learning_rate=.05, max_depth=3,
                                          min_samples_leaf=20, l2_regularization=1.0, random_state=0)
        g.fit(M[tr], y[tr])
        p = g.predict(M[te])
        sc.append(1 - ((y[te] - p) ** 2).sum() / ((y[te] - y[te].mean()) ** 2).sum())
    return float(np.mean(sc)), float(np.std(sc) / np.sqrt(5))


r2, se = cv(X)
r2s, _ = cv(np.random.default_rng(SEED + 1).permutation(X))
print('%s n=%d N=%d seed=%d meetings=%d rate=%.3f R2=%.3f se=%.3f shuffle=%.3f'
      % (RULE, N, NW, SEED, len(X), len(X) / NW, r2, se, r2s))

import os, json
os.makedirs('_data', exist_ok=True)
tagf = '_data/cv2_%s_n%d_N%d_s%d.jsonl' % (RULE.replace(':', '_'), N, NW, SEED)
with open(tagf, 'w', encoding='utf-8') as fh:
    keys = list(F[0].keys())
    for i in range(len(X)):
        rec = {k: float(F[i][k]) for k in keys}
        rec['delta'] = float(X[i])
        fh.write(json.dumps(rec) + '\n')
