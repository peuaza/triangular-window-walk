# -*- coding: utf-8 -*-
"""游走型规则：内点切割

与规范规则相同的比特与移动表，只改切割：
当窗口内出现「连续三个 1」时，切在该段第一个 1 的位置，保留左侧、丢弃右侧。
新窗口不是旧窗口的删减，累加器必须重建，代价 = 重建窗口长度。

用法:
  python _walk_rule.py bench <n> <tapes> [seed]   读计数
  python _walk_rule.py cv    <n> <tapes> [seed]   特征 + CV + 打乱对照 + 落盘
"""
import json
import os
import sys

import numpy as np
from sklearn.ensemble import HistGradientBoostingRegressor
from sklearn.model_selection import KFold

MODE = sys.argv[1] if len(sys.argv) > 1 else 'bench'
N = int(sys.argv[2]) if len(sys.argv) > 2 else 256
NW = int(sys.argv[3]) if len(sys.argv) > 3 else 500
SEED = int(sys.argv[4]) if len(sys.argv) > 4 else 2026
h = N // 2


def walk(a, n, counter=None):
    """规范比特与移动 + 内点切割。counter 是单元素列表，累计读次数。"""
    reads = 0

    def rd(k=1):
        nonlocal reads
        reads += k
        return reads

    L, R = 0, n - 1
    ell, r = 0, n - 1
    seen = set()
    steps = 0
    # 初始窗口统计
    rd(R - L + 1)
    while True:
        if ell >= r:
            return ell / (n - 1), reads
        if R - L + 1 < 3:
            return None, reads
        steps += 1
        if steps > 200 * n:
            return None, reads
        if (ell, r) in seen:
            return None, reads
        seen.add((ell, r))

        z = a[ell:r + 1]
        rd(len(z))
        m = len(z)
        k = max(1, m // 2)
        b2 = z.mean() > z[:k].mean()
        if m == 1:
            b1 = z[0] > 0.5
        elif m == 2:
            b1 = z[1] > z[0]
        else:
            u, v = z[:-1], z[1:]
            b1 = (np.corrcoef(u, v)[0, 1] > 0) if (u.std() > 1e-12 and v.std() > 1e-12) else False

        # 内点切割条件：窗口超过 n/2 且后半段有 1，切在那个 1 处
        cand = -1
        if m > n // 2:
            half = m // 2
            seg = z[half:]
            pos = np.argmax(seg > 0.5)
            if seg[pos] > 0.5 and half + pos > 0:
                cand = half + pos
        if cand > 0:
            # 切在 cand 处：新窗口 = [ell, ell+cand-1]
            newR = ell + cand - 1
            if newR - L + 1 >= 3:
                R = newR
                ell, r = L, R
                rd(R - L + 1)      # 重建
                seen.clear()
                continue
            return None, reads

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
            rd(2)
            continue
        if ny >= nx:
            return ny / (n - 1), reads
        if nx < L or ny > R:
            L, R = L + 1, R - 1
            ell, r = L, R
            seen.clear()
            rd(2)
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
    f['maxrun'] = float(max((len(s) for s in ''.join('1' if x > 0.5 else '0' for x in a).split('0')), default=0))
    for sc in (4, 8, 16, 32):
        step = max(1, n // sc)
        bl = np.array([a[i:i + step].mean() for i in range(0, n, step)])
        f['blk%d_mean' % sc] = float(bl.mean())
        f['blk%d_sd' % sc] = float(bl.std())
        f['blk%d_trend' % sc] = float(np.polyfit(np.arange(len(bl)), bl, 1)[0]) if len(bl) > 1 else 0.0
    return f


rng = np.random.default_rng(SEED)
X, F, RD = [], [], []
for _ in range(NW):
    B = np.zeros((N, 2), dtype=np.int8)
    for c in range(2):
        B[rng.choice(N, h, replace=False), c] = 1
    a = B.mean(1)
    d, rr = walk(a, N)
    RD.append(rr)
    if d is None:
        continue
    X.append(d)
    F.append(feats(a))

rd = np.array(RD)
print('MODE=%s n=%d tapes=%d seed=%d meet=%d rate=%.3f reads_mean=%.0f reads/n=%.2f reads/n^2=%.4f'
      % (MODE, N, NW, SEED, len(X), len(X) / NW, rd.mean(), rd.mean() / N, rd.mean() / N ** 2))

if MODE != 'cv' or len(X) < 50:
    raise SystemExit

X = np.array(X)
keys = list(F[0].keys())
M = np.array([[r[k] for k in keys] for r in F])
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
print('CV R2=%.3f se=%.3f shuffle=%.3f' % (r2, se, r2s))
os.makedirs('_data', exist_ok=True)
with open('_data/walk2_n%d_N%d_s%d.jsonl' % (N, NW, SEED), 'w', encoding='utf-8') as fh:
    for i in range(len(X)):
        rec = {k: float(F[i][k]) for k in keys}
        rec['delta'] = float(X[i])
        fh.write(json.dumps(rec) + '\n')
