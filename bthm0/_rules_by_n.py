# -*- coding: utf-8 -*-
"""多规则 × 多 n：R² 随 n 的漂移是规范规则特有的，还是普遍的？

每条规则用同一个预测器（ac1 的 cubic 与 GB），5 折 CV。
"""
import sys

import numpy as np
from sklearn.ensemble import HistGradientBoostingRegressor
from sklearn.model_selection import KFold

N = int(sys.argv[1]) if len(sys.argv) > 1 else 256
NW = int(sys.argv[2]) if len(sys.argv) > 2 else 1500
RULES = ['base', 'jump2', 'centre:0.6', 'pos:0.5']


def walk(a, n, cfg):
    kind, _, arg = cfg.partition(':')
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
        if kind == 'pos':
            b2 = ((ell + r) / 2 / n) > float(arg)
        else:
            b2 = z.mean() > z[:k].mean()
        if m == 1:
            b1 = a[ell] > 0.5
        elif m == 2:
            b1 = a[r] > a[ell]
        else:
            u, v = z[:-1], z[1:]
            b1 = (np.corrcoef(u, v)[0, 1] > 0) if (u.std() > 1e-12 and v.std() > 1e-12) else False
        o = ("B" if b1 else "D") if b2 else ("C" if b1 else "A")
        step = 2 if kind == 'jump2' else 1
        ny, nx = {"A": (ell - step, r), "B": (ell + step, r),
                  "C": (ell, r + step), "D": (ell, r - step)}[o]
        if ny < L or nx > R:
            if kind == 'centre' and abs((L + R) / 2 / n - float(arg)) < 0.10:
                L, R = L + 1, R - 1
            elif ny < L:
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


print('%-12s %6s %8s %9s %9s %8s %8s' % ('rule', 'n', '相遇', 'cubic', 'GB-ac1', 'SE3', 'SEgb'))
for cfg in RULES:
    for n in (128, 256, 512):
        rng = np.random.default_rng(1000 + n)
        h = n // 2
        X, A = [], []
        for _ in range(NW):
            B = np.zeros((n, 2), dtype=np.int8)
            for c in range(2):
                B[rng.choice(n, h, replace=False), c] = 1
            a = B.mean(1)
            d = walk(a, n, cfg)
            if d is None:
                continue
            u, v = a[:-1], a[1:]
            ac = float(np.corrcoef(u, v)[0, 1]) if (u.std() > 1e-12 and v.std() > 1e-12) else 0.0
            X.append(d); A.append(ac)
        X = np.array(X); A = np.array(A).reshape(-1, 1)
        kf = KFold(5, shuffle=True, random_state=0)
        c3, gb = [], []
        for tr, te in kf.split(X):
            c = np.polyfit(A[tr, 0], X[tr], 3)
            c3.append(1 - ((X[te] - np.polyval(c, A[te, 0])) ** 2).sum() / ((X[te] - X[tr].mean()) ** 2).sum())
            m = HistGradientBoostingRegressor(max_iter=250, learning_rate=.06, max_depth=5, random_state=0)
            m.fit(A[tr], X[tr]); p = m.predict(A[te])
            gb.append(1 - ((X[te] - p) ** 2).sum() / ((X[te] - X[tr].mean()) ** 2).sum())
        se3 = np.std(c3, ddof=1) / np.sqrt(5)
        seg = np.std(gb, ddof=1) / np.sqrt(5)
        print('%-12s %6d %8d %9.4f %9.4f %8.4f %8.4f' % (cfg, n, len(X), np.mean(c3), np.mean(gb), se3, seg))
