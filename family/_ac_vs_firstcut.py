# -*- coding: utf-8 -*-
"""首切方向与 ac 是否冗余：X ~ ac vs X ~ ac + 1[首切=列]"""
import numpy as np

def ac1(a):
    u, v = a[:-1], a[1:]
    if u.std() < 1e-12 or v.std() < 1e-12:
        return 0.0
    return float(np.corrcoef(u, v)[0, 1])

def walk(a0):
    a = a0
    n = len(a)
    y0, x1 = 0, n - 1
    first = None
    while x1 - y0 >= 3:
        y, x = y0, x1
        seen = set()
        while True:
            if (y, x) in seen:
                if first is None:
                    first = 'both'
                y0, x1 = y0 + 1, x1 - 1
                break
            seen.add((y, x))
            m = x - y + 1
            s = a[y:x + 1]
            k = max(1, m // 2)
            b2 = s.mean() > s[:k].mean()
            if m == 1:
                b1 = a[y] > 0.5
            elif m == 2:
                b1 = a[x] > a[y]
            else:
                u, v = s[:-1], s[1:]
                b1 = (np.corrcoef(u, v)[0, 1] > 0) if (u.std() > 1e-12 and v.std() > 1e-12) else False
            o = ("B" if b1 else "D") if b2 else ("C" if b1 else "A")
            ny, nx = {"A": (y - 1, x), "B": (y + 1, x), "C": (y, x + 1), "D": (y, x - 1)}[o]
            if ny < y0 or nx > x1:
                if first is None:
                    first = 'top' if ny < y0 else 'right'
                if ny < y0:
                    x1 = x - 1
                else:
                    y0 = y + 1
                break
            if ny >= nx:
                return ny / (n - 1), first if first else 'none'
            if nx < y0 or ny > x1:
                if first is None:
                    first = 'both'
                y0, x1 = y0 + 1, x1 - 1
                break
            y, x = ny, nx
    return None, None

n, h, NW = 256, 128, 800
rng = np.random.default_rng(2718)
X, A, T = [], [], []
for _ in range(NW):
    B = np.zeros((n, 2), dtype=np.int8)
    for c in range(2):
        B[rng.choice(n, h, replace=False), c] = 1
    a = B.mean(1)
    x, t = walk(a)
    if x is None:
        continue
    X.append(x); A.append(ac1(a)); T.append(t)

X = np.array(X); A = np.array(A)
col = np.array([1.0 if t == 'top' else 0.0 for t in T])
print('N=%d ｜ 首切=列 占 %.1f%% ｜ corr(ac, 1[列]) = %.3f'
      % (len(X), 100 * col.mean(), np.corrcoef(A, col)[0, 1]))

def fit(M, name):
    beta, *_ = np.linalg.lstsq(M, X, rcond=None)
    res = X - M @ beta
    print('%-24s R2 %.4f ｜ 残差 %.4f ｜ 系数 %s' % (name, 1 - res.var() / X.var(), res.std(ddof=M.shape[1]), np.round(beta[1:], 3)))

fit(np.column_stack([np.ones(len(X)), A]), 'ac')
fit(np.column_stack([np.ones(len(X)), col]), '1[首切=列]')
fit(np.column_stack([np.ones(len(X)), A, col]), 'ac + 1[首切=列]')
