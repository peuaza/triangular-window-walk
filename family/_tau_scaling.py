# -*- coding: utf-8 -*-
"""tau 族的 δ 随 n 的增长率：拟合 δ ~ n^alpha"""
import numpy as np

def walk(a0, tau):
    a = a0
    n = len(a)
    y0, x1 = 0, n - 1
    while x1 - y0 >= 3:
        y, x = y0, x1
        seen = set()
        while True:
            if (y, x) in seen:
                y0, x1 = y0 + 1, x1 - 1
                break
            seen.add((y, x))
            m = x - y + 1
            s = a[y:x + 1]
            b2 = s.mean() > s[:max(1, m // 2)].mean()
            if m == 1:
                b1 = a[y] > 0.5
            elif m == 2:
                b1 = a[x] > a[y]
            else:
                u, v = s[:-1], s[1:]
                b1 = (np.corrcoef(u, v)[0, 1] > tau) if (u.std() > 1e-12 and v.std() > 1e-12) else (tau < 0)
            o = ("B" if b1 else "D") if b2 else ("C" if b1 else "A")
            ny, nx = {"A": (y - 1, x), "B": (y + 1, x), "C": (y, x + 1), "D": (y, x - 1)}[o]
            if ny < y0:
                x1 = x - 1
                break
            if nx > x1:
                y0 = y + 1
                break
            if ny >= nx:
                return ny
            if nx < y0 or ny > x1:
                y0, x1 = y0 + 1, x1 - 1
                break
            y, x = ny, nx
    return None

print('%-6s %6s %8s %8s %8s' % ('n', 'N', 'mean δ', 'median δ', 'δ/n'))
rows = []
for n in (128, 256, 512, 1024):
    h = (n + 1) // 2
    rng = np.random.default_rng(4000 + n)
    xs = []
    NW = 500 if n <= 512 else 300
    for _ in range(NW):
        B = np.zeros((n, 2), dtype=np.int8)
        for c in range(2):
            B[rng.choice(n, h, replace=False), c] = 1
        d = walk(B.mean(1), 0.2)
        if d is not None:
            xs.append(d)
    x = np.array(xs)
    rows.append((n, len(x), x.mean(), np.median(x)))
    se = x.std(ddof=1) / np.sqrt(len(x))
    print('%-6d %6d %8.2f %8.2f %8.4f  SE %.2f'
          % (n, len(x), x.mean(), np.median(x), x.mean() / n, se), flush=True)

ns = np.array([r[0] for r in rows], float)
mus = np.array([r[2] for r in rows], float)
alpha, logc = np.polyfit(np.log(ns), np.log(mus), 1)
print()
print('拟合 mean δ = c·n^alpha ｜ alpha = %.3f ｜ c = %.3f' % (alpha, np.exp(logc)))
print('（alpha=0.5 表示扩散型增长，alpha=1 表示线性）')
