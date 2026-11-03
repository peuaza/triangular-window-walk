# -*- coding: utf-8 -*-
"""centre 族的条件律：分别拟合 Beta / 正态，看能否显式写出混合"""
import numpy as np
from scipy import stats, optimize

def walk(a0, c=0.6):
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
            b2 = s.mean() > s[:max(1, m // 2)].mean()
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
                if abs((y0 + x1) / 2 / n - c) < 0.10:
                    y0, x1 = y0 + 1, x1 - 1
                elif ny < y0:
                    x1 = x - 1
                else:
                    y0 = y + 1
                break
            if ny >= nx:
                return first, ny / (n - 1)
            if nx < y0 or ny > x1:
                y0, x1 = y0 + 1, x1 - 1
                break
            y, x = ny, nx
    return first, None

n, h, N = 512, 256, 1500
rng = np.random.default_rng(1717)
G = {}
for _ in range(N):
    B = np.zeros((n, 2), dtype=np.int8)
    for c in range(2):
        B[rng.choice(n, h, replace=False), c] = 1
    f, d = walk(B.mean(1))
    if d is not None:
        G.setdefault(f if f else 'none', []).append(d)

for k in ('top', 'right', 'both'):
    if k not in G or len(G[k]) < 40:
        continue
    x = np.array(G[k])
    xin = x[(x > 0) & (x < 1)]
    nll = lambda p: -np.sum(stats.beta.logpdf(xin, p[0], p[1]))
    ab = optimize.minimize(nll, [2.0, 2.0], method='Nelder-Mead').x
    ks_b = stats.kstest(x, lambda z: stats.beta.cdf(z, ab[0], ab[1])).statistic
    ks_n = stats.kstest(x, lambda z: stats.norm.cdf(z, x.mean(), x.std(ddof=1))).statistic
    print('%-6s N=%4d ｜ mean %.3f sd %.3f ｜ Beta(%.2f,%.2f) KS %.3f ｜ Normal KS %.3f'
          % (k, len(x), x.mean(), x.std(ddof=1), ab[0], ab[1], ks_b, ks_n), flush=True)
