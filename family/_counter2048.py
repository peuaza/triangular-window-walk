# -*- coding: utf-8 -*-
"""n=2048：counter2 的 Beta 拟合（后台跑，带进度）"""
import json
import time
import numpy as np
from scipy import stats, optimize

def walk(a0):
    a = a0
    n = len(a)
    y0, x1 = 0, n - 1
    cuts = 0
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
                b1 = (np.corrcoef(u, v)[0, 1] > 0) if (u.std() > 1e-12 and v.std() > 1e-12) else False
            o = ("B" if b1 else "D") if b2 else ("C" if b1 else "A")
            ny, nx = {"A": (y - 1, x), "B": (y + 1, x), "C": (y, x + 1), "D": (y, x - 1)}[o]
            if ny < y0 or nx > x1:
                if cuts % 2 == 0:
                    y0, x1 = y0 + 1, x1 - 1
                elif ny < y0:
                    x1 = x - 1
                else:
                    y0 = y + 1
                cuts += 1
                break
            if ny >= nx:
                return ny
            if nx < y0 or ny > x1:
                y0, x1 = y0 + 1, x1 - 1
                break
            y, x = ny, nx
    return None

n, h, NW = 2048, 1024, 600
rng = np.random.default_rng(2048)
xs = []
t0 = time.time()
for i in range(NW):
    B = np.zeros((n, 2), dtype=np.int8)
    for c in range(2):
        B[rng.choice(n, h, replace=False), c] = 1
    d = walk(B.mean(1))
    if d is not None:
        xs.append(d / (n - 1))
    if (i + 1) % 50 == 0:
        print('  %d/%d walks ｜ %d δ ｜ %.0fs' % (i + 1, NW, len(xs), time.time() - t0), flush=True)

x = np.array(xs)
xin = x[(x > 0) & (x < 1)]
nll = lambda p: -np.sum(stats.beta.logpdf(xin, p[0], p[1]))
ab = optimize.minimize(nll, [2.5, 2.5], method='Nelder-Mead').x
ks = stats.kstest(x, lambda z: stats.beta.cdf(z, ab[0], ab[1])).statistic
print('n=2048 ｜ N=%d ｜ a=%.3f b=%.3f ｜ KS=%.4f ｜ sd %.4f κ %+.3f ｜ 用时 %.0fs'
      % (len(x), ab[0], ab[1], ks, x.std(ddof=1), stats.kurtosis(x), time.time() - t0))
json.dump(dict(n=n, N=len(x), a=float(ab[0]), b=float(ab[1]), ks=float(ks),
               sd=float(x.std(ddof=1)), kurt=float(stats.kurtosis(x))),
          open('_counter2048.json', 'w'), indent=1)
