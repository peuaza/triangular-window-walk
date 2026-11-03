# -*- coding: utf-8 -*-
"""centre 族：二成分 Beta 混合拟合 + 矩比较"""
import numpy as np
from scipy import stats, optimize

def walk(a0, c=0.6):
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
                return ny / (n - 1)
            if nx < y0 or ny > x1:
                y0, x1 = y0 + 1, x1 - 1
                break
            y, x = ny, nx
    return None

n, h, N = 512, 256, 1000
rng = np.random.default_rng(31337)
xs = []
for _ in range(N):
    B = np.zeros((n, 2), dtype=np.int8)
    for c in range(2):
        B[rng.choice(n, h, replace=False), c] = 1
    d = walk(B.mean(1))
    if d is not None:
        xs.append(d)
x = np.array(xs)
xin = x[(x > 0) & (x < 1)]
print('N=%d ｜ mean %.4f ｜ sd %.4f ｜ skew %+.3f ｜ kurt %+.3f'
      % (len(x), x.mean(), x.std(ddof=1), stats.skew(x), stats.kurtosis(x)))

def nll_single(p):
    return -np.sum(stats.beta.logpdf(xin, p[0], p[1]))

def nll_mix(p):
    w, a1, b1, a2, b2 = p[0], p[1], p[2], p[3], p[4]
    if not (0 < w < 1) or min(a1, b1, a2, b2) <= 0.02:
        return 1e12
    d = w * stats.beta.pdf(xin, a1, b1) + (1 - w) * stats.beta.pdf(xin, a2, b2)
    return -np.sum(np.log(np.maximum(d, 1e-300)))

best = None
for w0 in (0.3, 0.5, 0.7):
    for a0, b0 in ((0.8, 1.5), (2, 3), (1.5, 2.5)):
        p0 = [w0, a0, b0, 2.0, 2.0]
        r = optimize.minimize(nll_mix, p0, method='Nelder-Mead',
                              options=dict(maxiter=8000, xatol=1e-5, fatol=1e-6))
        if best is None or r.fun < best.fun:
            best = r
p = best.x
print('二成分 Beta 混合：w=%.3f ｜ (a1,b1)=(%.2f,%.2f) ｜ (a2,b2)=(%.2f,%.2f) ｜ nll %.2f'
      % (p[0], p[1], p[2], p[3], p[4], best.fun))

cdf = lambda z: p[0] * stats.beta.cdf(z, p[1], p[2]) + (1 - p[0]) * stats.beta.cdf(z, p[3], p[4])
print('KS vs 二成分混合: %.4f ｜ 临界约 %.3f' % (stats.kstest(x, cdf).statistic, 1.36 / np.sqrt(len(x))))
h,edges = np.histogram(x, bins=20, range=(0,1))
print('二十箱直方图(%):', (h/h.sum()*100).round(1))
print('落在 [0.47,0.53] 的比例: %.3f' % np.mean((x>0.47)&(x<0.53)))
print('落在 [0.55,0.65] 的比例: %.3f' % np.mean((x>0.55)&(x<0.65)))
