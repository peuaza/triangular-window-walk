# -*- coding: utf-8 -*-
"""相图：横轴 k（周期 n/k），纵轴 p（起点双切概率），颜色 = 峰度"""
import io, json
import numpy as np
from scipy import stats

def walk(a0, p, rng):
    a = a0
    n = len(a)
    y0, x1 = 0, n - 1
    while x1 - y0 >= 3:
        y, x = y0, x1
        first = True
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
                if first and rng.random() < p:
                    y0, x1 = y0 + 1, x1 - 1
                elif ny < y0:
                    x1 = x - 1
                else:
                    y0 = y + 1
                break
            if ny >= nx:
                return ny
            if nx < y0 or ny > x1:
                y0, x1 = y0 + 1, x1 - 1
                break
            y, x = ny, nx
            first = False
    return None

n, h, N = 256, 128, 300
KS = (1, 2, 4, 8)
PS = (0.0, 0.25, 0.5, 0.75, 1.0)
grid = {}
for k in KS:
    nb, hb = n // k, (n // k) // 2
    for p in PS:
        rng = np.random.default_rng(hash((k, int(p * 100))) % (2**31))
        xs = []
        for _ in range(N):
            base = np.zeros(nb)
            for _ in range(2):
                base[rng.choice(nb, hb, replace=False)] += 0.5
            a = np.tile(base, k)
            d = walk(a, p, rng)
            if d is not None:
                xs.append(d / (n - 1))
        if len(xs) < 25:
            grid[(k, p)] = dict(rate=len(xs) / N, sd=None, kurt=None, N=len(xs))
        else:
            x = np.array(xs)
            grid[(k, p)] = dict(rate=len(xs) / N, sd=round(float(x.std(ddof=1)), 3),
                                kurt=round(float(stats.kurtosis(x)), 2), N=len(xs))
        g = grid[(k, p)]
        print('k=%d p=%.2f ｜ rate %.3f ｜ sd %s ｜ kurt %s ｜ N %d'
              % (k, p, g['rate'], g['sd'], g['kurt'], g['N']), flush=True)

json.dump({('%d|%.2f' % kk): vv for kk, vv in grid.items()}, open('_phase_grid.json', 'w'), indent=1)
print('saved _phase_grid.json')
