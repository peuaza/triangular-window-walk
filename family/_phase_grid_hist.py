# -*- coding: utf-8 -*-
"""相图（带直方图）：每格存 20 箱直方图"""
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

n, h, N = 256, 128, 400
KS = (1, 2, 4, 8)
PS = (0.0, 0.25, 0.5, 0.75, 1.0)
out = {}
for k in KS:
    nb = n // k
    for p in PS:
        rng = np.random.default_rng(1000 * k + int(p * 100) + 7)
        xs = []
        for _ in range(N):
            base = np.zeros(nb)
            for _ in range(2):
                base[rng.choice(nb, nb // 2, replace=False)] += 0.5
            a = np.tile(base, k)
            d = walk(a, p, rng)
            if d is not None:
                xs.append(d / (n - 1))
        rec = dict(rate=len(xs) / N, N=len(xs))
        if len(xs) >= 25:
            x = np.array(xs)
            hist, _ = np.histogram(x, bins=20, range=(0, 1))
            rec.update(sd=round(float(x.std(ddof=1)), 3),
                       kurt=round(float(stats.kurtosis(x)), 2),
                       hist=(hist / hist.sum()).round(3).tolist())
        out['%d|%.2f' % (k, p)] = rec
        print('k=%d p=%.2f ｜ %s' % (k, p, json.dumps(rec)[:120]), flush=True)

json.dump(out, open('_phase_grid_hist.json', 'w'))
print('saved _phase_grid_hist.json')
