# -*- coding: utf-8 -*-
"""记忆层级 M0–M3 的密度对比图"""
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import stats

def walk(a0, mode, par=1):
    a = a0
    n = len(a)
    y0, x1 = 0, n - 1
    cuts = 0
    hist = []
    while x1 - y0 >= 3:
        y, x = y0, x1
        seen = set()
        while True:
            if (y, x) in seen:
                y0, x1 = y0 + 1, x1 - 1
                cuts += 1; hist.append('B')
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
                sd_ = 'T' if ny < y0 else 'R'
                if mode == 'M1' and y0 > 0.5 * n:
                    y0, x1 = y0 + 1, x1 - 1
                elif mode == 'M2' and cuts % 2 == 0:
                    y0, x1 = y0 + 1, x1 - 1
                elif mode == 'M3' and len(hist) >= 2 and hist[-1] == hist[-2] == sd_:
                    if sd_ == 'T':
                        y0 = y + 1
                    else:
                        x1 = x - 1
                elif sd_ == 'T':
                    x1 = x - 1
                else:
                    y0 = y + 1
                hist.append(sd_); cuts += 1
                break
            if ny >= nx:
                return ny
            if nx < y0 or ny > x1:
                y0, x1 = y0 + 1, x1 - 1
                break
            y, x = ny, nx
    return None

n, h, N = 256, 128, 400
rng0 = np.random.default_rng(8888)
BASE = []
for _ in range(N):
    B = np.zeros((n, 2), dtype=np.int8)
    for c in range(2):
        B[rng0.choice(n, h, replace=False), c] = 1
    BASE.append(B.mean(1))

out = {}
for mode, label in (('M0', 'M0  no memory'), ('M1', 'M1  1 bit (position)'),
                    ('M2', 'M2  counter (cuts mod 2)'), ('M3', 'M3  history (last 2 sides)')):
    xs = []
    for a in BASE:
        d = walk(a, mode)
        if d is not None:
            xs.append(d / (n - 1))
    x = np.array(xs)
    hist, _ = np.histogram(x, bins=20, range=(0, 1))
    out[mode] = dict(label=label, sd=float(x.std(ddof=1)), kurt=float(stats.kurtosis(x)),
                     rate=len(xs) / N, hist=(hist / hist.sum()).round(3).tolist())
    print('%s ｜ sd %.3f κ %+.2f rate %.2f ｜ %s' % (label, out[mode]['sd'], out[mode]['kurt'],
                                                    out[mode]['rate'], (hist / hist.sum() * 100).round(0).astype(int)), flush=True)
json.dump(out, open('_memory_levels.json', 'w'), indent=1)

edges = np.linspace(0, 1, 21)
fig, axes = plt.subplots(1, 4, figsize=(13, 3.6), dpi=150, sharey=True)
for ax, (mode, rec) in zip(axes, out.items()):
    hh = np.array(rec['hist'])
    ax.fill_between(edges[:-1], 0, hh, step='post', color='#4a7fb5', alpha=0.85, lw=0)
    ax.plot(edges[:-1], hh, drawstyle='steps-post', color='black', lw=1.1)
    ax.set_title('%s\nsd %.3f, κ %+.2f' % (rec['label'], rec['sd'], rec['kurt']), fontsize=10)
    ax.set_xlim(0, 1); ax.set_xticks([0, .5, 1]); ax.set_yticks([])
axes[0].set_ylabel('density', fontsize=10)
fig.suptitle('Memory hierarchy: the shape changes only when the walk remembers something'
             ' (kappa = excess kurtosis; n=256, 400 matrices per panel)', fontsize=12)
fig.tight_layout()
fig.savefig('figs/memory_hierarchy_v2.png', bbox_inches='tight')
print('saved figs/memory_hierarchy.png')
