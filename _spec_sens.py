# -*- coding: utf-8 -*-
"""规格敏感性（表 2）：四个变体，n=256

沿用 code/walk_bin.py 的精确几何与比特定义，只改切割/成环策略：
  A 出界切一侧、重访双切（论文主规则）
  B 出界切一侧、重访即停
  C 出界切双侧、重访双切
  D 出界切双侧、重访即停
"""
import sys

import numpy as np

N = int(sys.argv[1]) if len(sys.argv) > 1 else 256
NW = int(sys.argv[2]) if len(sys.argv) > 2 else 1500
h = N // 2


def walk(B, n, cut_both=False, revisit_stop=False):
    a = B.mean(1)
    y0, x1 = 0, n - 1
    while x1 - y0 >= 3:
        y, x = y0, x1
        seen = set()
        hit = None
        while True:
            if (y, x) in seen:
                hit = "cyc"
                break
            seen.add((y, x))
            m = x - y + 1
            b2 = a[y:x + 1].mean() > a[y:y + max(1, m // 2)].mean()
            if m == 1:
                b1 = a[y] > 0.5
            elif m == 2:
                b1 = a[x] > a[y]
            else:
                seg = a[y:x + 1]
                u, v = seg[:-1], seg[1:]
                b1 = False if (u.std() < 1e-12 or v.std() < 1e-12) else (np.corrcoef(u, v)[0, 1] > 0)
            o = ("B" if b1 else "D") if b2 else ("C" if b1 else "A")
            ny, nx = {"A": (y - 1, x), "B": (y + 1, x), "C": (y, x + 1), "D": (y, x - 1)}[o]
            if ny < y0:
                hit = "top"
                break
            if nx > x1:
                hit = "right"
                break
            if ny >= nx:
                return ny / (n - 1)
            if nx < y0 or ny > x1:
                hit = "bad"
                break
            y, x = ny, nx
        if hit == "cyc":
            if revisit_stop:
                return None
            y0, x1 = y0 + 1, x1 - 1
        elif hit in ("top", "right"):
            if cut_both:
                # 同时施加列切与行切：切掉一大块
                x1 = x - 1
                y0 = y + 1
            elif hit == "top":
                x1 = x - 1
            else:
                y0 = y + 1
        else:
            y0, x1 = y0 + 1, x1 - 1
    return None


def ks_arcsine(v):
    v = np.sort(v)
    n = len(v)
    F = 2 * np.arcsin(np.sqrt(np.clip(v, 0, 1))) / np.pi
    emp = np.arange(1, n + 1) / n
    return float(max(np.abs(emp - F).max(), np.abs(emp - 1 / n - F).max()))


rng = np.random.default_rng(2026)
print('n=%d ｜ %d 条/变体' % (N, NW))
print('%-36s %8s %8s %8s %8s' % ('variant', 'P(delta)', 'sd', 'median', 'KS'))
for label, cb, rs in (('exit one side, revisit shrinks', False, False),
                      ('exit one side, revisit stops', False, True),
                      ('exit both, revisit shrinks', True, False),
                      ('exit both, revisit stops', True, True)):
    vals = []
    for _ in range(NW):
        B = np.zeros((N, 2), dtype=np.int8)
        for c in range(2):
            B[rng.choice(N, h, replace=False), c] = 1
        d = walk(B, N, cb, rs)
        if d is not None:
            vals.append(d)
    v = np.array(vals)
    print('%-36s %8.3f %8.3f %8.3f %8.3f'
          % (label, len(v) / NW, v.std(ddof=1) if len(v) > 1 else 0,
             np.median(v) if len(v) else 0, ks_arcsine(v) if len(v) > 5 else float('nan')))
