# -*- coding: utf-8 -*-
"""按首次切割类型拆分计数：N(n,k) = N_top + N_right + N_both"""
import itertools
import json
import numpy as np
from math import comb

def walk(a):
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
            if ny < y0:
                if first is None:
                    first = 'top'
                x1 = x - 1
                break
            if nx > x1:
                if first is None:
                    first = 'right'
                y0 = y + 1
                break
            if ny >= nx:
                return first, ny
            if nx < y0 or ny > x1:
                if first is None:
                    first = 'both'
                y0, x1 = y0 + 1, x1 - 1
                break
            y, x = ny, nx
    return first, None

out = {}
for n in (8, 9, 10, 11, 12):
    h = (n + 1) // 2
    acc = {t: [0] * n for t in ('top', 'right', 'both', 'none')}
    for tup in itertools.product((0.0, 0.5, 1.0), repeat=n):
        m = sum(1 for v in tup if v == 0.5)
        n1 = sum(1 for v in tup if v == 1.0)
        if 2 * n1 + m != 2 * h:
            continue
        f, d = walk(np.array(tup))
        if d is None:
            continue
        acc[f if f else 'none'][d] += comb(m, m // 2)
    out[n] = {t: v for t, v in acc.items()}
    tot = {t: sum(v) for t, v in acc.items()}
    print('n=%2d ｜ 按首切的矩阵计数：%s' % (n, tot), flush=True)
    for t in ('top', 'right', 'both'):
        if sum(acc[t]):
            print('      %-6s k=1..n-1: %s' % (t, [v for v in acc[t][1:]]))

json.dump({str(k): v for k, v in out.items()}, open('_counts_by_first.json', 'w'))
print('saved _counts_by_first.json')
