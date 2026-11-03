# -*- coding: utf-8 -*-
"""回文距离旋钮：从回文出发做 k 次随机换位，k↑ 越远离回文"""
import numpy as np
from scipy import stats

def walk(a):
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

n, h, N = 256, 128, 300
print('%-6s %7s %8s %8s  %s' % ('swaps', 'rate', 'sd', 'kurt', '二十箱直方图(%)'))
for k_sw in (0, 16, 64, 256, 1024, 4096):
    rng = np.random.default_rng(500 + k_sw % 97)
    xs = []
    for _ in range(N):
        b = np.zeros(n // 2)
        for _ in range(2):
            b[rng.choice(n // 2, n // 4, replace=False)] += 0.5
        a = np.concatenate([b, b[::-1]])          # 回文，满足 2n1+m=2h
        for _ in range(k_sw):                      # 随机换位（保持计数）
            i, j = rng.integers(0, n, 2)
            if a[i] != a[j]:
                a[i], a[j] = a[j], a[i]
        d = walk(a)
        if d is not None:
            xs.append(d / (n - 1))
    if len(xs) < 25:
        print('%-6d  样本不足 (%d)' % (k_sw, len(xs)))
        continue
    x = np.array(xs)
    hist, _ = np.histogram(x, bins=20, range=(0, 1))
    print('%-6d %7.3f %8.4f %8.3f  %s'
          % (k_sw, len(xs) / N, x.std(ddof=1), stats.kurtosis(x),
             (hist / hist.sum() * 100).round(0).astype(int)))
