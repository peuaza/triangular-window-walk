# -*- coding: utf-8 -*-
"""新动作：跳跃（B/D 跨 2 或 3 格）、滑动（第三比特选 (y+1,x-1)）、局部改写"""
import numpy as np
from scipy import stats

def walk(a0, mode):
    a = a0.copy()
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

            step = 2 if mode == 'jump2' else (3 if mode == 'jump3' else 1)
            if mode == 'slide':
                b3 = np.median(s) > 0.5
                if b3:                                  # 滑动：保持窗口长度，向右平移
                    ny, nx = y + 1, x - 1
                else:
                    o = ("B" if b1 else "D") if b2 else ("C" if b1 else "A")
                    ny, nx = {"A": (y - 1, x), "B": (y + 1, x), "C": (y, x + 1), "D": (y, x - 1)}[o]
            else:
                o = ("B" if b1 else "D") if b2 else ("C" if b1 else "A")
                ny, nx = {"A": (y - 1, x), "B": (y + step, x), "C": (y, x + 1), "D": (y, x - step)}[o]

            if mode == 'edit':                          # 局部改写：翻转当前位置的格子
                a[y] = 1.0 - a[y] if a[y] != 0.5 else 0.5

            if ny < y0:
                x1 = x - 1
                break
            if nx > x1:
                y0 = y + 1
                break
            if ny >= nx:
                return min(ny, n - 1)
            if nx < y0 or ny > x1:
                y0, x1 = y0 + 1, x1 - 1
                break
            y, x = ny, nx
    return None

n, h, N = 64, 32, 700
rng = np.random.default_rng(606)
base = []
for _ in range(N):
    B = np.zeros((n, 2), dtype=np.int8)
    for c in range(2):
        B[rng.choice(n, h, replace=False), c] = 1
    base.append(B.mean(1))

print('%-7s %7s %8s %8s %8s  %s' % ('mode', 'rate', 'sd', 'skew', 'kurt', '十箱直方图(%)'))
for mode in ('base', 'jump2', 'jump3', 'slide', 'edit'):
    xs = []
    for a in base:
        d = walk(a, mode)
        if d is not None:
            xs.append(d / (n - 1))
    if len(xs) < 30:
        print('%-7s  样本不足 (%d)' % (mode, len(xs)))
        continue
    x = np.array(xs)
    hist, _ = np.histogram(x, bins=10, range=(0, 1))
    print('%-7s %7.3f %8.4f %8.3f %8.3f  %s'
          % (mode, len(xs) / N, x.std(ddof=1), stats.skew(x), stats.kurtosis(x),
             (hist / hist.sum() * 100).round(0).astype(int)))
