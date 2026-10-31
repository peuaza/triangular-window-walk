# -*- coding: utf-8 -*-
"""excursion.py —— |ac| 低于阈值的停留段与间隔"""
import io, os, sys
import numpy as np

os.chdir(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ".")

def ac_seq(B):
    n = B.shape[0]
    a = B.mean(1)
    y0, x1 = 0, n - 1
    out = []
    while x1 - y0 >= 3:
        sub = a[y0:x1 + 1]
        v = (np.corrcoef(sub[:-1], sub[1:])[0, 1]
             if sub.std() > 1e-9 and len(sub) > 2 else np.nan)
        if np.isfinite(v):
            out.append(abs(v))
        y, x = y0, x1
        seen = set(); hh = None
        while True:
            if (y, x) in seen:
                hh = "cyc"; break
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
                b1 = (np.corrcoef(u, v)[0, 1] > 0) if (u.std() > 1e-12 and v.std() > 1e-12) else False
            o = ("B" if b1 else "D") if b2 else ("C" if b1 else "A")
            ny, nx = {"A": (y - 1, x), "B": (y + 1, x), "C": (y, x + 1), "D": (y, x - 1)}[o]
            if ny < y0: hh = "top"; break
            if nx > x1: hh = "right"; break
            if ny >= nx: hh = "hyp"; break
            if nx < y0 or ny > x1: hh = "bad"; break
            y, x = ny, nx
        if hh == "hyp": break
        if hh == "top": x1 = x - 1
        elif hh == "right": y0 = y + 1
        else: y0, x1 = y0 + 1, x1 - 1
    return np.array(out)

rng = np.random.default_rng(3737)
n = 512
eps = 0.013
lens, gaps, tot = [], [], 0
for _ in range(250):
    h = (n + 1) // 2
    B = np.zeros((n, 2))
    for j in range(2):
        B[rng.permutation(n)[:h], j] = 1.0
    s = ac_seq(B)
    if len(s) < 50:
        continue
    low = s < eps
    # 游程编码
    i = 0
    while i < len(low):
        j = i
        while j < len(low) and low[j] == low[i]:
            j += 1
        (lens if low[i] else gaps).append(j - i)
        i = j
    tot += len(s)
lens = np.array(lens); gaps = np.array(gaps)
print("阈值 eps=%.3f ｜总步 %d" % (eps, tot))
print("停留段（|ac|<eps）：n=%d ｜中位 %.1f ｜均值 %.2f" % (len(lens), np.median(lens), lens.mean()))
print("间隔段（|ac|>=eps）：n=%d ｜中位 %.1f ｜均值 %.2f" % (len(gaps), np.median(gaps), gaps.mean()))
print("处于停留段的比例：%.3f" % (lens.sum() / (lens.sum() + gaps.sum())))
print("Fano(停留段长度) = %.2f" % (lens.var() / max(lens.mean(), 1e-9)))
