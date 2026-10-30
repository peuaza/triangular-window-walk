# -*- coding: utf-8 -*-
"""first_switch.py —— 第一次变号时刻的分布与作用"""
import io, os, sys
import numpy as np

os.chdir(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ".")

def walk_dirs(B):
    n = B.shape[0]
    a = B.mean(1)
    y0, x1 = 0, n - 1
    dirs = []
    while x1 - y0 >= 3:
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
                b1 = (np.corrcoef(seg[:-1], seg[1:])[0, 1] > 0) if seg.std() > 1e-9 else (a[x] > a[y])
            o = ("B" if b1 else "D") if b2 else ("C" if b1 else "A")
            ny, nx = {"A": (y - 1, x), "B": (y + 1, x), "C": (y, x + 1), "D": (y, x - 1)}[o]
            if ny < y0: hh = "top"; break
            if nx > x1: hh = "right"; break
            if ny >= nx: hh = "hyp"; break
            if nx < y0 or ny > x1: hh = "bad"; break
            y, x = ny, nx
        if hh == "hyp": break
        if hh == "top": x1 = x - 1; dirs.append(0)
        elif hh == "right": y0 = y + 1; dirs.append(1)
        else: y0, x1 = y0 + 1, x1 - 1
    return np.array(dirs)

rng = np.random.default_rng(3535)
n = 512
rows = []; fano = []
for _ in range(400):
    h = (n + 1) // 2
    B = np.zeros((n, 2))
    for j in range(2):
        B[rng.permutation(n)[:h], j] = 1.0
    d = walk_dirs(B)
    if len(d) < 50:
        continue
    sw = np.where(d[1:] != d[:-1])[0]
    t1 = (sw[0] + 1) / len(d) if len(sw) else 1.0
    if len(sw) >= 2:
        gaps = np.diff(sw)
        fano.append(gaps.var()/max(gaps.mean(),1e-9))
    rows.append((t1, len(sw), d.mean(), len(d)))
R = np.array(rows)
print("游走 %d ｜首次变号时刻（归一化）分位：" % len(R),
      " ".join("%.2f" % x for x in np.percentile(R[:, 0], [5, 25, 50, 75, 95])))
print("无切换（t=1.0）占比：%.2f" % np.mean(R[:, 0] >= 0.999))
print("corr(首变号时刻, 比例)  = %.3f" % np.corrcoef(R[:, 0], R[:, 2])[0, 1])
print("corr(1−时刻, 比例)     = %.3f" % np.corrcoef(1 - R[:, 0], R[:, 2])[0, 1])
print("切换间隔 Fano 因子：中位 %.2f ｜均值 %.2f（Poisson=1）" % (np.median(fano), np.mean(fano)))
print("直方图（10 格）：", " ".join("%.2f" % x for x in np.histogram(R[:, 0], bins=10, range=(0, 1))[0] / len(R)))
