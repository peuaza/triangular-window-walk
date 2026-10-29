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
rows = []
for _ in range(400):
    h = (n + 1) // 2
    B = np.zeros((n, 2))
    for j in range(2):
        B[rng.permutation(n)[:h], j] = 1.0
    d = walk_dirs(B)
    if len(d) < 50:
        continue
    sw = np.where(d[1:] != d[:-1])[0]
    p = (sw.mean() + 1) / len(d) if len(sw) else 0.5
    p = (sw.mean() + 1) / len(d) if len(sw) else 0.5
    rows.append((p, len(sw), d.mean(), len(d)))
R = np.array(rows)
sel = R[:, 1] > 0
P = R[sel, 0]; F = R[sel, 2]
print("有切换的游走 %d ｜切换重心分位：" % sel.sum(),
      " ".join("%.2f" % x for x in np.percentile(P, [5, 25, 50, 75, 95])))
print("corr(p, 比例) = %.3f ｜corr(1−p, 比例) = %.3f" % (np.corrcoef(P, F)[0,1], np.corrcoef(1-P, F)[0,1]))
pred = np.where(np.abs(F - P) < np.abs(F - (1-P)), P, 1-P)
print("|比例 − 最近端点(p 或 1−p)| 中位 %.3f ｜均值 %.3f" % (np.median(np.abs(F-pred)), np.mean(np.abs(F-pred))))
print("直方图（10 格）：", " ".join("%.2f" % x for x in np.histogram(R[:, 0], bins=10, range=(0, 1))[0] / len(R)))
