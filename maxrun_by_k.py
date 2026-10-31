# -*- coding: utf-8 -*-
"""maxrun_by_k.py —— 最长段占比按切换数分组"""
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
        if hh == "top": x1 = x - 1; dirs.append(0)
        elif hh == "right": y0 = y + 1; dirs.append(1)
        else: y0, x1 = y0 + 1, x1 - 1
    return np.array(dirs)

rng = np.random.default_rng(3434)
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
    k = int(np.sum(d[1:] != d[:-1]))
    runs = []; r = 1
    for i in range(1, len(d)):
        if d[i] == d[i - 1]:
            r += 1
        else:
            runs.append(r); r = 1
    runs.append(r)
    rows.append((k, max(runs) / len(d), d.mean(), len(d)))
R = np.array(rows)
print("游走 %d" % len(R))
print("%-9s %-6s %-14s %-12s %-10s" % ("组", "n", "最长段占比中位", "比例 sd", "长度中位"))
for lo, hi, tag in ((0, 0, "k=0"), (1, 3, "k=1–3"), (4, 9, "k=4–9"), (10, 999, "k≥10")):
    sel = (R[:, 0] >= lo) & (R[:, 0] <= hi)
    if sel.sum() > 5:
        print("%-9s %-6d %-14.2f %-12.3f %-10.0f"
              % (tag, sel.sum(), np.median(R[sel, 1]), np.std(R[sel, 2]), np.median(R[sel, 3])))
print()
print("全体：最长段占比中位 %.2f ｜比例 sd %.3f" % (np.median(R[:, 1]), np.std(R[:, 2])))
