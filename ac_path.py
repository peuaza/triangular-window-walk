# -*- coding: utf-8 -*-
"""ac_path.py —— 沿路径的局部 ac 自相关随 lag 的衰减"""
import io, os, sys
import numpy as np

os.chdir(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ".")

def walk_acseq(B):
    n = B.shape[0]
    a = B.mean(1)
    y0, x1 = 0, n - 1
    seq = []
    while x1 - y0 >= 3:
        sub = a[y0:x1 + 1]
        if sub.std() > 1e-9 and len(sub) > 2:
            v = np.corrcoef(sub[:-1], sub[1:])[0, 1]
            if np.isfinite(v):
                seq.append(v)
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
    return np.array(seq)

rng = np.random.default_rng(1717)
n = 512
seqs = []
for _ in range(300):
    h = (n + 1) // 2
    B = np.zeros((n, 2))
    for j in range(2):
        B[rng.permutation(n)[:h], j] = 1.0
    s = walk_acseq(B)
    if len(s) >= 60:
        seqs.append(s)
print("游走 %d ｜平均序列长 %.0f" % (len(seqs), np.mean([len(s) for s in seqs])))
print("%-6s %-9s" % ("lag", "自相关"))
for lag in (1, 2, 3, 5, 10, 20, 40):
    vals = []
    for s in seqs:
        if len(s) > lag + 5:
            vals.append(np.corrcoef(s[:-lag], s[lag:])[0, 1])
    if vals:
        print("%-6d %-9.3f" % (lag, np.nanmean(vals)), flush=True)
