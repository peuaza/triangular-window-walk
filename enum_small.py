# -*- coding: utf-8 -*-
"""enum_small.py —— 有限样本空间的精确枚举

输入空间：n×2 的 0/1 矩阵，每列恰有 h=⌈n/2⌉ 个 1
总数：C(n,h)^2
对每个矩阵跑一次游走，精确统计 P(δ)、位置分布
"""
import io, json, os
from itertools import combinations
import numpy as np

os.chdir(os.path.dirname(os.path.abspath(__file__)))

def b1_of(a, B, y, x):
    m = x - y + 1
    if m == 1: return a[y] > 0.5
    if m == 2: return a[x] > a[y]
    seg = a[y:x + 1]
    if seg.std() < 1e-9: return a[x] > a[y]
    return np.corrcoef(seg[:-1], seg[1:])[0, 1] > 0

def walk(B):
    """B: n×2 的 0/1 矩阵（原样使用，不再二值化）"""
    n = B.shape[0]
    a = B.mean(1)
    y0, x1 = 0, n - 1; d = None
    while x1 - y0 >= 3:
        y, x = y0, x1; hit = None; seen = set()
        while True:
            if (y, x) in seen: hit = "cyc"; break
            seen.add((y, x))
            m = x - y + 1
            b2 = a[y:x + 1].mean() > a[y:y + max(1, m // 2)].mean()
            b1 = b1_of(a, B, y, x)
            o = ("B" if b1 else "D") if b2 else ("C" if b1 else "A")
            ny, nx = {"A": (y - 1, x), "B": (y + 1, x), "C": (y, x + 1), "D": (y, x - 1)}[o]
            if ny < y0: hit = "top"; break
            if nx > x1: hit = "right"; break
            if ny >= nx: hit = "hyp"; d = ny; break
            if nx < y0 or ny > x1: hit = "bad"; break
            y, x = ny, nx
        if hit == "hyp": break
        if hit == "top": x1 = x - 1
        elif hit == "right": y0 = y + 1
        else: y0, x1 = y0 + 1, x1 - 1
    return None if d is None else d / (n - 1)

def enumerate_n(n):
    h = (n + 1) // 2
    idx = list(combinations(range(n), h))
    tot = 0; hit = 0; pos = []; means = []
    for i1 in idx:
        c1 = np.zeros(n); c1[list(i1)] = 1.0
        for i2 in idx:
            c2 = np.zeros(n); c2[list(i2)] = 1.0
            B = np.c_[c1, c2]
            v = walk(B); tot += 1
            if v is not None: hit += 1; pos.append(v)
            means.append(B.mean())
    return tot, hit, np.array(pos), np.array(means), h

print("%-4s %-6s %-10s %-9s %-9s %s" % ("n", "h", "矩阵总数", "δ 数", "δ 率", "位置中位"))
out = []
for n in (6, 8, 10, 12):
    tot, hit, pos, means, h = enumerate_n(n)
    out.append({"n": n, "h": h, "total": tot, "hit": hit, "rate": hit / tot,
                "pos_med": float(np.median(pos)) if len(pos) else None,
                "pos_sd": float(pos.std()) if len(pos) else None,
                "mean_all": float(means.mean())})
    print("%-4d %-6d %-10d %-9d %-9.4f %s"
          % (n, h, tot, hit, hit / tot, ("%.3f" % np.median(pos)) if len(pos) else "-"), flush=True)
json.dump(out, io.open("_enum_small.json", "w"))
print("→ _enum_small.json")
