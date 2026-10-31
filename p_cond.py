# -*- coding: utf-8 -*-
"""p_cond.py —— p 对 ac 的条件均值与条件 sd（分箱）"""
import io, os, sys
import numpy as np

os.chdir(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ".")

def walk_info(B):
    n = B.shape[0]
    a = B.mean(1)
    ac = np.corrcoef(a[:-1], a[1:])[0, 1] if a.std() > 1e-9 else 0.0
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
    return ac, (np.mean(dirs) if len(dirs) >= 30 else np.nan)

rng = np.random.default_rng(2323)
n = 512
acs, ps = [], []
for _ in range(1500):
    h = (n + 1) // 2
    B = np.zeros((n, 2))
    for j in range(2):
        B[rng.permutation(n)[:h], j] = 1.0
    ac, p = walk_info(B)
    if np.isfinite(p):
        acs.append(ac); ps.append(p)
acs = np.array(acs); ps = np.array(ps)
q = np.quantile(acs, np.linspace(0, 1, 9))
print("ac 分箱（8 档）：")
print("%-16s %-6s %-9s %-9s" % ("ac 区间", "样本", "p 均值", "p 的 sd"))
for i in range(8):
    m = (acs >= q[i]) & (acs <= q[i + 1] if i == 7 else acs < q[i + 1])
    if m.sum() < 30:
        continue
    print("%-16s %-6d %-9.3f %-9.3f"
          % ("[%.3f,%.3f]" % (q[i], q[i + 1]), m.sum(), ps[m].mean(), ps[m].std()), flush=True)
