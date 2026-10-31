# -*- coding: utf-8 -*-
"""cluster_src.py —— 切换率是否由 |ac| 的大小决定"""
import io, os, sys
import numpy as np

os.chdir(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ".")

def walk_seq(B):
    n = B.shape[0]
    a = B.mean(1)
    y0, x1 = 0, n - 1
    rec = []   # (|ac|, 是否在下一步切换)
    prev = None
    while x1 - y0 >= 3:
        sub = a[y0:x1 + 1]
        cur = (np.corrcoef(sub[:-1], sub[1:])[0, 1]
               if sub.std() > 1e-9 and len(sub) > 2 else np.nan)
        y, x = y0, x1
        seen = set(); hh = None; di = None
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
            if ny < y0: hh = "top"; di = 0; break
            if nx > x1: hh = "right"; di = 1; break
            if ny >= nx: hh = "hyp"; break
            if nx < y0 or ny > x1: hh = "bad"; break
            y, x = ny, nx
        if hh == "hyp": break
        if np.isfinite(cur) and prev is not None:
            rec.append((abs(prev), 1 if di != prev_di else 0))
        prev = cur; prev_di = di
        if hh == "top": x1 = x - 1
        elif hh == "right": y0 = y + 1
        else: y0, x1 = y0 + 1, x1 - 1
    return rec

rng = np.random.default_rng(3636)
n = 512
allrec = []
for _ in range(250):
    h = (n + 1) // 2
    B = np.zeros((n, 2))
    for j in range(2):
        B[rng.permutation(n)[:h], j] = 1.0
    allrec.extend(walk_seq(B))
R = np.array(allrec)
print("记录 %d 步" % len(R))
q = np.quantile(R[:, 0], np.linspace(0, 1, 6))
print("%-16s %-8s %-10s" % ("|ac| 区间", "样本", "切换率"))
for i in range(5):
    sel = (R[:, 0] >= q[i]) & (R[:, 0] <= q[i + 1] if i == 4 else R[:, 0] < q[i + 1])
    if sel.sum() > 50:
        print("%-16s %-8d %-10.3f" % ("[%.3f,%.3f]" % (q[i], q[i + 1]), sel.sum(), R[sel, 1].mean()))
