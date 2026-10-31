# -*- coding: utf-8 -*-
"""ac_corr3.py —— 相邻段起点的 ac 相关（修正版）"""
import io, os, sys
import numpy as np
from concurrent.futures import ProcessPoolExecutor

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

def _b(args):
    n, count, seed = args
    from walk_bin import walk_bin
    rng = np.random.default_rng(seed)
    h = (n + 1) // 2
    out = []
    for _ in range(count):
        B = np.zeros((n, 2))
        for j in range(2):
            B[rng.permutation(n)[:h], j] = 1.0
        a = B.mean(1)
        y0, x1 = 0, n - 1
        prev = None
        while x1 - y0 >= 3:
            m0 = x1 - y0 + 1
            sub = a[y0:x1 + 1]
            cur = (np.corrcoef(sub[:-1], sub[1:])[0, 1]
                   if sub.std() > 1e-9 and len(sub) > 2 else np.nan)
            if prev is not None and np.isfinite(cur) and np.isfinite(prev):
                out.append((m0, prev, cur))
            prev = cur
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
    return out

if __name__ == "__main__":
    n = 512
    recs = []
    with ProcessPoolExecutor(max_workers=8) as ex:
        futs = [ex.submit(_b, (n, 120, 140000 + 151 * j)) for j in range(8)]
        for f in futs:
            recs.extend(f.result())
    R = np.array(recs)
    m = R[:, 0]; x1 = R[:, 1]; x2 = R[:, 2]
    print("相邻段 ac 对：%d" % len(R))
    print("%-10s %-8s %-9s %-10s %-12s" % ("m 区间", "样本", "rho", "实际切换率", "arccos(rho)/pi"))
    edges = [(4, 8), (8, 16), (16, 32), (32, 64), (64, 128), (128, 256), (256, 512)]
    for a, b in edges:
        sel = (m >= a) & (m < b)
        if sel.sum() < 100:
            continue
        rho = np.corrcoef(x1[sel], x2[sel])[0, 1]
        sw = np.mean(np.sign(x1[sel]) != np.sign(x2[sel]))
        pred = np.arccos(np.clip(rho, -1, 1)) / np.pi
        print("%-10s %-8d %-9.3f %-10.4f %-12.4f"
              % ("%d–%d" % (a, b - 1), sel.sum(), rho, sw, pred), flush=True)
