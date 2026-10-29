# -*- coding: utf-8 -*-
"""quarter.py —— 落回约定下 P(δ) 的高精度测量"""
import io, os, sys
import numpy as np
from concurrent.futures import ProcessPoolExecutor

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

def _b(args):
    n, count, seed = args
    from walk_bin import walk_bin
    rng = np.random.default_rng(seed)
    h = (n + 1) // 2
    hit = 0
    for _ in range(count):
        B = np.zeros((n, 2))
        for j in range(2):
            B[rng.permutation(n)[:h], j] = 1.0
        a = B.mean(1)
        y0, x1 = 0, n - 1
        won = False
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
                    if seg.std() > 1e-9:
                        c = np.corrcoef(seg[:-1], seg[1:])[0, 1]
                        b1 = (a[x] > a[y]) if np.isnan(c) else (c > 0)
                    else:
                        b1 = a[x] > a[y]
                o = ("B" if b1 else "D") if b2 else ("C" if b1 else "A")
                ny, nx = {"A": (y - 1, x), "B": (y + 1, x), "C": (y, x + 1), "D": (y, x - 1)}[o]
                if ny < y0: hh = "top"; break
                if nx > x1: hh = "right"; break
                if ny >= nx: won = True; hh = "hyp"; break
                if nx < y0 or ny > x1: hh = "bad"; break
                y, x = ny, nx
            if hh == "hyp": break
            if hh == "top": x1 = x - 1
            elif hh == "right": y0 = y + 1
            else: y0, x1 = y0 + 1, x1 - 1
        hit += won
    return hit, count

if __name__ == "__main__":
    print("%-6s %-8s %-9s %-9s %s" % ("n", "尝试", "δ", "P(δ)", "距 1/4"))
    for n, per in ((256, 6000), (512, 6000), (1024, 4000), (2048, 2400)):
        with ProcessPoolExecutor(max_workers=8) as ex:
            futs = [ex.submit(_b, (n, per // 8, 95000 + 131 * j)) for j in range(8)]
            hit = tot = 0
            for f in futs:
                a, c = f.result(); hit += a; tot += c
        p = hit / tot; se = np.sqrt(p * (1 - p) / tot)
        print("%-6d %-8d %-9d %-9.4f %-9.4f（%.1f SE）"
              % (n, tot, hit, p, se, (p - 0.25) / se), flush=True)
