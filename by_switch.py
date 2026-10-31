# -*- coding: utf-8 -*-
"""by_switch.py —— 按切换次数分解 δ 的分布"""
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
        dirs = []
        won = None
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
                if ny >= nx: won = ny / (n - 1); hh = "hyp"; break
                if nx < y0 or ny > x1: hh = "bad"; break
                y, x = ny, nx
            if hh == "hyp": break
            if hh == "top": x1 = x - 1; dirs.append(0)
            elif hh == "right": y0 = y + 1; dirs.append(1)
            else: y0, x1 = y0 + 1, x1 - 1
        if won is not None and len(dirs) >= 10:
            d = np.array(dirs)
            out.append((int(np.sum(d[1:] != d[:-1])), won))
    return out

if __name__ == "__main__":
    n = 512
    recs = []
    with ProcessPoolExecutor(max_workers=8) as ex:
        futs = [ex.submit(_b, (n, 250, 160000 + 163 * j)) for j in range(8)]
        for f in futs:
            recs.extend(f.result())
    R = np.array(recs)
    k = R[:, 0].astype(int); X = R[:, 1]
    print("n=%d ｜δ %d ｜平均切换 %.1f" % (n, len(X), k.mean()))
    print("%-8s %-7s %-8s %-9s %-9s" % ("切换数", "样本", "占比", "δ 均值", "δ 的 sd"))
    for kk in (0, 1, 2, 3, 4, 5, 6, 8, 10):
        sel = k == kk if kk < 10 else k >= 10
        if sel.sum() < 10:
            continue
        print("%-8s %-7d %-8.3f %-9.3f %-9.3f"
              % (("≥10" if kk == 10 else str(kk)), sel.sum(), sel.mean(), X[sel].mean(),
                 X[sel].std()), flush=True)
