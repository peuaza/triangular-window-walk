# -*- coding: utf-8 -*-
"""segments.py —— 把游走分解成「段」，量每段成功率 p 与段数 S

一段 = 两次切分之间（或从起点到首次切分/抵达对角线）
"""
import io, os, sys, time
import numpy as np
from concurrent.futures import ProcessPoolExecutor

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

def walk_seg(A):
    n = A.shape[0]
    B = (A > np.median(A, 0)).astype(float)
    a = B.mean(1)
    y0, x1 = 0, n - 1
    d = None
    segs = 0      # 总段数
    wins = 0      # 抵达对角线的段数（0 或 1）
    while x1 - y0 >= 3:
        y, x = y0, x1
        seen = set(); hit = None
        while True:
            if (y, x) in seen:
                hit = "cyc"; break
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
            if ny < y0: hit = "top"; break
            if nx > x1: hit = "right"; break
            if ny >= nx: hit = "hyp"; d = ny; break
            if nx < y0 or ny > x1: hit = "bad"; break
            y, x = ny, nx
        segs += 1
        if hit == "hyp":
            wins += 1
            break
        if hit == "top": x1 = x - 1
        elif hit == "right": y0 = y + 1
        else: y0, x1 = y0 + 1, x1 - 1
    return (d is not None), segs, wins

def _batch(args):
    n, count, seed = args
    rng = np.random.default_rng(seed)
    tot = win = segsum = 0
    for _ in range(count):
        ok, s, w = walk_seg(rng.random((n, 2)))
        tot += 1; segsum += s; win += w
    return tot, win, segsum

if __name__ == "__main__":
    print("%-6s %-6s %-8s %-9s %-10s %-9s" % ("n", "样本", "P(δ)", "平均段数 S", "每段 p", "p·S"))
    for n in (256, 512, 1024, 2048, 4096):
        workers, chunk = 8, 80
        tot = win = segs = 0
        with ProcessPoolExecutor(max_workers=workers) as ex:
            futs = [ex.submit(_batch, (n, chunk, 500 + 37 * j)) for j in range(workers)]
            for f in futs:
                t, w, s = f.result()
                tot += t; win += w; segs += s
        S = segs / tot
        p = win / segs
        print("%-6d %-6d %-8.4f %-9.1f %-10.5f %-9.4f"
              % (n, tot, win / tot, S, p, p * S), flush=True)
