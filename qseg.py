# -*- coding: utf-8 -*-
"""qseg.py —— 每段胜率 q(m) 与 (m, 是否赢) 的配对

对每场游走，记录每段的起始 m 与该段是否赢
"""
import io, os, sys
import numpy as np
from concurrent.futures import ProcessPoolExecutor

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

def walk_pairs(B):
    n = B.shape[0]
    a = B.mean(1)
    y0, x1 = 0, n - 1
    pairs = []
    while x1 - y0 >= 3:
        m0 = x1 - y0 + 1
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
            if ny >= nx: hit = "hyp"; break
            if nx < y0 or ny > x1: hit = "bad"; break
            y, x = ny, nx
        pairs.append((m0, hit == "hyp"))
        if hit == "hyp": break
        if hit == "top": x1 = x - 1
        elif hit == "right": y0 = y + 1
        else: y0, x1 = y0 + 1, x1 - 1
    return pairs

def _batch(args):
    n, count, seed = args
    rng = np.random.default_rng(seed)
    out = []
    for _ in range(count):
        h = (n + 1) // 2
        B = np.zeros((n, 2))
        for j in range(2):
            B[rng.permutation(n)[:h], j] = 1.0
        out.extend(walk_pairs(B))
    return out

if __name__ == "__main__":
    n = 1024
    recs = []
    with ProcessPoolExecutor(max_workers=8) as ex:
        futs = [ex.submit(_batch, (n, 20, 2500 + 29 * j)) for j in range(8)]
        for f in futs:
            recs.extend(f.result())
    m = np.array([r[0] for r in recs]); w = np.array([r[1] for r in recs])
    print("n=%d ｜段 %d ｜总胜 %d ｜P(δ)=%.4f" % (n, len(m), w.sum(), w.mean()))
    print("Σ q = %.4f（−ln 0.685 = 0.378）" % w.sum())
    print()
    print("%-10s %-7s %-9s %-9s" % ("m 区间", "段数", "q", "q·段数"))
    edges = [4, 8, 12, 20, 32, 64, 128, 256, 512, 1100]
    for a, b in zip(edges[:-1], edges[1:]):
        sel = (m >= a) & (m < b)
        if sel.sum() < 50:
            continue
        print("%-10s %-7d %-9.4f %-9.2f" % ("%d–%d" % (a, b - 1), sel.sum(), w[sel].mean(), w[sel].sum()), flush=True)
