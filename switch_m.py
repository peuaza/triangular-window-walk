# -*- coding: utf-8 -*-
"""switch_m.py —— 切换率随区间大小 m 的变化"""
import io, os, sys
import numpy as np
from concurrent.futures import ProcessPoolExecutor
from collections import defaultdict

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
        dirs = []      # (m_at_cut, direction)
        while x1 - y0 >= 3:
            m0 = x1 - y0 + 1
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
                    b1 = (np.corrcoef(seg[:-1], seg[1:])[0, 1] > 0) if seg.std() > 1e-9 else (a[x] > a[y])
                o = ("B" if b1 else "D") if b2 else ("C" if b1 else "A")
                ny, nx = {"A": (y - 1, x), "B": (y + 1, x), "C": (y, x + 1), "D": (y, x - 1)}[o]
                if ny < y0: hh = "top"; break
                if nx > x1: hh = "right"; break
                if ny >= nx: hh = "hyp"; break
                if nx < y0 or ny > x1: hh = "bad"; break
                y, x = ny, nx
            if hh == "hyp": break
            if hh == "top": x1 = x - 1; dirs.append((m0, 0))
            elif hh == "right": y0 = y + 1; dirs.append((m0, 1))
            else: y0, x1 = y0 + 1, x1 - 1
        out.append(dirs)
    return out

if __name__ == "__main__":
    n = 512
    seqs = []
    with ProcessPoolExecutor(max_workers=8) as ex:
        futs = [ex.submit(_b, (n, 200, 90000 + 127 * j)) for j in range(8)]
        for f in futs:
            seqs.extend(f.result())
    # 按 m 分箱统计切换率
    sw = defaultdict(lambda: [0, 0])
    for s in seqs:
        for k in range(1, len(s)):
            m = s[k][0]
            b = min(9, int(np.log2(max(4, m)) * 1.5))
            chg = 1 if s[k][1] != s[k - 1][1] else 0
            sw[b][0] += chg; sw[b][1] += 1
    print("游走 %d ｜切换率按 m 分箱：" % len(seqs))
    print("%-14s %-9s %s" % ("m 范围", "样本", "切换率"))
    edges = [(4, 8), (8, 16), (16, 32), (32, 64), (64, 128), (128, 256), (256, 512)]
    for a, b in edges:
        hit = tot = 0
        for s in seqs:
            for k in range(1, len(s)):
                if a <= s[k][0] < b:
                    tot += 1
                    hit += (s[k][1] != s[k - 1][1])
        if tot >= 100:
            print("%-14s %-9d %.4f" % ("%d–%d" % (a, b - 1), tot, hit / tot), flush=True)
