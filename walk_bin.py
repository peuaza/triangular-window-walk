# -*- coding: utf-8 -*-
"""walk_bin.py —— 直接吃 0/1 矩阵的游走（不重新二值化）"""
import io, os, sys
import numpy as np

os.chdir(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ".")

def walk_bin(B):
    """B: n×2 的 0/1 矩阵，直接用"""
    n = B.shape[0]
    a = B.mean(1)
    y0, x1 = 0, n - 1
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
                u, v = seg[:-1], seg[1:]
                if u.std() < 1e-12 or v.std() < 1e-12:
                    # 约定 A（论文 Remark 1）：自相关无定义时取 b1 = 0
                    b1 = False
                else:
                    b1 = np.corrcoef(u, v)[0, 1] > 0
            o = ("B" if b1 else "D") if b2 else ("C" if b1 else "A")
            ny, nx = {"A": (y - 1, x), "B": (y + 1, x), "C": (y, x + 1), "D": (y, x - 1)}[o]
            if ny < y0: hit = "top"; break
            if nx > x1: hit = "right"; break
            if ny >= nx: hit = "hyp"; return ny / (n - 1)
            if nx < y0 or ny > x1: hit = "bad"; break
            y, x = ny, nx
        if hit == "top": x1 = x - 1
        elif hit == "right": y0 = y + 1
        else: y0, x1 = y0 + 1, x1 - 1
    return None
