# -*- coding: utf-8 -*-
"""wm_exact.py —— 小 m 的精确胜率 w(m)

枚举所有列和为 h=⌈m/2⌉ 的 m×2 的 0/1 矩阵，对每个跑 walk_bin
"""
import io, os, sys, time
from itertools import combinations
import numpy as np

os.chdir(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ".")
from walk_bin import walk_bin

def enum_w(m):
    h = (m + 1) // 2
    idx = list(combinations(range(m), h))
    tot = 0; hit = 0
    for i1 in idx:
        c1 = np.zeros(m); c1[list(i1)] = 1.0
        for i2 in idx:
            c2 = np.zeros(m); c2[list(i2)] = 1.0
            if walk_bin(np.c_[c1, c2]) is not None:
                hit += 1
            tot += 1
    return tot, hit

if __name__ == "__main__":
    print("%-4s %-8s %-9s %-9s" % ("m", "矩阵数", "δ 数", "w(m)"))
    for m in (12,):
        t0 = time.time()
        tot, hit = enum_w(m)
        print("%-4d %-8d %-9d %-9.4f ｜%.0fs" % (m, tot, hit, hit / tot, time.time() - t0), flush=True)
