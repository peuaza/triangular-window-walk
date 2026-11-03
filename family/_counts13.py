# -*- coding: utf-8 -*-
"""类级与矩阵级计数：M(n,k) = #{a: δ=k}，N(n,k) = Σ_{δ(a)=k} C(m,m/2)，n=4..12"""
import itertools
import json
import time
import numpy as np
from math import comb

def delta_of(a):
    n = len(a)
    y0, x1 = 0, n - 1
    while x1 - y0 >= 3:
        y, x = y0, x1
        seen = set()
        while True:
            if (y, x) in seen:
                y0, x1 = y0 + 1, x1 - 1
                break
            seen.add((y, x))
            m = x - y + 1
            s = a[y:x + 1]
            b2 = s.mean() > s[:max(1, m // 2)].mean()
            if m == 1:
                b1 = a[y] > 0.5
            elif m == 2:
                b1 = a[x] > a[y]
            else:
                u, v = s[:-1], s[1:]
                b1 = (np.corrcoef(u, v)[0, 1] > 0) if (u.std() > 1e-12 and v.std() > 1e-12) else False
            o = ("B" if b1 else "D") if b2 else ("C" if b1 else "A")
            ny, nx = {"A": (y - 1, x), "B": (y + 1, x), "C": (y, x + 1), "D": (y, x - 1)}[o]
            if ny < y0:
                x1 = x - 1
                break
            if nx > x1:
                y0 = y + 1
                break
            if ny >= nx:
                return ny
            if nx < y0 or ny > x1:
                y0, x1 = y0 + 1, x1 - 1
                break
            y, x = ny, nx
    return None

out = {}
for n in [13]:
    h = (n + 1) // 2
    t0 = time.time()
    M = [0] * n
    N = [0] * n
    tot_seq = 0
    for tup in itertools.product((0.0, 0.5, 1.0), repeat=n):
        m = sum(1 for v in tup if v == 0.5)
        n1 = sum(1 for v in tup if v == 1.0)
        if 2 * n1 + m != 2 * h:
            continue
        tot_seq += 1
        d = delta_of(np.array(tup))
        if d is not None:
            M[d] += 1
            N[d] += comb(m, m // 2)
    out[n] = dict(M=M, N=N, seq=tot_seq, sec=round(time.time() - t0, 1))
    print('n=%2d ｜ 合法序列 %6d ｜ M(k)=%s ｜ N(k)=%s ｜ %.0fs'
          % (n, tot_seq, M, N, time.time() - t0), flush=True)

json.dump({str(k): v for k, v in out.items()}, open('_counts13.json', 'w'), indent=1)
print('saved _counts13.json')
