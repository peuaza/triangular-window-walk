# -*- coding: utf-8 -*-
"""相变点实测：单遍建表 + k 检查点

建表：一遍扫过纸带，在每个检查点位置存 (S1,S2,S3) ⇒ O(n) 读
查询：从 cp(x) 扫到 x，≤ step = n/k 格 ⇒ O(n/k)
总读 ≈ n + cuts × n/k     （cuts ≈ n/2）
"""
import sys

import numpy as np

N = int(sys.argv[1]) if len(sys.argv) > 1 else 128
NW = int(sys.argv[2]) if len(sys.argv) > 2 else 20
h = N // 2


def walk(a, n, k):
    step = max(1, n // k)
    # ---- 单遍建表 ----
    pref = {}
    s1 = s2 = s3 = 0.0
    prev = None
    reads = 0
    nxt = 0
    for i in range(n):
        if i == nxt:
            pref[i] = (s1, s2, s3)
            nxt += step
        v = a[i]
        reads += 1
        s1 += v; s2 += v * v
        s3 += 0.0 if prev is None else prev * v
        prev = v
    pref[n] = (s1, s2, s3)

    def pref_at(x):
        nonlocal reads
        c = (x // step) * step
        if c not in pref:
            c = max([p for p in pref if p <= x])
        b1, b2, b3 = pref[c]
        prev = a[c - 1] if c > 0 else None
        for i in range(c, x):
            v = a[i]; reads += 1
            b1 += v; b2 += v * v
            b3 += 0.0 if prev is None else prev * v
            prev = v
        return b1, b2, b3

    def seg(lo, hi):
        if hi < lo:
            return 0.0, 0.0, 0.0
        a1, a2, _ = pref_at(hi + 1)
        b1, b2, _ = pref_at(lo)
        _, _, c3 = pref_at(hi)
        _, _, d3 = pref_at(lo)
        return a1 - b1, a2 - b2, c3 - d3

    L, R = 0, n - 1
    ell, r = 0, n - 1
    cuts = 0
    steps = 0
    while True:
        if ell >= r:
            return ell / (n - 1), reads, cuts
        if R - L + 1 < 3:
            return None, reads, cuts
        steps += 1
        if steps > 20 * n:
            return None, reads, cuts
        m = r - ell + 1
        kk = max(1, m // 2)
        s1w, s2w, s3w = seg(ell, r)
        sh, _, _ = seg(ell, ell + kk - 1)
        b2 = (s1w / m) > (sh / kk)
        if m >= 3:
            nu = m - 1
            mu = (s1w - a[r]) / nu
            mv = (s1w - a[ell]) / nu
            vu = (s2w - a[r] ** 2) / nu - mu * mu
            vv = (s2w - a[ell] ** 2) / nu - mv * mv
            b1 = False if (vu <= 1e-24 or vv <= 1e-24) else ((s3w / nu - mu * mv) / np.sqrt(vu * vv) > 0)
        elif m == 1:
            b1 = a[ell] > 0.5
        else:
            b1 = a[r] > a[ell]
        mv_ = {(True, True): (1, 0), (True, False): (0, -1),
               (False, True): (0, 1), (False, False): (-1, 0)}[(b2, b1)]
        ny, nx = ell + mv_[0], r + mv_[1]
        if ny < L or nx > R:
            if ny < L:
                R = r - 1
            else:
                L = ell + 1
            cuts += 1
            ell, r = L, R
            if R - L + 1 < 3:
                return None, reads, cuts
            continue
        if ny >= nx:
            return ny / (n - 1), reads, cuts
        if nx < L or ny > R:
            L, R = L + 1, R - 1
            cuts += 1
            ell, r = L, R
            if R - L + 1 < 3:
                return None, reads, cuts
            continue
        ell, r = ny, nx


rng = np.random.default_rng(31)
print('n=%d' % N)
print('%-6s %10s %8s %8s %8s %10s' % ('k', '总读', '读/n', '读/n²', '切割', 'n+n²/(2k)'))
for k in (1, 2, 4, 8, 16, 32, 64, 128):
    if k > N:
        continue
    R, C = [], []
    for _ in range(NW):
        B = np.zeros((N, 2), dtype=np.int8)
        for c in range(2):
            B[rng.choice(N, h, replace=False), c] = 1
        d, rd, cu = walk(B.mean(1), N, k)
        if d is not None:
            R.append(rd); C.append(cu)
    if not R:
        print('%-6d 无输出' % k); continue
    mr, mc = np.mean(R), np.mean(C)
    print('%-6d %10.0f %8.1f %8.3f %8.1f %10.0f'
          % (k, mr, mr / N, mr / N / N, mc, N + N * N / (2 * k)))
