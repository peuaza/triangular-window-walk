# -*- coding: utf-8 -*-
"""arcsine_conv.py —— δ 分布向反正弦收敛？"""
import io, os, sys
import numpy as np
from scipy import stats
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
        v = walk_bin(B)
        if v is not None:
            out.append(v)
    return out

def arcsine_cdf(x):
    return 2 / np.pi * np.arcsin(np.sqrt(np.clip(x, 0, 1)))

if __name__ == "__main__":
    print("%-6s %-6s %-9s %-9s %-9s %-9s" % ("n", "δ", "sd", "中位", "KS(arcsine)", "KS(Beta .6)"))
    for n, per in ((48, 2400), (64, 2400), (96, 2400), (128, 2400), (256, 1600), (512, 1200)):
        recs = []
        with ProcessPoolExecutor(max_workers=8) as ex:
            futs = [ex.submit(_b, (n, per // 8, 41000 + 101 * j)) for j in range(8)]
            for f in futs:
                recs.extend(f.result())
        X = np.array(recs)
        d1 = stats.kstest(X, arcsine_cdf).statistic
        d2 = stats.kstest(X, "beta", args=(0.6, 0.6, 0, 1)).statistic
        print("%-6d %-6d %-9.4f %-9.3f %-9.4f %-9.4f"
              % (n, len(X), X.std(), np.median(X), d1, d2), flush=True)
