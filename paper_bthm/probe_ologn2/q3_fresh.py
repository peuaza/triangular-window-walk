# -*- coding: utf-8 -*-
"""q3: cost law of a single *fresh* window (no trajectory), b1=0 scan.

For a window [L,R] of length M the epoch cost is
    cost = M - m* + 1   (m* = largest m with A_{floor(m/2)} >= A_m)
    cost = M - 1        if no such m ("meet": the run would output L).
This script measures E[cost]/M, P(cost/M >= c), P(meet), and the
dependence on the *top* comparison
    d = sum(a[R-k+1..R]) - sum(a[L..R-k]),  k = floor(M/2)
which is the quantity the martingale picture says controls the depth:
the scan walks k downwards and the test at m=2k is  sum[k,m) > sum[0,k).

Windows are taken from long random ensemble tapes, so their content is a
uniform random sub-multiset arrangement (limit of the ensemble for length M).
"""
import sys
import numpy as np
from q0_walk import prefix, epoch, sample_tapes


def main():
    n = 8192
    a = next(sample_tapes(n, 1, 7))
    P = prefix(a)
    rng = np.random.default_rng(5)
    print('%8s %8s %9s %9s %9s %9s %9s %9s %9s'
          % ('M', 'windows', 'E[cost]/M', 'med', 'P(.05)', 'P(.1)', 'P(.25)', 'P(meet)', 'E[d^2]/M'))
    for M in (16, 32, 64, 128, 256, 512, 1024, 2048, 4096):
        reps = max(200, int(2e6 / M))
        ratios, meets, d2, costs = [], [], [], []
        for _ in range(reps):
            L = int(rng.integers(0, n - M + 1))
            R = L + M - 1
            ms, cost = epoch(P, L, R)
            ratios.append(cost / M)
            costs.append(cost)
            meets.append(ms is None)
            k = M // 2
            s1 = int(P[L + k] - P[L])
            s2 = int(P[L + M] - P[L + k])
            d2.append((s2 - s1) ** 2)
        ratios = np.array(ratios)
        print('%8d %8d %9.4f %9.4f %9.3f %9.3f %9.3f %9.4f %9.3f'
              % (M, reps, ratios.mean(), np.median(ratios),
                 np.mean(ratios >= .05), np.mean(ratios >= .1),
                 np.mean(ratios >= .25), np.mean(np.array(meets)),
                 np.mean(d2) / M))
        sys.stdout.flush()


if __name__ == '__main__':
    main()
