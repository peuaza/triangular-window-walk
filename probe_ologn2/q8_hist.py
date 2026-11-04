# -*- coding: utf-8 -*-
"""q8: (a) recursion check T(a_n) = cost_0 + T(a[L'..R']); (b) T/n^2 histogram
(bimodality) at n=1024; (c) a tight large-n value of E[T]/n^2 with many tapes.
"""
import sys
import numpy as np
from q0_walk import prefix, epoch, walk_b0_fast, sample_tapes


def check_recursion(ns=(200,)):
    for n in ns:
        bad = 0
        tot = 0
        for a in sample_tapes(n, 200, 11):
            P = prefix(a)
            d, T, diag = walk_b0_fast(a, want_diag=True)
            L0, R0, M0, ms, c0 = diag[0]
            if ms is None:
                ok = (T == c0)
            else:
                R1 = R0 - 1 if ms == M0 else R0
                d2, T2, _ = walk_b0_fast(a[1:R1 + 1])
                ok = (T == c0 + T2)
            tot += 1
            bad += (not ok)
        print('(a) recursion T = cost_0 + T(sub-tape): %d tapes, violations %d' % (tot, bad))


def hist(n=1024, ntapes=2000, nb=20):
    Ts = []
    for a in sample_tapes(n, ntapes, 424242):
        _, T, _ = walk_b0_fast(a)
        Ts.append(T)
    Ts = np.array(Ts, dtype=float) / n ** 2
    print('\n(b) T/n^2 histogram, n=%d, %d tapes (mean %.4f, sd %.4f, se %.4f)'
          % (n, ntapes, Ts.mean(), Ts.std(ddof=1), Ts.std(ddof=1) / np.sqrt(ntapes)))
    edges = np.array([0, 0.0005, 0.001, 0.002, 0.005, 0.01, 0.02, 0.04, 0.06, 0.08,
                      0.10, 0.13, 0.16, 0.20, 0.25, 1.01])
    h, _ = np.histogram(Ts, bins=edges)
    for i in range(len(edges) - 1):
        print('    [%.4f,%.4f): %5d  %5.3f  %s' % (edges[i], edges[i + 1], h[i], h[i] / ntapes,
                                                   '#' * int(60 * h[i] / ntapes)))
    print('    P(T<=n)=%.3f  P(T>=0.05n^2)=%.3f  P(T>=0.10n^2)=%.3f  median=%.4f'
          % (np.mean(Ts <= 1.0 / n), np.mean(Ts >= 0.05), np.mean(Ts >= 0.10), np.median(Ts)))


if __name__ == '__main__':
    check_recursion()
    hist()
