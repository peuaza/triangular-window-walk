# -*- coding: utf-8 -*-
"""q1: step-count statistics of the b1=0 (concavity-scan) walk, random tapes.

Reports, per n: mean/median/quantiles of T (b2 tests), T/n^2, number of epochs,
and the distribution of the per-epoch cost normalised by the window length
M (only epochs with M >= 32 to avoid small-M noise).
"""
import sys
import numpy as np
from q0_walk import walk_b0_fast, sample_tapes


def run(n, ntapes, seed, mmin=32):
    Ts, eps_, ratios, Ms = [], [], [], []
    epfrac = []  # fraction of tests per tape
    for a in sample_tapes(n, ntapes, seed):
        d, T, diag = walk_b0_fast(a, want_diag=True)
        Ts.append(T)
        eps_.append(len(diag))
        for (L, R, M, ms, cost) in diag:
            if M >= mmin:
                ratios.append(cost / M)
                Ms.append(M)
        if T:
            epfrac.append(sum(c for (_, _, _, _, c) in diag) / T)
    Ts = np.array(Ts, dtype=float)
    ratios = np.array(ratios)
    Ms = np.array(Ms, dtype=float)
    q = np.quantile(Ts / n ** 2, [0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.99])
    print('n=%5d tapes=%5d  T: mean %9.1f  T/n^2 mean %.4f med %.4f  '
          'epochs mean %.1f (=%.3f n)'
          % (n, ntapes, Ts.mean(), (Ts / n ** 2).mean(), np.median(Ts / n ** 2),
             np.mean(eps_), np.mean(eps_) / n))
    print('        T/n^2 quantiles 1/5/10/25/50/75/90/99%%: '
          + ' '.join('%.4f' % x for x in q))
    print('        P(T >= 0.02 n^2)=%.3f  P(T>=0.05 n^2)=%.3f  P(T>=0.10 n^2)=%.3f  P(T<=n)=%.3f'
          % (np.mean(Ts >= 0.02 * n ** 2), np.mean(Ts >= 0.05 * n ** 2),
             np.mean(Ts >= 0.10 * n ** 2), np.mean(Ts <= n)))
    if ratios.size:
        print('        epochs with M>=%d: %d  cost/M mean %.4f med %.4f  '
              'P(cost/M>=.05)=%.3f P(>=.1)=%.3f  mean cost %.1f (=%.3f M)'
              % (mmin, ratios.size, ratios.mean(), np.median(ratios),
                 np.mean(ratios >= .05), np.mean(ratios >= .1),
                 (ratios * Ms).mean(), ratios.mean()))
    return Ts, eps_, ratios, Ms


if __name__ == '__main__':
    ns = [int(x) for x in sys.argv[1:]] or [64, 128, 256, 512, 1024]
    nt = 2000
    for n in ns:
        run(n, nt if n <= 512 else 400, 2026)
