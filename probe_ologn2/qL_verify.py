# -*- coding: utf-8 -*-
"""qL_verify.py -- complete numerical verification of the frontier route.

Deterministic chain (Lemma B, R''=n-1):  cost_i >= cost(L_i,n-1) - delta_i,
delta_i = n-1-R_i <= i.  Hence for every N with all first N epochs realised

    T >= LB(N) := sum_{i<N} cost(i,n-1) - N(N-1)/2 ,

and since a[i..n-1] is an *unconditioned* exchangeable sub-tape of length n-i,

    E[T] >= sum_{i<N} c0 (n-i) - N(N-1)/2 ~= c0 n N - (1+c0) N^2/2
    with c0 = 0.1365 (sharp, measured) or 0.008 (rigorous, REPORT_EPOCH_BOUND).

Run:  $P qL_verify.py <n> <ntapes> <seed>
"""
import sys
import numpy as np
from q0_walk import prefix, epoch, sample_tapes, walk_b0_fast


def main(n, nt, seed):
    NMAX = int(0.5 * n)
    viol_ep = tot_ep = 0
    meet_M, nreal, Tvals = [], [], []
    SUF = []                      # (nt, NMAX) suffix costs
    for a in sample_tapes(n, nt, seed):
        Pf = prefix(a)
        d, T, diag = walk_b0_fast(a, want_diag=True)
        Tvals.append(T); nreal.append(len(diag))
        if diag and diag[-1][3] is None:
            meet_M.append(diag[-1][2])
        suf = np.empty(NMAX, dtype=np.int64)
        for i in range(NMAX):
            suf[i] = epoch(Pf, i, n - 1)[1]
        SUF.append(suf)
        for (L, R, M, ms, cost) in diag:
            i = L
            if i >= NMAX:
                continue
            tot_ep += 1
            if cost < suf[i] - (n - 1 - R):
                viol_ep += 1
    SUF = np.array(SUF)
    T = np.array(Tvals, float); nreal = np.array(nreal)
    print('n=%-5d tapes=%-5d  E[T]/n^2=%.4f  E[#ep]/n=%.3f  P(meet ever)=%.4f'
          % (n, nt, T.mean() / n ** 2, nreal.mean() / n, len(meet_M) / nt))
    print('   per-epoch frontier inequality: %d epochs, %d violations' % (tot_ep, viol_ep))
    if meet_M:
        mm = np.array(meet_M)
        print('   meet window lengths: min=%d median=%d max=%d  (n=%d)' % (mm.min(), int(np.median(mm)), mm.max(), n))
    print('   E[cost(i,n-1)]/(n-i): ', end='')
    for frac in (0.0, 0.02, 0.05, 0.1, 0.2, 0.3):
        i = int(frac * n)
        if i < NMAX:
            print('i=%.2f:%.3f ' % (frac, SUF[:, i].mean() / (n - i)), end='')
    print()

    print('\n   %6s %8s %9s %11s %11s %12s' % ('N', 'N/n', 'P(all N)', 'E[LB]/n^2', 'E[T]/n^2', 'min(T-LB)'))
    best = (0, -1)
    for alpha in (0.02, 0.05, 0.08, 0.10, 0.12, 0.15, 0.20, 0.30, 0.40):
        N = max(1, int(round(alpha * n)))
        if N > NMAX:
            continue
        ok = nreal >= N
        LB = SUF[:, :N].sum(1) - N * (N - 1) // 2
        print('   %6d %8.3f %9.4f %11.5f %11.5f %12.1f'
              % (N, N / n, ok.mean(), LB.mean() / n ** 2, T.mean() / n ** 2, (T - LB).min()))
        if LB.mean() / n ** 2 > best[1]:
            best = (N, LB.mean() / n ** 2)
    print('   best empirical LB: N=%d (%.3f n), LB/n^2=%.5f' % (best[0], best[0] / n, best[1]))
    for c0, nm in ((0.1365, 'sharp c0=0.1365'), (0.008, 'rigorous c0=0.008')):
        Ns = c0 * n / (1 + c0)
        val = c0 * (n * Ns - Ns * Ns / 2) - Ns * Ns / 2
        print('   %s: N*=%.1f (%.4f n)  predicted bound = %.4e n^2' % (nm, Ns, Ns / n, val / n ** 2))


if __name__ == '__main__':
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 512,
         int(sys.argv[2]) if len(sys.argv) > 2 else 300,
         int(sys.argv[3]) if len(sys.argv) > 3 else 2026)
