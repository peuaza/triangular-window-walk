# -*- coding: utf-8 -*-
"""qL_frontier.py -- test the frontier route to E[T] = Omega(n^2).

Deterministic chain (Lemma B with R''=n-1), for every realised epoch i:
      cost_i  >=  cost(i, n-1) - delta_i ,      delta_i = n-1-R_i <= i
so with N epochs realised
      T  >=  sum_{i<N} cost(i, n-1)  -  sum_{i<N} i .

cost(i,n-1) is an *unconditioned* scan of the suffix a[i..n-1] (length n-i), i.e.
exactly a fresh window, so E[cost(i,n-1)] >= c0 (n-i) with c0 the single-epoch
constant (0.008 rigorous / 0.13 sharp).  Hence

      E[T] >= c0 * sum_{i<N}(n-i) - N^2/2  ~=  c0*n*N - (1+c0)N^2/2 ,
maximised at N* = c0 n/(1+c0), value c0^2 n^2/(2(1+c0)).

This script checks the chain pathwise and measures everything.
"""
import sys
import numpy as np
from q0_walk import prefix, epoch, sample_tapes


def run(n, nt, seed):
    rng = np.random.default_rng(seed)
    rows = []
    q_meet = []          # did the walk meet (terminate early)?
    for a in sample_tapes(n, nt, seed):
        P = prefix(a)
        d, T, diag = __import__('q0_walk').walk_b0_fast(a, want_diag=True)
        realized = diag
        # did it end by meeting?  (last record has mstar None)
        met = bool(realized and realized[-1][3] is None)
        nreal = len(realized)
        # suffix costs  cost(i, n-1) for i = 0..n-3
        suf = np.full(n, np.nan)
        for i in range(0, n - 2):
            _, c = epoch(P, i, n - 1)
            suf[i] = c
        rows.append((T, nreal, met, suf, realized))
    return rows


def main(n, nt, seed):
    rows = run(n, nt, seed)
    T = np.array([r[0] for r in rows], float)
    nreal = np.array([r[1] for r in rows], float)
    met = np.array([r[2] for r in rows], bool)
    print('n=%d tapes=%d  E[T]/n^2=%.4f  E[#epochs]/n=%.3f  P(meet ever)=%.4f (%d tapes)'
          % (n, len(T), T.mean() / n**2, nreal.mean() / n, met.mean(), met.sum()))
    print('   E[cost(i,n-1)]/(n-i) for i = 0, n/100, n/50, n/20, n/10:')
    for frac in (0.0, 0.01, 0.02, 0.05, 0.1):
        i = int(frac * n)
        v = np.array([r[3][i] for r in rows], float)
        print('      i=%6d  mean cost=%.2f  /(n-i)=%.4f' % (i, np.nanmean(v), np.nanmean(v) / (n - i)))

    print('\n   chain check:  LB(N) = sum_{i<N} cost(i,n-1) - N(N-1)/2   vs  T')
    print('   %6s %10s %10s %10s %10s %10s' % ('N', 'P(all N ep)', 'E[LB]', 'E[T]', 'min(T-LB)', 'E[LB]/n^2'))
    for alpha in (0.005, 0.008, 0.02, 0.05, 0.1, 0.2, 0.4):
        N = max(1, int(round(alpha * n)))
        ok = nreal >= N
        if ok.sum() == 0:
            print('   %6d  -- no tape reaches N epochs' % N)
            continue
        LB = np.array([r[3][:N].sum() - N * (N - 1) / 2 for r in rows], float)
        diff = T - LB
        print('   %6d %10.4f %10.1f %10.1f %10.1f %10.5f'
              % (N, ok.mean(), LB[ok].mean(), T[ok].mean(), diff[ok].min(), LB[ok].mean() / n**2))
        print('            (violations T<LB: %d of %d)' % (np.sum(diff < -1e-9), len(diff)))

    # the two optima
    for c0, name in ((0.008, 'rigorous c0=0.008'), (0.13, 'sharp c0=0.13')):
        Ns = c0 * n / (1 + c0)
        val = c0 * (n * Ns - Ns * Ns / 2) - Ns * Ns / 2
        print('\n   %s: N*=%.1f (alpha=%.4f)  bound = %.3e n^2' % (name, Ns, Ns / n, val / n**2))


if __name__ == '__main__':
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 512,
         int(sys.argv[2]) if len(sys.argv) > 2 else 400,
         int(sys.argv[3]) if len(sys.argv) > 3 else 2026)
