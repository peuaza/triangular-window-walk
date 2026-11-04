# -*- coding: utf-8 -*-
"""qL_run.py -- is E[cost_i | F_{i-1}] >= c M_i with c uniform?  (Lemma Q, stated form)

Proxy for F_{i-1}: the run length of consecutive double cuts ending at i-1
(a double cut is m*==M, cost 1).  If the conditional mean cost/M decays along a
run, Lemma Q is false as stated and the correct statement must be amortised.

Also reports the *joint* table with M/n, because long runs force small M.
Run:  $P qL_run.py <n> <ntapes> <seed>
"""
import sys
import numpy as np
from q0_walk import prefix, epoch, sample_tapes


def trace(a):
    n = len(a)
    P = prefix(a)
    L, R = 0, n - 1
    out = []
    while True:
        if L >= R:
            break
        M = R - L + 1
        if M < 3:
            break
        mstar, cost = epoch(P, L, R)
        out.append((L, R, M, mstar, cost))
        if mstar is None:
            break
        if mstar == M:
            L, R = L + 1, R - 1
        else:
            L += 1
    return out


def main(n, nt, seed):
    runs, cost, Mv, Mn, idx = [], [], [], [], []
    nreal = []
    for a in sample_tapes(n, nt, seed):
        tr = trace(a)
        nreal.append(len(tr))
        run = 0
        for i, (L, R, M, ms, c) in enumerate(tr):
            runs.append(run); cost.append(c); Mv.append(M); Mn.append(M / n); idx.append(i / n)
            run = run + 1 if (ms is not None and ms == M) else 0
    runs = np.array(runs); cost = np.array(cost, float); Mv = np.array(Mv, float)
    Mn = np.array(Mn); idx = np.array(idx)
    print('n=%d nt=%d epochs=%d  E[cost/M]=%.4f  P(cost=1)=%.4f  E[#ep]/n=%.3f'
          % (n, nt, len(cost), (cost / Mv).mean(), (cost == 1).mean(), np.mean(nreal) / n))

    print('\n(A) run-length table (all M)')
    print('%10s %9s %11s %11s %11s' % ('run', 'epochs', 'mean cost/M', 'P(cost>=.05M)', 'P(cost>=2)'))
    edges = [0, 1, 2, 3, 4, 6, 8, 12, 16, 24, 32, 48, 64, 96, 128, 10**9]
    for lo, hi in zip(edges[:-1], edges[1:]):
        s = (runs >= lo) & (runs < hi)
        if s.sum() >= 100:
            print('%4d-%-5d %9d %11.4f %11.4f %11.4f'
                  % (lo, min(hi, 999), s.sum(), (cost[s] / Mv[s]).mean(),
                     np.mean(cost[s] >= .05 * Mv[s]), np.mean(cost[s] >= 2)))

    print('\n(B) run-length table within M/n bands (joint control)')
    for mlo, mhi in [(0.5, 0.6), (0.6, 0.7), (0.7, 0.8), (0.8, 0.9), (0.9, 1.01)]:
        sm = (Mn >= mlo) & (Mn < mhi)
        print('  M/n in [%.1f,%.1f)  epochs %d' % (mlo, mhi, sm.sum()))
        for lo, hi in [(0, 1), (1, 2), (2, 4), (4, 8), (8, 16), (16, 32), (32, 64), (64, 10**9)]:
            s = sm & (runs >= lo) & (runs < hi)
            if s.sum() >= 60:
                print('      run %3d-%-4d  n=%6d  mean cost/M %.4f  P(cost>=2) %.3f'
                      % (lo, min(hi, 999), s.sum(), (cost[s] / Mv[s]).mean(), np.mean(cost[s] >= 2)))

    print('\n(C) run-length distribution: fraction of epochs with run >= k')
    for k in (1, 2, 4, 8, 16, 32, 64, 128, 256):
        print('    k=%4d  P(run>=k)=%.4f   E[cost/M | run>=k]=%.4f'
              % (k, np.mean(runs >= k), (cost[runs >= k] / Mv[runs >= k]).mean() if (runs >= k).sum() else np.nan))


if __name__ == '__main__':
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 512,
         int(sys.argv[2]) if len(sys.argv) > 2 else 300,
         int(sys.argv[3]) if len(sys.argv) > 3 else 2026)
