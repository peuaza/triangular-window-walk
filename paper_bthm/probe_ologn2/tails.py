# -*- coding: utf-8 -*-
"""tails.py -- check the two concentration inputs of the proof numerically.

Input 1 (endpoint):  P(min(D_K, -2x) >= 3 + lambda) for lambda = c sqrt(K).
Input 2 (bridges):   distribution of max_s |dP_s| and max_t |dC_t|,
                     compared with the Serfling/union-bound form
                         B(a) = (N/(2h)) * 2 exp(-a'^2/(N/2)),  a'=a-2h, h=0.2 sqrt(N)
                     and with the sharp Brownian-bridge form 2 exp(-2a^2/(N v)).
"""
import sys
import numpy as np

sys.path.insert(0, '.')
from q0_walk import sample_tapes


def run(M=1024, nt=6000, seed=23):
    K = M // 4
    A = np.array(list(sample_tapes(M, nt, seed)))
    B = A - 1
    S = np.zeros((nt, M + 1), dtype=np.int64)
    np.cumsum(B, axis=1, out=S[:, 1:])
    x = S[:, 2 * K]
    DK = S[:, 2 * K] - 2 * S[:, K]
    p = S[:, K:2 * K + 1] - S[:, K:K + 1]
    pK = p[:, -1]
    s = np.arange(K + 1)
    dP = p - (s[None, :] / K) * pK[:, None]
    q = S[:, 2 * K:] - S[:, 2 * K:2 * K + 1]
    R = q[:, -1]
    t = np.arange(2 * K + 1)
    dC = q - (t[None, :] / (2 * K)) * R[:, None]
    mdP = np.abs(dP).max(1)
    mdC = np.abs(dC).max(1)
    print('M=%d K=%d  sqrtK=%.1f' % (M, K, np.sqrt(K)))
    print('  quantiles max|dP|/sqrtK:', np.round(np.quantile(mdP / np.sqrt(K), [.5, .75, .9, .95, .99]), 2))
    print('  quantiles max|dC|/sqrtK:', np.round(np.quantile(mdC / np.sqrt(K), [.5, .75, .9, .95, .99]), 2))
    print('  sd(D_K)/sqrtK = %.2f   sd(-2x)/sqrtK = %.2f' % (DK.std() / np.sqrt(K), (-2 * x).std() / np.sqrt(K)))
    for c in [0, 0.5, 1, 1.5, 2, 2.5, 3]:
        lam = c * np.sqrt(K)
        e = (DK >= 3 + lam) & (-2 * x >= 3 + lam)
        print('   c=%.1f  P(both endpoints >= 3+lam)=%.5f' % (c, e.mean()))
    # empirical bridge tails vs theory
    print('  a/sqrtK   P(max|dC|>=a)   serfling-grid-bound   sharp-bridge-form')
    N = 2 * K
    for c in [1.0, 1.5, 2.0, 2.5, 3.0]:
        a = c * np.sqrt(K)
        h = 0.2 * np.sqrt(N)
        ap = max(a - 2 * h, 0.01)
        grid = (N / (2 * h)) * 2 * np.exp(-2 * ap ** 2 / N)
        sharp = 2 * np.exp(-2 * a ** 2 / (N * 0.5))
        print('   %.1f      %.4f            %.4g                %.4g'
              % (c, np.mean(mdC >= a), min(grid, 1), min(sharp, 1)))
    # also: does the assembled bound work?  P(success) >= P(G)*(1-tail)
    for c in [2.0, 2.5, 3.0]:
        lam = c * np.sqrt(K)
        e1 = (DK >= 3 + lam) & (-2 * x >= 3 + lam)
        e2 = mdP <= lam / 2
        G = e1 & e2
        tail = np.mean(mdC[G] >= lam / 2) if G.sum() > 20 else float('nan')
        print('   assembled c=%.1f: P(G)=%.5f  nG=%4d  P(maxdC>=lam/2|G)=%.3f  bound=%.2e'
              % (c, G.mean(), G.sum(), tail, G.mean() * (1 - min(1, (N / (0.4 * np.sqrt(N))) * 2 * np.exp(-2 * (lam / 2 - 0.4 * np.sqrt(N)) ** 2 / N)))))


if __name__ == '__main__':
    run(int(sys.argv[1]) if len(sys.argv) > 1 else 1024,
        int(sys.argv[2]) if len(sys.argv) > 2 else 6000)
