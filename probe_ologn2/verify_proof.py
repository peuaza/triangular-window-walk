# -*- coding: utf-8 -*-
"""verify_proof.py -- numerical verification of every quantitative ingredient of
the discrete proof of  E[cost_0] >= c M.

Ingredients
  (L1)  {D_u >= 3 for all u in [M/4, M/2]}  =>  cost_0 >= M/2          (exact)
  (L2)  D_{K+s} = chord(s) + dC_{2s} - 2 dP_s                            (exact)
  (E1)  P(min(D_K, -2x) >= 3 + lam)                                     (table)
  (E2)  P(max_s |dP_s| <= lam/4)                                        (table)
  (E3)  P(max_t |dC_t| >= a | X)  vs the bridge bound B(a)=2exp(-2a^2/(N v))
  (A)   assembled  inf over X  bound:  P(G) * (1 - B(lam/2))
"""
import sys
import numpy as np

sys.path.insert(0, '.')
from q0_walk import sample_tapes


def one(M, nt, seed):
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

    # (L1) exactness: D-event implies cost
    m_ = np.arange(2, M + 1)
    k_ = m_ // 2
    Cc = (S[:, k_] * m_[None, :] >= S[:, m_] * k_[None, :])
    has = Cc.any(1)
    mstar = np.where(has, M - np.argmax(Cc[:, ::-1], axis=1), 1)
    cost = M - mstar + 1
    u = np.arange(K, 2 * K + 1)
    D = S[:, 2 * u] - 2 * S[:, u]
    surv = (D >= 3).all(1)

    # (L2) identity
    x_ = x
    chord = DK[:, None] + (s[None, :] / K) * (-2 * x_[:, None] - DK[:, None])
    lhs = S[:, 2 * (K + s)] - 2 * S[:, K + s]
    ident = int(np.abs(lhs - (chord + dC[:, ::2] - 2 * dP)).max())

    print('M=%4d K=%3d nt=%5d  E[cost]/M=%.4f  P(cost>=M/2)=%.4f  P(D-event)=%.4f  '
          'P(cost>=M/2|D-event)=%.3f  identity-err=%d'
          % (M, K, nt, cost.mean() / M, np.mean(cost >= M / 2), surv.mean(),
             np.mean(cost[surv] >= M / 2), ident))
    print('     sd(D_K)/sqrtK=%.2f  sd(-2x)/sqrtK=%.2f  median max|dP|/sqrtK=%.2f  '
          'median max|dC|/sqrtK=%.2f'
          % (DK.std() / np.sqrt(K), (-2 * x).std() / np.sqrt(K),
             np.median(mdP) / np.sqrt(K), np.median(mdC) / np.sqrt(K)))
    # (E3) bridge tail vs sharp bound (v = fraction of nonzero cells)
    v = float(np.mean(B != 0))
    print('     v=%.3f   a/sqrtK:  ' % v, end='')
    for c in [0.5, 1.0, 1.5, 2.0]:
        print('%.1f:%.3f/%.3f  ' % (c, np.mean(mdC >= c * np.sqrt(K)),
                                    2 * np.exp(-2 * (c * np.sqrt(K)) ** 2 / (2 * K * v))), end='')
    print()
    # assembled bound with the sharp bridge inequality, sup over lam = c sqrtK
    best = (0, None)
    for c in np.arange(0.25, 3.01, 0.25):
        lam = c * np.sqrt(K)
        e1 = (DK >= 3 + lam) & (-2 * x >= 3 + lam)
        e2 = mdP <= lam / 4
        G = e1 & e2
        Bp = min(1.0, 2 * np.exp(-2 * (lam / 2) ** 2 / (2 * K * v)))
        val = G.mean() * (1 - Bp)
        if val > best[0]:
            best = (val, c)
    print('     assembled:  sup_lam  P(G)(1-B) = %.3e   at c=%.2f' % best)


if __name__ == '__main__':
    for M, nt, sd in [(256, 20000, 41), (512, 20000, 42), (1024, 12000, 43), (4096, 4000, 44)]:
        one(M, nt, sd)
