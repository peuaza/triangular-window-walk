# -*- coding: utf-8 -*-
"""decomp.py -- verify the four-block decomposition used by the proof attempt.

M = 4K.  b = a-1, S_m = sum_{i<m} b_i,  D_u = S_{2u} - 2S_u.
Q1=[0,K), Q2=[K,2K), C=[2K,4K).
For u = K+s (0<=s<=K):
    D_{K+s} = chord(s) + dC_{2s} - 2 dP_s
with  chord(s) linear from D_K to D_{2K}=-2x,
      dP_s = p_s - (s/K) p_K     (Q2 walk deviation from its chord),
      dC_t = q_t - (t/2K) R      (C walk deviation from its chord, R=sum C=-x).
Checks the identity, the implication  {D_u>=3 for u in [K,2K]} => cost>=M/2,
and the numerical sizes of the endpoint margins and bridge deviations.
"""
import sys
import numpy as np

sys.path.insert(0, '.')
from q0_walk import sample_tapes


def analyse(M=512, nt=4000, seed=5):
    K = M // 4
    A = np.array(list(sample_tapes(M, nt, seed)))
    B = A - 1
    S = np.zeros((nt, M + 1), dtype=np.int64)
    np.cumsum(B, axis=1, out=S[:, 1:])

    # full scan, cost
    m = np.arange(2, M + 1)
    k = m // 2
    Cc = (S[:, k] * m[None, :] >= S[:, m] * k[None, :])
    has = Cc.any(1)
    mstar = np.where(has, M - np.argmax(Cc[:, ::-1], axis=1), 1)
    cost = M - mstar + 1

    # D_u for u in [K,2K]
    u = np.arange(K, 2 * K + 1)
    D = S[:, 2 * u] - 2 * S[:, u]
    surv = (D >= 3).all(1)
    print('M=%d K=%d nt=%d' % (M, K, nt))
    print('E[cost]/M = %.4f   P(cost>=M/2)=%.4f   P(D_u>=3 all u)=%.4f' %
          (cost.mean() / M, np.mean(cost >= M / 2), surv.mean()))
    print('P(cost>=M/2 | D-event) = %.4f   P(D-event | cost>=M/2) = %.4f' %
          (np.mean(cost[surv] >= M / 2), np.mean(surv[cost >= M / 2])))
    print('min over u of D_u: mean %.2f  (sqrt K=%.2f)' % (D.min(1).mean(), np.sqrt(K)))

    # identity check
    x = S[:, 2 * K]
    DK = S[:, 2 * K] - 2 * S[:, K]
    p = S[:, K:K + K + 1] - S[:, K:K + 1]          # p_s, s=0..K
    pK = p[:, -1]
    q = S[:, 2 * K:4 * K + 1] - S[:, 2 * K:2 * K + 1]  # q_t, t=0..2K
    R = q[:, -1]
    s = np.arange(K + 1)
    chord = DK[:, None] + (s[None, :] / K) * (-2 * x[:, None] - DK[:, None])
    dP = p - (s[None, :] / K) * pK[:, None]
    t = np.arange(2 * K + 1)
    dC = q - (t[None, :] / (2 * K)) * R[:, None]
    lhs = S[:, 2 * (K + s)] - 2 * S[:, K + s]
    rhs = chord + dC[:, ::2] - 2 * dP
    print('identity  max|lhs-rhs| =', int(np.abs(lhs - rhs).max()))

    # endpoint margins
    print('D_K/sqrt(K): mean %.2f sd %.2f ; -2x/sqrt(K): mean %.2f sd %.2f'
          % (DK.mean() / np.sqrt(K), DK.std() / np.sqrt(K),
             (-2 * x).mean() / np.sqrt(K), (-2 * x).std() / np.sqrt(K)))
    print('P(D_K>=3)=%.3f  P(-2x>=3)=%.3f  P(both)=%.3f'
          % (np.mean(DK >= 3), np.mean(-2 * x >= 3), np.mean((DK >= 3) & (-2 * x >= 3))))
    for c in [0.5, 1.0, 1.5, 2.0]:
        lam = c * np.sqrt(K)
        g1 = (DK >= 3 + lam) & (-2 * x >= 3 + lam)
        g2 = np.abs(dP).max(1) <= lam / 4
        print('  c=%.1f lam=%.1f  P(endpoints)=%.3f  P(max|dP|<=lam/4)=%.3f  P(G)=%.3f'
              % (c, lam, g1.mean(), g2.mean(), (g1 & g2).mean()))
    for c in [0.5, 1.0, 1.5, 2.0]:
        lam = c * np.sqrt(K)
        print('  c=%.1f  P(max|dC|>=lam/2)=%.3f   P(max|dC|>=lam/2 | G)=%.3f'
              % (c, np.mean(np.abs(dC).max(1) >= lam / 2),
                 np.mean(np.abs(dC).max(1)[(DK >= 3 + lam) & (-2 * x >= 3 + lam)] >= lam / 2)))
    print('typical sizes: sd(dP max) %.2f sqrtK, sd(dC max) %.2f sqrtK'
          % (np.abs(dP).max(1).std() / np.sqrt(K), np.abs(dC).max(1).std() / np.sqrt(K)))


if __name__ == '__main__':
    analyse(int(sys.argv[1]) if len(sys.argv) > 1 else 512,
            int(sys.argv[2]) if len(sys.argv) > 2 else 4000)
