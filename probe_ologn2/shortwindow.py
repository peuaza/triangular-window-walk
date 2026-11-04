# -*- coding: utf-8 -*-
"""shortwindow.py -- the fix: use a window of length 2J near the top, J = delta*M.

Lemma 1 is local in u, so {D_u >= 3 for u in [M/2-J, M/2]} already gives
cost_0 >= 2J.  In the four-block decomposition with J in place of K the two
bridges have length J and 2J (wiggle sqrt J) while the chord endpoints
D_{M/2-J}, -2x both have sd Theta(sqrt M).  For J = delta*M with delta small
the margins therefore dominate the bridge wiggles by a factor delta^{-1/2},
and the union bound over the two bridge events is no longer lossy.

Prints:
  P(D_u>=3 on the window)                    (the true sufficient event)
  P(min(D_{u0},-2x) >= 3+c sqrt J)           (endpoint, chord)
  P(max|dP| <= c sqrt J /4),  P(max|dC| <= c sqrt J /2)
  assembled bound  P(E1) - P(dP bad) - P(dC bad),  times 2J/M
"""
import sys
import numpy as np

sys.path.insert(0, '.')
from q0_walk import sample_tapes


def run(M, deltas, nt, seed):
    A = np.array(list(sample_tapes(M, nt, seed)))
    B = A - 1
    S = np.zeros((nt, M + 1), dtype=np.int64)
    np.cumsum(B, axis=1, out=S[:, 1:])
    m_ = np.arange(2, M + 1)
    k_ = m_ // 2
    Cc = (S[:, k_] * m_[None, :] >= S[:, m_] * k_[None, :])
    has = Cc.any(1)
    mstar = np.where(has, M - np.argmax(Cc[:, ::-1], axis=1), 1)
    cost = M - mstar + 1
    v = float(np.mean(B != 0))
    print('M=%d  v=%.3f  E[cost]/M=%.4f  P(cost>=M/2)=%.4f' % (M, v, cost.mean() / M, np.mean(cost >= M / 2)))
    print('%7s %7s %8s %10s %10s %10s %10s %10s %10s'
          % ('delta', 'J', 'P(Dwin)', 'c*', 'P(E1)', 'P(dP ok)', 'P(dC ok)', 'P(G)', '2J/M*P(G)'))
    for delta in deltas:
        J = int(round(delta * M))
        if J < 4:
            continue
        u0 = M // 2 - J
        u = np.arange(u0, M // 2 + 1)
        D = S[:, 2 * u] - 2 * S[:, u]
        surv = (D >= 3).all(1)
        # blocks
        p = S[:, u0:M // 2 + 1] - S[:, u0:u0 + 1]
        pJ = p[:, -1]
        s = np.arange(J + 1)
        dP = p - (s[None, :] / J) * pJ[:, None]
        q = S[:, M - 2 * J:] - S[:, M - 2 * J:M - 2 * J + 1]
        t = np.arange(2 * J + 1)
        dC = q - (t[None, :] / (2 * J)) * q[:, -1][:, None]
        mdP = np.abs(dP).max(1)
        mdC = np.abs(dC).max(1)
        Du0 = S[:, 2 * u0] - 2 * S[:, u0]
        x = S[:, M // 2]
        best = (0, 0, 0, 0, 0)
        for c in np.arange(0.5, 6.01, 0.5):
            lam = c * np.sqrt(J)
            e1 = (Du0 >= 3 + lam) & (-2 * x >= 3 + lam)
            okP = np.mean(mdP <= lam / 4)
            okC = np.mean(mdC <= lam / 2)
            val = max(0.0, e1.mean() - (1 - okP) - (1 - okC))
            if val > best[0]:
                best = (val, c, e1.mean(), okP, okC)
        val, c, pe1, okp, okc = best
        print('%7.3f %7d %8.4f %10.1f %10.4f %10.4f %10.4f %10.4f %10.4f'
              % (delta, J, surv.mean(), c, pe1, okp, okc, val, 2 * J / M * val))


if __name__ == '__main__':
    for M, nt, sd in [(1024, 8000, 61), (4096, 3000, 62), (16384, 600, 63)]:
        run(M, [0.005, 0.01, 0.02, 0.04, 0.08, 0.125], nt, sd)
        print()
