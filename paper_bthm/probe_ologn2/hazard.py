# -*- coding: utf-8 -*-
"""hazard.py -- one-step hazard of the top-of-scan concavity tests, plus tail shape.

For the first epoch on [0,n-1]:
  conc[m] = 1 iff A_floor(m/2) >= A_m  (test fails / scan stops here)
  E_m     = no concavity for m' in [m, n]
  cost    = n - m* + 1
Studies:
  P(E_m) for m in a grid
  hazard h_m = P(conc at m-1 | E_m)
  tail P(cost >= c n)
"""
import sys
import numpy as np

sys.path.insert(0, '.')
from q0_walk import sample_tapes


def conc_matrix(A):
    nt, n = A.shape
    S = np.zeros((nt, n + 1), dtype=np.int64)
    np.cumsum(A, axis=1, out=S[:, 1:])
    m = np.arange(2, n + 1, dtype=np.int64)
    k = m // 2
    # conc[i,j] for m=m[j] : S[i,k]*m >= S[i,m]*k
    return (S[:, k] * m[None, :] >= S[:, m] * k[None, :]), S


def main(n=512, nt=4000, seed=99):
    A = np.array(list(sample_tapes(n, nt, seed)))
    C, S = conc_matrix(A)
    ms = np.arange(2, n + 1)
    # mstar = largest m with concavity ; index into ms
    has = C.any(1)
    idx = np.where(has, (n - 2) - np.argmax(C[:, ::-1], axis=1), -1)  # reversed argmax
    # careful: column j corresponds to m = j+2 ; reversed index r -> m = n - r
    mstar = np.where(has, n - np.argmax(C[:, ::-1], axis=1), 1)
    cost = n - mstar + 1
    print('n=%d nt=%d  E[cost]/n=%.4f  P(cost=1)=%.3f' % (n, nt, cost.mean() / n, np.mean(cost == 1)))
    # E_m : no concavity at any m' >= m
    suffix_any = np.zeros((nt, n), dtype=bool)  # suffix_any[:, j] = any concavity at m >= j+2
    suffix_any[:, -1] = C[:, -1]
    for j in range(n - 3, -1, -1):
        suffix_any[:, j] = suffix_any[:, j + 1] | C[:, j]
    print('%8s %10s %10s %10s' % ('m', 'P(E_m)', 'P(conc_{m-1})', 'hazard'))
    for m in [int(x) for x in np.linspace(2, n, 17)]:
        j = m - 2
        if j <= 0:
            continue
        PE = np.mean(~suffix_any[:, j])
        pc = np.mean(C[:, j - 1])
        hz = np.mean(C[:, j - 1] & ~suffix_any[:, j]) / max(PE, 1e-12)
        print('%8d %10.4f %10.4f %10.4f' % (m, PE, pc, hz))
    print()
    for c in [0.02, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8]:
        print('P(cost >= %.2f n) = %.4f' % (c, np.mean(cost >= c * n)))
    # hazard of the *scan* as a function of remaining window
    print('\ncost/n quantiles', np.round(np.quantile(cost / n, [0.1, 0.25, 0.5, 0.75, 0.9]), 3))


if __name__ == '__main__':
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 512,
         int(sys.argv[2]) if len(sys.argv) > 2 else 4000)
