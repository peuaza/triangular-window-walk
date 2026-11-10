# -*- coding: utf-8 -*-
"""structure.py -- what do long-first-epoch tapes look like at the quarter scale?

Prints, conditional on cost_0 >= n/2 (and on other slices), the means of the
quarter sums / sqrt(n), the first-half sum, D_{n/4}, etc.
"""
import sys
import numpy as np

sys.path.insert(0, '.')
from q0_walk import sample_tapes


def stats(n=512, nt=6000, seed=3):
    A = np.array(list(sample_tapes(n, nt, seed)))
    L = n // 2
    S = np.zeros((nt, n + 1), dtype=np.int64)
    np.cumsum(A, axis=1, out=S[:, 1:])
    m = np.arange(2, n + 1)
    k = m // 2
    C = (S[:, k] * m[None, :] >= S[:, m] * k[None, :])
    has = C.any(1)
    mstar = np.where(has, n - np.argmax(C[:, ::-1], axis=1), 1)
    cost = n - mstar + 1
    Q = np.stack([S[:, n // 4], S[:, n // 2] - S[:, n // 4],
                  S[:, 3 * n // 4] - S[:, n // 2], S[:, n] - S[:, 3 * n // 4]], 1)
    D_K = S[:, n // 2] - 2 * S[:, n // 4]
    print('n=%d  E[cost]/n=%.4f  P(cost>=n/2)=%.3f' % (n, cost.mean() / n, np.mean(cost >= n / 2)))
    print('Q = quarter sums / sqrt(n):')
    for name, sel in [('all', np.ones(nt, bool)),
                      ('cost>=n/2', cost >= n / 2),
                      ('cost>=n/4', cost >= n / 4),
                      ('cost=1', cost == 1)]:
        if sel.sum() == 0:
            continue
        print('%-12s frac=%.3f  E[Q]/sqrt(n)=%s  E[x]/sqrt(n)=%+.2f  P(Q1<0)=%.2f P(Q3>0)=%.2f P(Q4>0)=%.2f'
              % (name, sel.mean(), np.round(Q[sel].mean(0) / np.sqrt(n), 2),
                 S[sel, n // 2].mean() / np.sqrt(n), np.mean(Q[sel, 0] < 0),
                 np.mean(Q[sel, 2] > 0), np.mean(Q[sel, 3] > 0)))
    # condition on the endpoint conditions only
    e1 = D_K >= 3
    e2 = S[:, n // 2] <= -3
    print('P(D_{n/4}>=3)=%.3f  P(x<=-3)=%.3f  P(both)=%.3f' % (e1.mean(), e2.mean(), (e1 & e2).mean()))
    sel = e1 & e2
    print('given both: P(cost>=n/2)=%.3f  (lift %.2f)' % (np.mean(cost[sel] >= n / 2), np.mean(cost[sel] >= n / 2) / np.mean(cost >= n / 2)))
    # where do the failures happen? first t with concavity, among sel
    idx = np.where(sel)[0]
    firstfail = []
    for i in idx[:2000]:
        js = np.where(C[i])[0]
        mm = m[js]
        # only look at m >= n/2
        mm2 = mm[mm >= n // 2]
        firstfail.append(mm2.max() if mm2.size else 0)
    firstfail = np.array(firstfail)
    print('among sel: P(no concavity m>=n/2)=%.3f' % np.mean(firstfail == 0))


if __name__ == '__main__':
    stats(int(sys.argv[1]) if len(sys.argv) > 1 else 512)
