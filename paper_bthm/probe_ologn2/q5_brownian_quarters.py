# -*- coding: utf-8 -*-
"""q5: (a) Brownian-limit prediction for E[cost]/M of a fresh window;
(b) exhaustive verification of the four-quarter-block lemma.

(a) With A_j = 1 + X_j/j and D_k = X_{2k} - 2X_k, the b2=1 test at even m=2k is
D_k > 0.  Under diffusive scaling the window path converges to sigma*W, so
    Z(theta) = W(2 theta)/(2 theta) - W(theta)/theta ,  theta in (0, 1/2]
and the epoch cost satisfies cost/M -> 1/2 - theta*,  theta* = sup{theta: Z(theta) <= 0}.
This script estimates E[1/2 - theta*] by simulation (fine grid, exact BM covariance).

(b) Lemma (four quarters): if the four consecutive quarter blocks of a window
(lengths as equal as possible) have strictly increasing means, then b2 = 1 for
every m > M/2, hence the scan reaches m = ceil(M/2) and the epoch costs
>= M - ceil(M/2) + 1 = floor(M/2) + 1.
Checked here exhaustively for all 3^M tapes, M = 4..11.
"""
import itertools
import sys
import numpy as np


# ------------------------------------------------------------------ (a)
def brownian():
    rng = np.random.default_rng(3)
    N = 1024               # theta grid: i/(2N), i = 1..N  (theta <= 1/2)
    reps = 20000
    out = []
    for _ in range(reps):
        # BM on grid t_j = j/N, j = 0..N ; W(2 theta_i) with 2theta_i = i/N
        dW = rng.normal(0.0, 1.0 / np.sqrt(N), size=N)
        W = np.concatenate(([0.0], np.cumsum(dW)))     # W at j/N, j=0..N
        i = np.arange(1, N + 1)                        # i/N = 2 theta
        theta = i / (2.0 * N)
        Z = W[i] / (i / N) - W[i // 2] / (i // 2 / N)
        neg = np.flatnonzero(Z <= 0.0)
        th = theta[neg[-1]] if neg.size else 0.0
        out.append(2.0 * (0.5 - th))
    out = np.array(out)
    print('Brownian limit: E[1/2 - theta*] = %.4f  (median %.4f, P(>0)=%.3f)'
          % (out.mean(), np.median(out), np.mean(out > 0)))
    # variance of a: for the ensemble, V = Var(a) is about 0.5; D_k has step
    # variance 6V, so cost/M is scale free: 1/2 - theta* is the prediction.
    print('predicted E[cost]/M for a fresh long window: %.4f' % out.mean())


# ------------------------------------------------------------------ (b)
def quarters(M):
    """four consecutive blocks with lengths differing by <= 1"""
    q, r = divmod(M, 4)
    lens = [q + (1 if i < r else 0) for i in range(4)]
    b = [0]
    for x in lens:
        b.append(b[-1] + x)
    return lens, b


def b2_one(a, L, m):
    """1 if b2([L, L+m-1]) = 1, else 0 (exact integer arithmetic)."""
    k = max(1, m // 2)
    s = int(a[L:L + m].sum())
    sk = int(a[L:L + k].sum())
    return int(s * k > sk * m)


def check_quarters(Ms):
    for M in Ms:
        lens, b = quarters(M)
        bad = 0
        tested = 0
        for tup in itertools.product((0, 1, 2), repeat=M):
            a = np.array(tup, dtype=np.int64)
            qm = [a[b[i]:b[i + 1]].mean() for i in range(4)]
            if not (qm[0] < qm[1] < qm[2] < qm[3]):
                continue
            tested += 1
            for m in range(M // 2 + 1, M + 1):
                if b2_one(a, 0, m) != 1:
                    bad += 1
                    break
        print('M=%2d quarter lens %s : quarter-monotone tapes %6d / %d   '
              'violations of "b2=1 for all m > M/2" : %d'
              % (M, lens, tested, 3 ** M, bad))


if __name__ == '__main__':
    brownian()
    check_quarters(range(4, 12))
