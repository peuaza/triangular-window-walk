# -*- coding: utf-8 -*-
"""q7: three checks that pin down where the quadratic cost lives.

(1) Exact identity check: for every window and every m,
       b2 = 1 (even m=2k)  <=>  D_k > 0,   D_k = X_{2k} - 2 X_k
       b2 = 1 (odd  m=2k+1)<=>  k*(X_{2k+1} - 2 X_k) > X_k
    (X = prefix sums of b = a-1).  Verified on all 3^n tapes, n <= 9.

(2) Frontier inequality:  cost(L,R) >= cost(L,n-1) - (n-1-R).
    (the suffix scan shares every test at m <= R-L+1)  Verified exhaustively.

(3) Annealed vs quenched: along the trajectory, bucket epochs by M and compare
    mean cost/M with the fresh-window value for the same M; also compare the
    realised T with Sum_L chat(M_L) where chat is the fresh-window mean.
"""
import itertools
import sys
import numpy as np
from q0_walk import prefix, epoch, walk_b0_fast, sample_tapes


def check_identities(nmax=9):
    bad = 0
    tot = 0
    for n in range(3, nmax + 1):
        for tup in itertools.product((0, 1, 2), repeat=n):
            a = np.array(tup, dtype=np.int64)
            b = a - 1
            X = np.concatenate(([0], np.cumsum(b)))
            for m in range(2, n + 1):
                k = m // 2
                s = int(a[:m].sum())
                sk = int(a[:k].sum())
                b2 = int(s * k > sk * m)
                if m % 2 == 0:
                    pred = int(X[m] - 2 * X[k] > 0)
                else:
                    pred = int(k * (X[m] - 2 * X[k]) > X[k])
                tot += 1
                if pred != b2:
                    bad += 1
    print('(1) identity check: %d (tape,m) pairs, mismatches %d' % (tot, bad))


def check_frontier(ns=(10, 11, 12)):
    tot = bad = 0
    for n in ns:
        for a in sample_tapes(n, 400, 5):
            P = prefix(a)
            for L in range(0, n - 2):
                N = n - L
                ms, cs = epoch(P, L, n - 1)
                for R in range(L + 2, n):
                    M = R - L + 1
                    _, c = epoch(P, L, R)
                    tot += 1
                    if c < cs - (N - M):
                        bad += 1
    print('(2) frontier inequality cost(L,R) >= cost(L,n-1)-(n-1-R): %d pairs, violations %d'
          % (tot, bad))


def annealed_vs_quenched(n, ntapes, seed=2026):
    # fresh-window reference curve chat(M)
    n_ref = 8192
    aref = next(sample_tapes(n_ref, 1, 7))
    Pref = prefix(aref)
    rng = np.random.default_rng(5)
    Ms = np.array([16, 24, 32, 48, 64, 96, 128, 192, 256, 384, 512, 768, 1024, 1536, 2048, 4096])
    chat = {}
    for M in Ms:
        reps = max(150, int(4e5 / M))
        vals = []
        for _ in range(reps):
            L = int(rng.integers(0, n_ref - M + 1))
            _, c = epoch(Pref, L, L + M - 1)
            vals.append(c)
        chat[int(M)] = float(np.mean(vals))
    print('\n(3) fresh-window chat(M) = E[cost|M]:  '
          + '  '.join('%d:%.1f' % (M, chat[int(M)]) for M in Ms))
    print('    chat(M)/M:  ' + '  '.join('%d:%.3f' % (M, chat[int(M)] / M) for M in Ms))

    buckets = [(32, 64), (64, 128), (128, 256), (256, 512), (512, 1 << 30)]
    acc = {b: [0, 0.0] for b in buckets}          # count, sum cost
    Tot, Pred = [], []
    for a in sample_tapes(n, ntapes, seed):
        P = prefix(a)
        d, T, diag = walk_b0_fast(a, want_diag=True)
        pred = 0.0
        for (L, R, M, ms, cost) in diag:
            for b in buckets:
                if b[0] <= M < b[1]:
                    acc[b][0] += 1
                    acc[b][1] += cost
                    break
            # nearest reference M
            Mref = min(Ms, key=lambda x: abs(x - M))
            pred += chat[int(Mref)] * (M / Mref) if Mref else 0
        Tot.append(T)
        Pred.append(pred)
    print('    quenched (trajectory) mean cost/M by M-bucket vs fresh window:')
    for b in buckets:
        c, s = acc[b]
        if c:
            mid = (b[0] + min(b[1], n)) / 2
            print('      M in [%4d,%4d): epochs %6d  mean cost/M %.4f' % (b[0], b[1], c, s / c / mid))
    Tot = np.array(Tot, float)
    Pred = np.array(Pred, float)
    print('    per-tape: mean realised T = %.0f, mean predicted from chat(M_L) = %.0f, ratio %.3f'
          % (Tot.mean(), Pred.mean(), Tot.mean() / Pred.mean()))


if __name__ == '__main__':
    check_identities(8)
    check_frontier()
    annealed_vs_quenched(512, 250)
