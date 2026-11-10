# -*- coding: utf-8 -*-
"""qL_meet.py -- termination (meet) statistics as a function of window length.

(a) fresh windows: P(no concavity at all) vs M, large sample.
(b) along trajectories: for each realised epoch, M_i and whether it met.
    Gives the per-epoch meet hazard h(M) = P(meet | window length M) and hence
    P(meet at a window of length >= M0) = sum over epochs with M_i >= M0.
Run:  $P qL_meet.py <n> <ntapes> <seed>
"""
import sys
import numpy as np
from q0_walk import prefix, epoch, sample_tapes, walk_b0_fast


def fresh(M, reps, rng, n_ref=1 << 16, seed=7):
    a = next(sample_tapes(n_ref, 1, seed))
    P = prefix(a)
    cnt = 0
    for _ in range(reps):
        L = int(rng.integers(0, n_ref - M + 1))
        ms, c = epoch(P, L, L + M - 1)
        if ms is None:
            cnt += 1
    return cnt / reps


def main(n, nt, seed):
    rng = np.random.default_rng(seed)
    print('(a) fresh-window meet probability vs M   (n_ref=65536, same tape, %s draws)')
    for M in (8, 12, 16, 24, 32, 48, 64, 96, 128, 192, 256, 384, 512):
        reps = 400000 if M <= 64 else 150000
        p = fresh(M, reps, rng)
        se = np.sqrt(max(p * (1 - p), 1e-12) / reps)
        print('    M=%4d  P(meet)=%.6f +- %.6f   (1/M=%.5f)' % (M, p, se, 1 / M))

    print('\n(b) trajectory meet hazard by window length band')
    band = [(3, 8), (8, 16), (16, 32), (32, 64), (64, 128), (128, 256), (256, 512), (512, 10**9)]
    tot = {b: 0 for b in band}
    met = {b: 0 for b in band}
    for a in sample_tapes(n, nt, seed):
        Pf = prefix(a)
        d, T, diag = walk_b0_fast(a, want_diag=True)
        for (L, R, M, ms, cost) in diag:
            for b in band:
                if b[0] <= M < b[1]:
                    tot[b] += 1
                    if ms is None:
                        met[b] += 1
                    break
    for b in band:
        if tot[b]:
            print('    M in [%4d,%5d)  epochs %7d  meets %5d  hazard %.6f'
                  % (b[0], min(b[1], 10**9), tot[b], met[b], met[b] / tot[b]))


if __name__ == '__main__':
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 512,
         int(sys.argv[2]) if len(sys.argv) > 2 else 600,
         int(sys.argv[3]) if len(sys.argv) > 3 else 2026)
