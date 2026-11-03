# -*- coding: utf-8 -*-
"""Read cost of exact BTHM-0 simulation, on the meeting set, across n.

Reports, for each n: meeting rate, reads per tape (mean/median/p90/max) over tapes
that produce output, reads/n, reads/(n^2/ln n), the worst tape's fraction of the
budget, and whether the incremental machine reproduces the published walk exactly.

Usage:  python _bz_inc3.py <n> <tapes> [seed]
"""
import sys
import time

import numpy as np

from _bz_inc import walk_ref, walk_exact_bf, make_tape, Tap, Win, r2
from _bz_inc2 import walk_inc_code


def main():
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 256
    NW = int(sys.argv[2]) if len(sys.argv) > 2 else 400
    seed = int(sys.argv[3]) if len(sys.argv) > 3 else 2026
    rng = np.random.default_rng(seed)
    budget = n * n / np.log(n)
    rec = []
    agree = 0
    t0 = time.time()
    for t in range(NW):
        a = make_tape(rng, n)
        c = np.round(a * 2).astype(np.int32)
        T = Tap(c)
        d_inc, reads, steps, cuts, ties, tie_reads, br3, synci = walk_inc_code(T, n)
        d_ref = walk_ref(a, n)
        same = ((d_ref is None and d_inc is None) or
                (d_ref is not None and d_inc is not None and abs(d_ref - d_inc) < 1e-12))
        agree += same
        rec.append((d_ref is not None, reads, steps, cuts, ties, tie_reads, d_ref))
    rec = np.array([(r[0], r[1], r[2], r[3], r[4], r[5]) for r in rec], dtype=float)
    mt = rec[:, 0] > 0
    M = rec[mt]
    print('n=%4d  tapes=%4d  seed=%d   meeting=%d (%.3f)   budget n^2/ln n = %.0f   reproduces published walk: %.4f'
          % (n, NW, seed, mt.sum(), mt.mean(), budget, agree / NW))
    print('   reads: mean %8.1f  median %8.1f  p90 %8.1f  max %8.0f'
          % (M[:, 1].mean(), np.median(M[:, 1]), np.percentile(M[:, 1], 90), M[:, 1].max()))
    print('   reads/n = %.2f   reads/n^2 = %.4f   reads/(n^2/ln n): mean %.3f  max %.3f   steps/n = %.3f   cuts/tape = %.1f   ties/tape = %.1f'
          % (M[:, 1].mean() / n, M[:, 1].mean() / n / n, M[:, 1].mean() / budget,
             M[:, 1].max() / budget, M[:, 2].mean() / n, M[:, 3].mean(), M[:, 4].mean()))
    print('   tapes over budget: %d / %d' % ((M[:, 1] > budget).sum(), len(M)))
    print('   elapsed %.1fs' % (time.time() - t0))


if __name__ == '__main__':
    main()
