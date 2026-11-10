# -*- coding: utf-8 -*-
"""Final measurement: exact BTHM-0 simulation in O(n) reads with O(log n) bits.

Machine M_exact = canonical rule with b1 = [autocorrelation > 0] in exact integer
arithmetic (tie -> False).  Pure incremental path: no step ever reads more than a
handful of cells, no rebuild after a cut.

For each n and sample of tapes this prints
  * meeting rate, and agreement with the published implementation (_onetable.py)
  * reads per tape: mean / median / p90 / max  (over meeting tapes)
  * reads/n, reads/(n^2/ln n), worst tape's share of the budget
  * the largest number of cell reads consumed by any single step
  * R^2 of the simulated output against the published code's output

Usage: python _bz_final.py <n> <tapes> [seed]
"""
import sys
import time

import numpy as np

from _bz_inc import walk_ref, walk_exact_bf, make_tape, Tap, Win, r2


def run_one(T, n, want_step_stats=False):
    L, R = 0, n - 1
    w = Win(T, 0, n - 1)
    steps = 0
    cuts = 0
    ties = 0
    max_step = 0
    while True:
        if w.y >= w.x:
            out = w.y / (n - 1)
            break
        if R - L + 1 < 3:
            out = None
            break
        steps += 1
        if steps > 20 * n:
            out = None
            break
        r0 = T.reads
        m = w.x - w.y + 1
        k = max(1, m // 2)
        b2 = (w.S * k > w.Sh * m)
        if m == 1:
            b1 = w.cy > 1
        elif m == 2:
            b1 = w.cx > w.cy
        else:
            np_ = m - 1
            A = w.S - w.cx
            B = w.S - w.cy
            C_ = w.Q - w.cx * w.cx
            D_ = w.Q - w.cy * w.cy
            if np_ * C_ - A * A <= 0 or np_ * D_ - B * B <= 0:
                b1 = False
            else:
                d = np_ * w.Rr - A * B
                if d == 0:
                    ties += 1
                b1 = d > 0
        o = ("B" if b1 else "D") if b2 else ("C" if b1 else "A")
        dy, dx = {"A": (-1, 0), "B": (1, 0), "C": (0, 1), "D": (0, -1)}[o]
        ny, nx = w.y + dy, w.x + dx
        if ny < L or nx > R:
            cuts += 1
            if ny < L:
                a_ = (w.y == L)
                R = w.x - 1
                w.del_right()
            else:
                a_ = (w.x == R)
                L = w.y + 1
                w.del_left()
            assert a_
            w.sync_head()
            if want_step_stats:
                max_step = max(max_step, T.reads - r0)
            continue
        if ny >= nx:
            out = ny / (n - 1)
            break
        if nx < L or ny > R:            # never observed; would need a rebuild
            L, R = L + 1, R - 1
            w = Win(T, L, R)
            if want_step_stats:
                max_step = max(max_step, T.reads - r0)
            continue
        w.move(dy, dx)
        if want_step_stats:
            max_step = max(max_step, T.reads - r0)
    return out, T.reads, steps, cuts, ties, max_step


def main():
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 256
    NW = int(sys.argv[2]) if len(sys.argv) > 2 else 400
    seed = int(sys.argv[3]) if len(sys.argv) > 3 else 2026
    rng = np.random.default_rng(seed)
    budget = n * n / np.log(n)
    allrows = []
    meet_rows = []
    agree = 0
    max_step = 0
    t0 = time.time()
    for t in range(NW):
        a = make_tape(rng, n)
        c = np.round(a * 2).astype(np.int32)
        T = Tap(c)
        out_inc, reads, steps, cuts, ties, ms = run_one(T, n, want_step_stats=True)
        max_step = max(max_step, ms)
        d_ref = walk_ref(a, n)
        same = ((d_ref is None and out_inc is None) or
                (d_ref is not None and out_inc is not None and abs(d_ref - out_inc) < 1e-12))
        agree += same
        allrows.append((reads, steps, cuts))
        if out_inc is not None:
            meet_rows.append((reads, steps, cuts, ties, out_inc, d_ref))
    allrows = np.array(allrows, dtype=float)
    M = np.array([(r[0], r[1], r[2], r[3]) for r in meet_rows], dtype=float)
    print('n=%4d tapes=%4d seed=%-5d meeting=%3d (%.3f)  budget n^2/ln n=%8.0f  agreement with published walk: %.4f'
          % (n, NW, seed, len(M), len(M) / NW, budget, agree / NW))
    print('   per tape (meeting set): reads mean %7.1f  median %7.1f  p90 %7.1f  max %7.0f'
          % (M[:, 0].mean(), np.median(M[:, 0]), np.percentile(M[:, 0], 90), M[:, 0].max()))
    print('   reads/n=%.2f  reads/n^2=%.4f  reads/(n^2/ln n) mean=%.3f max=%.3f  steps/n=%.3f  cuts=%.1f  ties=%.1f'
          % (M[:, 0].mean() / n, M[:, 0].mean() / n / n, M[:, 0].mean() / budget,
             M[:, 0].max() / budget, M[:, 1].mean() / n, M[:, 2].mean(), M[:, 3].mean()))
    print('   meeting tapes over budget: %d/%d   max reads in any single step: %d   all tapes: reads/n=%.2f'
          % ((M[:, 0] > budget).sum(), len(M), max_step, allrows[:, 0].mean() / n))
    both = [(r[5], r[4]) for r in meet_rows if r[5] is not None]
    if both:
        tt = np.array([b[0] for b in both])
        pp = np.array([b[1] for b in both])
        print('   R^2 of simulator vs published X on the %d tapes both meet: %.6f   (identical on %.4f)'
              % (len(both), r2(pp, tt), (np.abs(tt - pp) < 1e-12).mean()))
    print('   elapsed %.1fs' % (time.time() - t0))


if __name__ == '__main__':
    main()
