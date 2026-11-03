# -*- coding: utf-8 -*-
"""Bit-exact reproduction of the PUBLISHED BTHM-0 machine at Theta(n) reads.

The published walk (_onetable.py) uses np.corrcoef; at an exactly-zero
autocorrelation (which happens on 1-14% of steps) the sign of the test is decided
by float roundoff.  Those steps are rare per tape, so we resolve them by scanning
the window in float exactly as the published code does (O(m) reads) and keep the
O(1)-per-step integer path everywhere else.

Everything else is incremental: a cut fires only when the moving header is on the
interval boundary, so the new window is the old window minus one end cell.

Usage:
    python _bz_inc2.py <n> <tapes> [seed]
"""
import sys
import time

import numpy as np

from _bz_inc import walk_ref, make_tape, Tap, Win, r2


def walk_inc_code(T, n, approx_wnorm=False):
    """Reproduce the published machine exactly.  Integer incremental state; at an
    exact tie in the autocorrelation, evaluate the published float expression."""
    L, R = 0, n - 1
    w = Win(T, 0, n - 1)
    c = T.c
    steps = 0
    cuts = 0
    ties = 0
    tie_reads = 0
    branch3 = 0
    sync_iters = 0
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
                if d != 0:
                    b1 = d > 0
                else:
                    # exact tie: replay the published float expression on the window
                    ties += 1
                    z = np.asarray(c[w.y:w.x + 1], dtype=float)
                    T.reads += m                     # honest: the window is read
                    tie_reads += m
                    u, v = z[:-1], z[1:]
                    b1 = bool(np.corrcoef(u, v)[0, 1] > 0) if (u.std() > 1e-12 and v.std() > 1e-12) else False
        o = ("B" if b1 else "D") if b2 else ("C" if b1 else "A")
        dy, dx = {"A": (-1, 0), "B": (1, 0), "C": (0, 1), "D": (0, -1)}[o]
        ny, nx = w.y + dy, w.x + dx
        if ny < L or nx > R:
            cuts += 1
            if ny < L:
                assert w.y == L
                R = w.x - 1
                w.del_right()
            else:
                assert w.x == R
                L = w.y + 1
                w.del_left()
            sync_iters += 2
            w.sync_head()
            continue
        if ny >= nx:
            out = ny / (n - 1)
            break
        if nx < L or ny > R:
            branch3 += 1
            L, R = L + 1, R - 1
            w = Win(T, L, R)
            continue
        w.move(dy, dx)
        sync_iters += 2
    return out, T.reads, steps, cuts, ties, tie_reads, branch3, sync_iters


def main():
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 256
    NW = int(sys.argv[2]) if len(sys.argv) > 2 else 250
    seed = int(sys.argv[3]) if len(sys.argv) > 3 else 2026
    rng = np.random.default_rng(seed)
    budget = n * n / np.log(n)
    rows = []
    agree = 0
    meet_code = 0
    d_pred, d_true = [], []
    t0 = time.time()
    bad = []
    for t in range(NW):
        a = make_tape(rng, n)
        c = np.round(a * 2).astype(np.int32)
        T = Tap(c)
        d_inc, reads, steps, cuts, ties, tie_reads, branch3, synci = walk_inc_code(T, n)
        d_ref = walk_ref(a, n)
        same = ((d_ref is None and d_inc is None) or
                (d_ref is not None and d_inc is not None and abs(d_ref - d_inc) < 1e-12))
        agree += same
        if not same and len(bad) < 5:
            bad.append((t, d_ref, d_inc))
        if d_ref is not None:
            meet_code += 1
            d_true.append(d_ref)
            d_pred.append(d_inc if d_inc is not None else np.nan)
        rows.append((reads, steps, cuts, ties, tie_reads, branch3, synci))
    rows = np.array(rows, dtype=float)
    print('n=%d tapes=%d seed=%d  budget n^2/ln n=%.0f' % (n, NW, seed, budget))
    print('published-code reproduction: %d/%d tapes identical (%.4f)' % (agree, NW, agree / NW))
    if bad:
        print('  mismatches:', bad)
    print('%-12s %10s %10s %10s' % ('stat', 'mean', 'median', 'max'))
    for i, nm in enumerate(('reads', 'steps', 'cuts', 'ties', 'tie_reads', 'branch3', 'sync_iters')):
        v = rows[:, i]
        print('%-12s %10.1f %10.1f %10.0f' % (nm, v.mean(), np.median(v), v.max()))
    print('reads/n = %.2f   reads/n^2 = %.4f   reads/(n^2/ln n) = %.3f   max/budget = %.3f'
          % (rows[:, 0].mean() / n, rows[:, 0].mean() / n / n,
             rows[:, 0].mean() / budget, rows[:, 0].max() / budget))
    print('all tapes under budget: %s   (n=%d, tapes=%d)'
          % (bool((rows[:, 0] < budget).all()), n, NW))
    if d_true:
        p = np.array(d_pred); tt = np.array(d_true)
        ok = ~np.isnan(p)
        print('R2 on published-meeting tapes (dropping %d non-reproduced): %.6f'
              % ((~ok).sum(), r2(p[ok], tt[ok])) if ok.any() else 'none')
    print('elapsed %.1fs' % (time.time() - t0))


if __name__ == '__main__':
    main()
