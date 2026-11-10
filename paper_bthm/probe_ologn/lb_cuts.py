# -*- coding: utf-8 -*-
"""Route (d) probe: aggregate route (a) over ALL cuts and charge the jumps.

(1) per-cut one-way CC = log2 #distinct rows (exact, every k, both directions)
(2) sum_k CC_k  -- the quantity in the jump-charging theorem
        T * Jbar * s >= sum_k CC_k ,   Jbar = mean_t |sigma(t+1)-sigma(t)|
(3) two-way lower bounds at the middle cut: log2 rank of a full
    (count-balanced) sub-rectangle and a greedy fooling set.  These are the only
    lower-bound tools available for CC_2way; the table shows how far they are
    from the min(k,n-k)*log2(3) upper bound, i.e. what route (d) still needs.

Deltas come from cert_n<n>.npy (row 0 = meeting position L, i.e. delta = L/(n-1),
or -1 for no halt), in p0_tapes.enum_tapes order.

Usage:  python lb_cuts.py 6 8 10 11 12 13
"""
import sys
import random
from collections import defaultdict

import numpy as np

sys.path.insert(0, '.')
from p0_tapes import enum_tapes, count_tapes   # noqa: E402

LOG3 = float(np.log2(3.0))


def load(n):
    A = np.array([a for a in enum_tapes(n)], dtype=np.int8)
    D = np.load('cert_n%d.npy' % n)[0].astype(np.int16)
    assert len(A) == len(D) == count_tapes(n)
    return A, D


def cut_profile(A, D, n):
    out = []
    for k in range(1, n):
        rows, cols = defaultdict(dict), defaultdict(dict)
        for i in range(len(A)):
            lh = tuple(int(v) for v in A[i, :k])
            rh = tuple(int(v) for v in A[i, k:])
            d = int(D[i])
            rows[lh][rh] = d
            cols[rh][lh] = d
        sig_r = {tuple(sorted(v.items())) for v in rows.values()}
        sig_c = {tuple(sorted(v.items())) for v in cols.values()}
        out.append((k, len(rows), len(sig_r), len(cols), len(sig_c)))
    return out


def bits(x):
    return float(np.log2(x)) if x and x > 1 else 0.0


def balanced_rect_rank(A, D, n):
    """log2 rank_Reals of the biggest full count-balanced rectangle at k=n//2.

    Rows = left halves with a fixed (c0,c1,c2), columns = right halves with a
    fixed (d0,d1,d2) such that c+d has #0-#2 = eps.  Every pair in such a block
    is a tape of the ensemble, so the block lies inside the promise and
    CC_2way >= log2 rank(block).
    """
    k = n // 2
    eps = 0 if n % 2 == 0 else 1
    buck = defaultdict(lambda: defaultdict(dict))
    for i in range(len(A)):
        L = A[i, :k].astype(np.int64)
        R = A[i, k:].astype(np.int64)
        lc = tuple(int(v) for v in np.bincount(L, minlength=3))
        rc = tuple(int(v) for v in np.bincount(R, minlength=3))
        buck[(lc, rc)][tuple(int(v) for v in L)][tuple(int(v) for v in R)] = int(D[i])
    best = (0, None, None)
    for (lc, rc), d in buck.items():
        if lc[0] + rc[0] - lc[2] - rc[2] != eps:
            continue
        xs = sorted(d)
        ys = sorted({y for dd in d.values() for y in dd})
        ys = [y for y in ys if all(y in d[x] for x in xs)]
        xs = [x for x in xs if all(y in d[x] for y in ys)]
        if len(xs) < 2 or len(ys) < 2 or len(xs) * len(ys) > 3_000_000:
            continue
        M = np.array([[d[x][y] for y in ys] for x in xs], dtype=np.float64)
        r = int(np.linalg.matrix_rank(M))
        if r > best[0]:
            best = (r, (lc, rc), (len(xs), len(ys)))
    return best


def fooling_set(A, D, n, value=None, trials=8, cap=150000):
    """greedy fooling set at the middle cut: all diagonal value v, cross != v.

    The cross pair (x_i,y_j) must itself be a tape of the ensemble, otherwise
    the fooling-set argument gives nothing (the protocol need not be correct off
    the promise).  (x_i,y_j) is a tape iff #0-#2 of the two halves sums to eps,
    i.e. iff the two rows have the same d = #0-#2 value -- so the whole fooling
    set is restricted to one d-class.
    """
    k = n // 2
    eps = 0 if n % 2 == 0 else 1
    idx = defaultdict(dict)
    dval = {}
    for i in range(len(A)):
        x = tuple(int(v) for v in A[i, :k])
        y = tuple(int(v) for v in A[i, k:])
        idx[x][y] = int(D[i])
        dval[x] = sum(1 for v in x if v == 0) - sum(1 for v in x if v == 2)
    if value is None:
        cnt = defaultdict(int)
        for i in range(len(D)):
            cnt[int(D[i])] += 1
        value = max(cnt.items(), key=lambda kv: kv[1])[0]
    pairs = [(x, y) for x, dd in idx.items() for y, v in dd.items() if v == value]
    best, bestd = 0, None
    for d in sorted(set(dval.values())):
        sub = [p for p in pairs if dval[p[0]] == d]
        if not sub:
            continue
        rng = random.Random(7)
        for _ in range(trials):
            rng.shuffle(sub)
            chosen = []
            for (x, y) in sub[:min(len(sub), cap)]:
                if all(idx[x].get(y2, -99) != value and idx[x2].get(y, -99) != value
                       for (x2, y2) in chosen):
                    chosen.append((x, y))
            if len(chosen) > best:
                best, bestd = len(chosen), d
    return value, best, len(pairs), bestd


if __name__ == '__main__':
    for n in [int(x) for x in sys.argv[1:]] or [6, 8, 10]:
        A, D = load(n)
        prof = cut_profile(A, D, n)
        tot = 0.0
        prob = 0.0
        print('n=%2d tapes=%6d  distinct deltas=%d' % (n, len(A), len({int(v) for v in D})))
        print('   k   #lh  #rows  cc1(L->R)   #rh  #cols  cc1(R->L)   min*log2(3)')
        for (k, nl, nr, nr2, nc) in prof:
            a, b = bits(nr), bits(nc)
            tot += max(a, b)
            prob += min(a, b)
            print('  %2d  %4d  %5d    %6.2f    %4d  %5d    %6.2f      %6.2f'
                  % (k, nl, nr, a, nr2, nc, b, min(k, n - k) * LOG3))
        cap = 0.25 * n * n * LOG3
        print('   sum_k max(dir) = %.2f  |  sum_k min(dir) = %.2f (= cc2 upper bd)'
              % (tot, prob))
        print('   ceiling 0.25*n^2*log2(3) = %.2f   ratio(min) = %.3f  ratio(max) = %.3f'
              % (cap, prob / cap, tot / cap))
        print('   sum_k min(k,n-k)*log2(3) = %.2f'
              % (LOG3 * sum(min(k, n - k) for k in range(1, n))))
        r, cls, shp = balanced_rect_rank(A, D, n)
        print('   cc2 LB rank : block %s shape %s -> rank %d -> %.2f bits'
              % (cls, shp, r, bits(r)))
        v, fs, npair, fsl = fooling_set(A, D, n)
        print('   cc2 LB fool : value %d, d-class %s, size %d / %d pairs -> %.2f bits'
              % (v, fsl, fs, npair, bits(fs)))
        sys.stdout.flush()
