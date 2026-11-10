# -*- coding: utf-8 -*-
"""sl2_controls -- the two traps, measured on this construction.

TRAP 1  a delta insensitive to the tape has e = O(1).  Control: `geo`, the same walk
        with the content bit removed (a fixed monotone collapse).  Its delta is the
        same on every tape -- measured spread 0 cells -- and a machine that outputs it
        without reading anything is exact, so e = 0, not Theta(n).  The ballast rule's
        delta spread is 0.36-0.40 n cells, so no such machine exists for it.
TRAP 2  a defect tested once per query fires somewhere in the walk.  Control:
        `perquery`, the same epoch branch reached whenever hash(w) % m == 0 for the
        CURRENT window (probability 1/m per query, ~ln n queries per walk).  Measured:
        the fraction of tapes on which it fires at least once, and where it first fires.
        If that fraction is ~1 there is no rare set and the rule is simply superlinear.
"""
import os
import sys

import numpy as np
from numba import njit

HERE = os.path.dirname(os.path.abspath(__file__))
PAPER = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(PAPER, 'probe_class_rerun'))
sys.path.insert(0, HERE)

import pcr_lib as P                                              # noqa: E402
import sl2_engine as E                                           # noqa: E402


@njit(cache=True)
def win_hash(a, ell, r):
    h = 0
    pw = 1
    for i in range(ell, r + 1):
        h = (h + pw * a[i]) % E.P31
        pw = (pw * 3) % E.P31
    return h


@njit(cache=True)
def geo_walk(a, Pa, Paa, PR, n):
    """geometry-only collapse: delta = 1/2 on every tape, by construction."""
    ell, r = 0, n - 1
    steps = 0
    while True:
        if ell >= r:
            return ell / (n - 1.0), steps
        steps += 1
        if steps % 2 == 0:
            ell += 1
        else:
            r -= 1


@njit(cache=True)
def perquery_first_fire(a, Pa, Paa, PR, n):
    """follow the ballast NORMAL walk; report (first window at which a per-query
    defect test would fire, m at that window); -1 if it never fires."""
    ell, r = 0, n - 1
    while True:
        if ell >= r:
            return -1, 0
        m = r - ell + 1
        if win_hash(a, ell, r) % m == 0:
            return ell, m
        if E.canon_b1(a, Pa, Paa, PR, ell, r) == 1:
            ell += 1
        else:
            r -= 1


def main():
    for n in (512, 1024, 2048):
        X = P.make_tapes(n, 200, 2026 + n)
        # trap 1 control
        d = []
        for t in range(200):
            S = P.St(X[t])
            out, steps = geo_walk(S.a, S.Pa, S.Paa, S.PR, n)
            d.append(out)
        print('n=%-5d TRAP1 control `geo`: delta takes %d distinct value(s), sd %.1f cells'
              % (n, len(set(np.round(d, 12))), float(np.std(d)) * n))
        # trap 2 control
        fire = 0
        pos = []
        for t in range(200):
            S = P.St(X[t])
            e0, m0 = perquery_first_fire(S.a, S.Pa, S.Paa, S.PR, n)
            if e0 >= 0:
                fire += 1
                pos.append(m0)
        print('n=%-5d TRAP2 control `perquery`: fires on %d/200 tapes (%.1f%%), '
              'first fire at window m: med %d (=%.2f n); hazard sum_m 1/m ~ ln n = %.2f'
              % (n, fire, 100.0 * fire / 200.0,
                 int(np.median(pos)) if pos else -1,
                 (np.median(pos) / float(n)) if pos else 0.0, np.log(n)))


if __name__ == '__main__':
    main()
