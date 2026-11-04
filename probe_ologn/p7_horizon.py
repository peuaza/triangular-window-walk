# -*- coding: utf-8 -*-
"""p7: the natural 'horizon' jump for the b1=0 process -- exact or not?

Claim tested: for consecutive left cuts with fixed R, g(L+1) <= g(L).
Algorithm: remember the last g as a horizon h; the next scan starts at h
(and always re-tests r=R first, which costs O(1)).  If the claim held, the
whole run would cost O(n) b2-tests.
"""
import sys
import numpy as np
sys.path.insert(0, '..')
from _arr_longrun_fast import bit_b2_ref
from p0_tapes import enum_tapes
from p2_b0 import walk_b0


def walk_b0_horizon(a, n):
    """O(n)-read candidate; returns (delta, b2_tests)"""
    L, R = 0, n - 1
    h = R
    tests = 0
    while True:
        if R - L + 1 < 3:
            return None, tests
        # a double cut is decided by the single test at r=R
        tests += 1
        if bit_b2_ref(a[L:R + 1]) == 0:
            L, R = L + 1, R - 1
            h = R
            continue
        # no double cut: scan down from the horizon
        g = None
        start = min(h, R - 1)
        for r in range(start, L, -1):
            tests += 1
            if bit_b2_ref(a[L:r + 1]) == 0:
                g = r
                break
        if g is None:
            return L / (n - 1), tests
        L += 1
        h = g


def main():
    for n in (int(x) for x in sys.argv[1:]):
        tapes = list(enum_tapes(n))
        bad = 0
        tot_tests = 0
        ex = None
        for a in tapes:
            d0, _ = walk_b0(a, n)
            d1, t = walk_b0_horizon(a, n)
            tot_tests += t
            same = (d0 is None and d1 is None) or (d0 is not None and d1 is not None and abs(d0 - d1) < 1e-12)
            if not same:
                bad += 1
                if ex is None:
                    ex = (a.tolist(), d0, d1)
        print('n=%d tapes=%d  horizon-jump mismatches: %d (%.3f)   b2-tests/tape=%.2f  (n=%.2f)'
              % (n, len(tapes), bad, bad / len(tapes), tot_tests / len(tapes), n))
        if ex:
            print('   first counterexample:', ex)


if __name__ == '__main__':
    main()
