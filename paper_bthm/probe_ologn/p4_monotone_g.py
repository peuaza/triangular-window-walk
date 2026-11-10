# -*- coding: utf-8 -*-
"""p4: is g(L) = max{r<=R : b2([L,r])=0} monotone in L while R is fixed?

If g(L+1) <= g(L), the largest-r scan for L+1 never needs to look above g(L),
so a single downward sweep costs O(n) total.  Test exhaustively.
"""
import sys, itertools
import numpy as np
sys.path.insert(0, '..')
from _arr_longrun_fast import bit_b2_ref
from p2_b0 import enum_tapes


def gfun(a, L, R, cache=None):
    for r in range(R, L, -1):
        k = (r - L + 1) // 2
        if k < 1:
            k = 1
        # b2 == 0 exactly when sum(first k)*m >= sum(window)*k
        if bit_b2_ref(a[L:r + 1]) == 0:
            return r
    return None


def main():
    for n in (int(x) for x in sys.argv[1:]):
        tapes = list(enum_tapes(n))
        bad = 0
        pairs = 0
        examples = []
        for a in tapes:
            L, R = 0, n - 1
            gs = []
            while R - L >= 1:
                g = gfun(a, L, R)
                gs.append((L, R, g))
                if g is None:
                    break
                if g == R:
                    L, R = L + 1, R - 1
                else:
                    L += 1
            # check monotonicity along consecutive left cuts with equal R
            for i in range(len(gs) - 1):
                (L1, R1, g1), (L2, R2, g2) = gs[i], gs[i + 1]
                if R1 == R2 and g1 is not None and g2 is not None and L2 == L1 + 1:
                    pairs += 1
                    if g2 > g1:
                        bad += 1
                        if len(examples) < 3:
                            examples.append((a.tolist(), L1, R1, g1, g2))
        print('n=%d tapes=%d  consecutive-left-cut pairs=%d  violations g(L+1)>g(L): %d'
              % (n, len(tapes), pairs, bad))
        for e in examples:
            print('   ex', e)


if __name__ == '__main__':
    main()
