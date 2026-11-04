# -*- coding: utf-8 -*-
"""p5: (a) violations of g(L+1)<=g(L) only when g(L+1)=R?  (b) double-cut census."""
import sys, itertools
import numpy as np
sys.path.insert(0, '..')
from _arr_longrun_fast import bit_b2_ref
from p2_b0 import enum_tapes
from p4_monotone_g import gfun


def main():
    for n in (int(x) for x in sys.argv[1:]):
        tapes = list(enum_tapes(n))
        viol = violR = 0
        ndouble = nleft = 0
        dbl_per_tape = []
        for a in tapes:
            L, R = 0, n - 1
            gs = []
            d = l = 0
            while R - L >= 1:
                g = gfun(a, L, R)
                gs.append((L, R, g))
                if g is None:
                    break
                if g == R:
                    d += 1
                    L, R = L + 1, R - 1
                else:
                    l += 1
                    L += 1
            ndouble += d
            nleft += l
            dbl_per_tape.append(d)
            for i in range(len(gs) - 1):
                (L1, R1, g1), (L2, R2, g2) = gs[i], gs[i + 1]
                if R1 == R2 and g1 is not None and g2 is not None and L2 == L1 + 1 and g2 > g1:
                    viol += 1
                    if g2 == R2:
                        violR += 1
        print('n=%d tapes=%d  violations=%d of which g2==R2: %d   double cuts/tape=%.3f left cuts/tape=%.3f'
              % (n, len(tapes), viol, violR, ndouble / len(tapes), nleft / len(tapes)))
        print('    double cuts per tape: mean %.2f median %d max %d ; frac of tapes with 0 doubles %.3f'
              % (np.mean(dbl_per_tape), np.median(dbl_per_tape), max(dbl_per_tape),
                 np.mean(np.array(dbl_per_tape) == 0)))


if __name__ == '__main__':
    main()
