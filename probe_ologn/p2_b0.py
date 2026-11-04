# -*- coding: utf-8 -*-
"""p2: structure of the b1:=0 variant.  Does it admit a closed form / jump?

b1:=0 move table:  b2=1 -> r-1 ;  b2=0 -> ell-1.
Starting an epoch at (L,R): b2([L,R])=0 gives a double cut; otherwise r decreases
until the first r with b2([L,r])=0, which gives a left cut; if no such r before
r=L, the run meets and outputs L/(n-1).
"""
import sys, itertools, json
import numpy as np
sys.path.insert(0, '..')
from _arr_longrun_fast import walk_naive, bit_b2_ref


def walk_b0(a, n):
    """b1 forced to 0; same control flow as walk_naive."""
    L, R = 0, n - 1
    ell, r = 0, n - 1
    steps = 0
    seen = set()
    while True:
        if ell >= r:
            return ell / (n - 1), steps
        if R - L + 1 < 3:
            return None, steps
        steps += 1
        if steps > 20 * n * n:
            return None, steps
        if (ell, r) in seen:
            return None, steps
        seen.add((ell, r))
        b2 = int(bit_b2_ref(a[ell:r + 1]))
        dy, dx = (0, -1) if b2 else (-1, 0)
        if ell + dy < L or r + dx > R:
            if ell == L and r == R:
                L, R = L + 1, R - 1
            elif ell + dy < L:
                L += 1
            else:
                R -= 1
            ell, r = L, R
            seen.clear()
            continue
        if ell + dy >= r + dx:
            return (ell + dy) / (n - 1), steps
        if r + dx < L or ell + dy > R:
            L, R = L + 1, R - 1
            ell, r = L, R
            seen.clear()
            continue
        ell, r = ell + dy, r + dx


from p0_tapes import enum_tapes


def scan_for_L(a, L, R):
    """the epoch starting at (L,R): return ('meet',L) or ('cut', side, r_first0)."""
    for r in range(R, L, -1):
        if bit_b2_ref(a[L:r + 1]) == 0:
            return ('cut', 'double' if r == R else 'left', r)
    return ('meet', L, None)


def main():
    for n in (int(x) for x in sys.argv[1:]):
        tapes = list(enum_tapes(n))
        agree = 0
        meetL_hist = {}
        for a in tapes:
            d0, s0 = walk_b0(a, n)
            # replay the (L,R) sequence and compare each epoch against scan_for_L
            L, R = 0, n - 1
            ok = True
            out = None
            while True:
                if L >= R:
                    out = L / (n - 1)
                    break
                if R - L + 1 < 3:
                    out = None
                    break
                kind = scan_for_L(a, L, R)
                if kind[0] == 'meet':
                    out = L / (n - 1)
                    break
                _, side, rf = kind
                if side == 'double':
                    L, R = L + 1, R - 1
                else:
                    L += 1
            if (out is None and d0 is None) or (out is not None and d0 is not None and abs(out - d0) < 1e-12):
                agree += 1
            if d0 is not None:
                meetL_hist[int(round(d0 * (n - 1)))] = meetL_hist.get(int(round(d0 * (n - 1))), 0) + 1
        print('n=%d tapes=%d  closed-form-vs-b0 agreement %d/%d  meet-at-L hist %s'
              % (n, len(tapes), agree, len(tapes), dict(sorted(meetL_hist.items()))))


if __name__ == '__main__':
    main()
