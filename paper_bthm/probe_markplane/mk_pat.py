# -*- coding: utf-8 -*-
"""mk_pat -- pcr_lib.walk for rule `pattern` with a one-bit mark plane.

The machine is pcr_lib.Runner/Walk under the INC charge ("1 read per cell taken
in or given up + re-span per cut"): the opening pass costs
n + 2(n-1) + 2 + max(1,n//2), every legal move 1 read, the head re-sync its
distance, and every cut its re-span plus head re-sync.  The walk, the table, the
cut convention, delta and end are pcr_lib's; mk_gate checks that policy 0
reproduces Runner(...).run('pattern') on every tape (reads, delta, end, steps,
cuts).

policies
  0  base      the machine, no plane.
  1  adj       m[i] = [a[i+1] == a[i]], written while the opening pass reads the
               tape cell by cell.  The adjacent-pair pass (the 2(n-1) term) then
               needs one read per pair instead of two: reading cell i returns
               a[i] and m[i], and m[i]=1 gives a[i+1] with no read.  The only
               place in either machine where one read can yield two cells.
  2  occ       m[i] = [a[i-2..i] == 2,1,1] ("an occurrence ends here"), written in
               the opening pass.  phi is then a range OR over marks, evaluated
               by reading them: the composable range summary a 1-bit plane can
               hold.  Consultations are charged at 1 read per mark.
  3  occ_null  the same marks, written but never consulted (phi stays on the
               machine's own prefix array, which the INC charge gives free).
               The best case of policy 2.
"""
import os
import sys

import numpy as np
from numba import njit

HERE = os.path.dirname(os.path.abspath(__file__))
PAPER = os.path.dirname(HERE)
for d in ('probe_class_rerun', 'probe_triple_nat'):
    p = os.path.join(PAPER, d)
    if p not in sys.path:
        sys.path.insert(0, p)

import pcr_lib as P                                             # noqa: E402

TABLE = np.array(P.TABLE, np.int8)
CUTL, CUTR = P.CUTL, P.CUTR


@njit(cache=True)
def walk_mk(policy, a, Pa, Ppat, n, cap, mk, stamp, genbox):
    """One tape.  Returns (out, steps, msum, inc, cuts, end, ink, consults)."""
    L, R_ = 0, n - 1
    ell, r = L, R_
    y, x = L, R_
    hy, hx = 0, max(1, n // 2) - 1
    gen = genbox[0] + 1
    steps = 0
    msum = 0
    cuts = 0
    ink = 0
    consults = 0
    kb = max(1, n // 2)

    if policy == 1:
        # opening: the S pass reads every cell and writes the adjacency bit,
        # the adjacent-pair pass then reads one cell per pair
        prev = a[0]
        for i in range(1, n):
            mk[i - 1] = 1 if a[i] == prev else 0
            ink += 1
            prev = a[i]
        rd = 0
        for i in range(n - 1):
            if mk[i] == 1:
                rd += 1
            else:
                rd += 2
        inc = n + rd + 2 + kb
    else:
        inc = n + 2 * (n - 1) + 2 + kb
        if policy == 2 or policy == 3:
            for i in range(2, n):
                if a[i - 2] == 2 and a[i - 1] == 1 and a[i] == 1:
                    mk[i] = 1
                    ink += 1

    out = -1.0
    end = 1
    while True:
        if ell >= r:
            out = ell / (n - 1.0)
            end = 0
            break
        if R_ - L + 1 < 3:
            end = 1
            break
        steps += 1
        if steps > cap:
            end = 3
            break
        key = ell * (n + 1) + r
        if stamp[key] == gen:
            end = 2
            break
        stamp[key] = gen
        m = r - ell + 1
        msum += m
        if policy == 2:
            b1 = 0
            if m >= 3:
                j = ell + 2
                while j <= r:
                    consults += 1
                    if mk[j] == 1:
                        b1 = 1
                        break
                    j += 1
        else:
            b1 = 0
            if m >= 3 and (Ppat[r - 1] - Ppat[ell]) > 0:
                b1 = 1
        k = m >> 1
        if k == 0:
            k = 1
        b2 = 1 if (Pa[r + 1] - Pa[ell]) * k > (Pa[ell + k] - Pa[ell]) * m else 0
        act = TABLE[2 * b2 + b1]
        if act >= CUTL:
            cuts += 1
            if act == CUTL:
                L += 1
            elif act == CUTR:
                R_ -= 1
            else:
                L += 1
                R_ -= 1
            if L > R_:
                end = 1
                break
            ell, r = L, R_
            inc += abs(y - L) + abs(x - R_)
            y, x = L, R_
            k = (x - y + 1) >> 1
            if k == 0:
                k = 1
            ty, tx = y, y + k - 1
            inc += abs(hy - ty) + abs(hx - tx)
            hy, hx = ty, tx
            continue
        dy = 1 if act == P.LPLUS else (-1 if act == P.LMINUS else 0)
        dx = 1 if act == P.RPLUS else (-1 if act == P.RMINUS else 0)
        if ell + dy < L or r + dx > R_:
            cuts += 1
            if ell == L and r == R_:
                L += 1
                R_ -= 1
                y += 1
                x -= 1
                inc += 2
            elif ell + dy < L:
                L += 1
            else:
                R_ -= 1
            if L > R_:
                end = 1
                break
            ell, r = L, R_
            inc += abs(y - L) + abs(x - R_)
            y, x = L, R_
            k = (x - y + 1) >> 1
            if k == 0:
                k = 1
            ty, tx = y, y + k - 1
            inc += abs(hy - ty) + abs(hx - tx)
            hy, hx = ty, tx
            continue
        if ell + dy >= r + dx:
            out = (ell + dy) / (n - 1.0)
            end = 0
            break
        ell += dy
        r += dx
        y, x = ell, r
        inc += 1
        k = (x - y + 1) >> 1
        if k == 0:
            k = 1
        ty, tx = y, y + k - 1
        inc += abs(hy - ty) + abs(hx - tx)
        hy, hx = ty, tx
    if policy == 2:
        inc += consults          # every mark consultation costs one read
    genbox[0] = gen
    return out, steps, msum, inc, cuts, end, ink, consults


class PatMarked(object):
    """Holds one tape's structures, marks and stamp, as pcr_lib.Runner does."""

    def __init__(self, X):
        S = self.S = P.St(X)
        self.n = S.n
        self.mk = np.zeros(S.n, np.int8)
        self.stamp = np.zeros((S.n + 1) * (S.n + 1), np.int32)
        self.genbox = np.zeros(1, np.int64)

    def run(self, policy, cap=None):
        S = self.S
        n = self.n
        if cap is None:
            cap = 40 * n * n + 100
        out, steps, msum, inc, cuts, end, ink, cons = walk_mk(
            policy, S.a, S.Pa, S.Ppat, n, cap, self.mk, self.stamp, self.genbox)
        return dict(reads=int(inc), ink=int(ink), delta=float(out), end=int(end),
                    steps=int(steps), cuts=int(cuts), scan=int(msum),
                    consults=int(cons))


POLICIES = ('none', 'adj', 'occ', 'occ_null')
PCODE = {'none': 0, 'adj': 1, 'occ': 2, 'occ_null': 3}


def run_pat(X, policy='none'):
    return PatMarked(X).run(PCODE[policy])


def run_base(X):
    out, steps, msum, inc, cuts, end, _ = P.Runner(X).run('pattern')
    return dict(reads=int(inc), ink=0, delta=float(out), end=int(end),
                steps=int(steps), cuts=int(cuts), scan=int(msum), consults=0)
