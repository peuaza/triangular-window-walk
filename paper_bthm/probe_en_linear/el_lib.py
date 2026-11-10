# -*- coding: utf-8 -*-
"""el_lib -- epoch trace + tape draws for the O(n)-on-the-meeting-domain probe.

WHY
---
CES_EXACT showed that for `pattern` and `runs` the walk's per-epoch outcome is
fixed by (b2,b1) at the epoch start, that the only error is one terminal epoch
per tape, and that on MEETING tapes the terminal is tiny (M0 <= 17) so a fixed
W2 buffer makes the machine exact at 2.4-2.8 n reads.  The same question is asked
here of the other fifteen rules: is there a machine that is exact on the meeting
domain at O(n) reads?

The trace is the evidence.  One walk per tape, recording per epoch

     M0   = R-L+1 at the epoch start (the region width)
     b2   = the majority bit at the epoch start
     b1   = the rule bit at the epoch start
     out  = how the epoch ended  (0 meet, 1 cutL, 2 cutR, 3 cutB, 4 narrow,
            5 cycle, 6 cap)
     ln   = steps inside the epoch

and the walk's own (delta, end, steps, cuts, inc, msum) so every cell can be
gated against pcr_lib.Runner -- the machine the record's exponents were
measured with.  The predicate evaluation is imported from pcr_lib unchanged
(b1_of / b2_of), so the walk is the record's walk, not a re-implementation.

REFERENCE AND GATE
------------------
pcr_lib.walk is the definition (PROBE: it is the excluded-middle four-move table
with conv 'a' for every rule but canonical, and the header reset after a cut).
el_gate.py checks trace_walk against pcr_lib.Runner on every rule and every tape
at n = 48, 64, 128: delta, end, steps, cuts, inc reads and msum must all agree.
Nothing else in this directory is allowed to run before that gate passes.
"""
import os
import sys

import numpy as np
from numba import njit

HERE = os.path.dirname(os.path.abspath(__file__))
PAPER = os.path.dirname(HERE)
for _d in ('probe_class_rerun', 'probe_e_uniform', 'probe_triple_nat', 'probe_hybrid3'):
    _p = os.path.join(PAPER, _d)
    if _p not in sys.path:
        sys.path.insert(0, _p)

import pcr_lib as P                                              # noqa: E402
from pcr_lib import (LPLUS, LMINUS, RPLUS, RMINUS, CUTL, CUTR, CUTB,   # noqa: E402
                     HALT, TABLE_NP, b1_of, b2_of)
import tr_bal4                                                   # noqa: E402

# BAL is one device at two thresholds: 2 is pcr_lib's BAL>m/2, 4 is the
# catalogue's BAL>m/4 (tr_bal4's own gate: thr=2 == pcr_lib).  The threshold is
# passed into the trace so both rules are one walk with one b1.
BALTHR = {'BAL>m/2': 2, 'BAL>m/4': 4}
# NBLK is one device at three saturations: the count of visits to the level
# Ps[ell] in (ell, r+1].  pcr_lib has NBLK=1 and NBLK>2 but not NBLK>3, so the
# same query is written once with a saturation parameter and gated against
# pcr_lib's two (el_gate: NBLK=1 and NBLK>2 rows).
NBLSAT = {'NBLK=1': 1, 'NBLK>2': 3, 'NBLK>3': 4}
# rule id for the trace's numba dispatch.  BAL>m/4 is the BAL device
# (R_BALH, threshold passed separately); NBLK>3 is the NBLK device
# (R_NBLK2, saturation passed separately).  Without this, BAL>m/4 would
# fall through to b1_of's NBLK branch and measure the wrong predicate.
RID_OF = {'BAL>m/4': P.R_BALH, 'NBLK>3': P.R_NBLK2}


@njit(cache=True)
def nblk_b1(Pa, Ps, psmin, bstart, bpos, sat, ell, r):
    """1 iff #{j in [ell+1, r+1] : Ps[j] == Ps[ell]} >= sat (sat=1: == 1)."""
    b = Ps[ell] - psmin
    j0 = bstart[b]
    j1 = bstart[b + 1]
    lo, hi = j0, j1
    while lo < hi:
        mid = (lo + hi) >> 1
        if bpos[mid] < ell + 1:
            lo = mid + 1
        else:
            hi = mid
    a0 = lo
    lo, hi = j0, j1
    while lo < hi:
        mid = (lo + hi) >> 1
        if bpos[mid] <= r + 1:
            lo = mid + 1
        else:
            hi = mid
    nb = lo - a0
    if sat == 1:
        return 1 if nb == 1 else 0
    return 1 if nb >= sat else 0

# rules of this probe, grouped by the device their statistic needs (see
# EN_LINEAR.md sec.1).  The grouping is a claim, gated in el_eval.py.
MONOID = ['dense2', 'peak', 'trip', 'oddrun', 'canonical', 'LAS', 'KCHG']
RANGE_MAX = ['MAXP', 'AMAX', 'AMIN']
RANGE_LEVEL = ['BAL>m/2', 'BAL>m/4', 'NBLK=1', 'NBLK>2', 'NBLK>3']
RULES = MONOID + RANGE_MAX + RANGE_LEVEL

# how b1 is evaluated at the region, under the O(log n)-bit state rule
EVAL = {}
for _r in MONOID:
    EVAL[_r] = 'sliding'          # one counter, 1 read per cell given up
for _r in RANGE_MAX:
    EVAL[_r] = 'rangemax'         # range max of the prefix walk
for _r in RANGE_LEVEL:
    EVAL[_r] = 'level'            # a level index: count / last occurrence

NS = (512, 1024, 2048, 4096)


def tape(n, i):
    """T_n draw i -- the record's family, byte for byte (eu_lib.make_tape)."""
    rng = np.random.default_rng(5000 + n + 100003 * i)
    h = n // 2
    X = np.zeros((n, 2), np.int8)
    for c in range(2):
        X[rng.choice(n, h, replace=False), c] = 1
    return X


@njit(cache=True)
def trace_walk(rule, conv, a, s, Pa, Paa, PR, Pd, Pn2, Ppk, Ptr, Ps, cpos, pnz, C2,
               bstart, bpos, psmin, SPs, SNs, LOG2, r2s, r2e, SP2, rcs, rce, rcodd, SPC,
               R, IP, HX, Ppat, d1, d2, stamp, genbox, n, cap, balthr, nbsat,
               epM, epb2, epb1, epout, eplen):
    """pcr_lib.walk, line for line, plus one record per epoch.

    Returns (out, steps, msum, inc, cuts, end, nep).
    """
    L, R_ = 0, n - 1
    ell, r = L, R_
    y, x = L, R_
    hy, hx = 0, max(1, n // 2) - 1
    gen = genbox[0] + 1
    steps = 0

    msum = 0
    cuts = 0
    inc = n + 2 * (n - 1) + 2 + max(1, n // 2)
    out = -1.0
    end = 1
    nep = 0
    curM = 0
    cur_b1 = -1
    cur_b2 = -1
    curl = 0
    fresh = 1

    while True:
        if ell >= r:
            out = ell / (n - 1.0)
            end = 0
            epM[nep] = curM if not fresh else (R_ - L + 1)
            epb2[nep] = cur_b2
            epb1[nep] = cur_b1
            epout[nep] = 0
            eplen[nep] = 0 if fresh else curl
            nep += 1
            break
        if R_ - L + 1 < 3:
            end = 1
            epM[nep] = curM if not fresh else (R_ - L + 1)
            epb2[nep] = cur_b2
            epb1[nep] = cur_b1
            epout[nep] = 4
            eplen[nep] = 0 if fresh else curl
            nep += 1
            break
        steps += 1
        if steps > cap:
            end = 3
            epM[nep] = curM if not fresh else (R_ - L + 1)
            epb2[nep] = cur_b2
            epb1[nep] = cur_b1
            epout[nep] = 6
            eplen[nep] = 0 if fresh else curl
            nep += 1
            break
        key = ell * (n + 1) + r
        if stamp[key] == gen:
            end = 2
            epM[nep] = curM if not fresh else (R_ - L + 1)
            epb2[nep] = cur_b2
            epb1[nep] = cur_b1
            epout[nep] = 5
            eplen[nep] = 0 if fresh else curl
            nep += 1
            break
        stamp[key] = gen
        m = r - ell + 1
        msum += m
        b1 = b1_of(rule, a, Pa, Paa, PR, Pd, Pn2, Ppk, Ptr, Ps, cpos, pnz, C2,
                   bstart, bpos, psmin, SPs, SNs, LOG2, r2s, r2e, SP2, rcs, rce, rcodd, SPC,
                   R, IP, HX, Ppat, d1, d2, ell, r)
        if rule == P.R_BALH:
            b1 = tr_bal4.bal_b1(a, Pa, Ps, psmin, bstart, bpos, balthr, ell, r)
        elif rule == P.R_NBLK2 or rule == P.R_NBLK1:
            b1 = nblk_b1(Pa, Ps, psmin, bstart, bpos, nbsat, ell, r)
        if rule == P.R_HASHMOD:
            b2 = 0
        elif rule == P.R_PC or rule == P.R_CX:
            b2 = 1 if (R_ - L + 1) - m < max(1, (R_ - L + 1) // 64) else 0
        else:
            b2 = b2_of(a, Pa, ell, r)
        if fresh == 1:
            curM = R_ - L + 1
            cur_b1 = b1
            cur_b2 = b2
            curl = 0
            fresh = 0
        curl += 1
        if rule == P.R_HASHMOD:
            act = CUTL if b1 == 1 else RMINUS
        elif rule == P.R_PC or rule == P.R_CX:
            W = R_ - L + 1
            if m <= 2:
                act = LPLUS
            elif W - m < max(1, W // 64):
                act = LPLUS if b1 == 1 else RMINUS
            else:
                act = CUTL if ((ell - L) & 1) else CUTR
        else:
            act = TABLE_NP[2 * b2 + b1]
        if act >= CUTL:
            cuts += 1
            if act == CUTL:
                L += 1
                oc = 1
            elif act == CUTR:
                R_ -= 1
                oc = 2
            else:
                L += 1
                R_ -= 1
                oc = 3
            if L > R_:
                end = 1
                epM[nep] = curM
                epb2[nep] = cur_b2
                epb1[nep] = cur_b1
                epout[nep] = oc
                eplen[nep] = curl
                nep += 1
                break
            ell, r = L, R_
            gen += 1
            inc += abs(y - L) + abs(x - R_)
            y, x = L, R_
            k = (x - y + 1) >> 1
            if k == 0:
                k = 1
            ty, tx = y, y + k - 1
            inc += abs(hy - ty) + abs(hx - tx)
            hy, hx = ty, tx
            epM[nep] = curM
            epb2[nep] = cur_b2
            epb1[nep] = cur_b1
            epout[nep] = oc
            eplen[nep] = curl
            nep += 1
            fresh = 1
            continue
        dy = 1 if act == LPLUS else (-1 if act == LMINUS else 0)
        dx = 1 if act == RPLUS else (-1 if act == RMINUS else 0)
        if ell + dy < L or r + dx > R_:
            cuts += 1
            if conv == 1:
                if ell + dy < L:
                    R_ -= 1
                    oc = 2
                else:
                    L += 1
                    oc = 1
            else:
                if ell == L and r == R_:
                    L += 1
                    R_ -= 1
                    y += 1
                    x -= 1
                    inc += 2
                    oc = 3
                elif ell + dy < L:
                    L += 1
                    oc = 1
                else:
                    R_ -= 1
                    oc = 2
            if L > R_:
                end = 1
                epM[nep] = curM
                epb2[nep] = cur_b2
                epb1[nep] = cur_b1
                epout[nep] = oc
                eplen[nep] = curl
                nep += 1
                break
            ell, r = L, R_
            gen += 1
            inc += abs(y - L) + abs(x - R_)
            y, x = L, R_
            k = (x - y + 1) >> 1
            if k == 0:
                k = 1
            ty, tx = y, y + k - 1
            inc += abs(hy - ty) + abs(hx - tx)
            hy, hx = ty, tx
            epM[nep] = curM
            epb2[nep] = cur_b2
            epb1[nep] = cur_b1
            epout[nep] = oc
            eplen[nep] = curl
            nep += 1
            fresh = 1
            continue
        if ell + dy >= r + dx:
            out = (ell + dy) / (n - 1.0)
            end = 0
            epM[nep] = curM
            epb2[nep] = cur_b2
            epb1[nep] = cur_b1
            epout[nep] = 0
            eplen[nep] = curl
            nep += 1
            break
        if r + dx < L or ell + dy > R_:
            L += 1
            R_ -= 1
            ell, r = L, R_
            gen += 1
            inc += abs(y - L) + abs(x - R_)
            y, x = L, R_
            k = (x - y + 1) >> 1
            if k == 0:
                k = 1
            ty, tx = y, y + k - 1
            inc += abs(hy - ty) + abs(hx - tx)
            hy, hx = ty, tx
            epM[nep] = curM
            epb2[nep] = cur_b2
            epb1[nep] = cur_b1
            epout[nep] = 3
            eplen[nep] = curl
            nep += 1
            fresh = 1
            continue
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
    genbox[0] = gen
    return out, steps, msum, inc, cuts, end, nep


class Mach(object):
    """One tape's structures plus a stamp array that is reused across tapes."""

    def __init__(self, X, n):
        self.S = P.St(X)
        self.n = n
        if not hasattr(Mach, '_stamp') or Mach._stamp_n != n:
            Mach._stamp = np.zeros((n + 1) * (n + 1), np.int32)
            Mach._stamp_n = n
            Mach._genbox = np.zeros(1, np.int64)
        self.stamp = Mach._stamp
        self.genbox = Mach._genbox

    def trace(self, rule, cap=None, balthr=None, nbsat=None):
        S = self.S
        n = self.n
        if cap is None:
            cap = 40 * n * n + 100
        if balthr is None:
            balthr = BALTHR.get(rule, 2)
        if nbsat is None:
            nbsat = NBLSAT.get(rule, 3)
        epM = np.zeros(n + 4, np.int32)
        epb2 = np.zeros(n + 4, np.int8)
        epb1 = np.zeros(n + 4, np.int8)
        epout = np.zeros(n + 4, np.int8)
        eplen = np.zeros(n + 4, np.int64)
        cv = P.CONV.get(rule, 0)
        rid_ = RID_OF[rule] if rule in RID_OF else P.RID[rule]
        r = trace_walk(rid_, cv, S.a, S.s, S.Pa, S.Paa, S.PR, S.Pd,
                       S.Pn2, S.Ppk, S.Ptr, S.Ps, S.cpos, S.pnz, S.C2,
                       S.bstart, S.bpos, S.psmin, S.SPs, S.SNs, S.LOG2,
                       S.r2s, S.r2e, S.SP2, S.rcs, S.rce, S.rcodd, S.SPC,
                       S.R, S.IP, S.HX, S.Ppat, S.d1, S.d2, self.stamp,
                       self.genbox, n, cap, int(balthr), int(nbsat),
                       epM, epb2, epb1, epout, eplen)
        out, steps, msum, inc, cuts, end, nep = r
        return dict(delta=out, end=int(end), steps=int(steps), msum=int(msum),
                    inc=int(inc), cuts=int(cuts), nep=int(nep),
                    epM=epM[:nep].copy(), epb2=epb2[:nep].copy(),
                    epb1=epb1[:nep].copy(), epout=epout[:nep].copy(),
                    eplen=eplen[:nep].copy())
