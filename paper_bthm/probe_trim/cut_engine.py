# -*- coding: utf-8 -*-
"""cyc_engine -- the model's walk with the state quadruple exposed.

WHY
---
tr_lib's machines carry a stamp on the WINDOW PAIR (ell, r) and stop the walk on
a repeat (end=2); nothing in the record ever sees the region, so nothing in the
record can tell a true state cycle from a projected one.  The state of the walk is
the quadruple (ell, r, L, R): the tape is fixed and b1, b2 are functions of
(ell, r) alone, so the next quadruple is a function of the current one.  This
module runs exactly that walk and, instead of stopping, asks the one-bit oracle
"has this quadruple occurred before in this run?".

The predicates are NOT re-implemented: b1 is pcr_lib.b1_of (the eighteen-rule
predicate module the whole record uses), b2 is pcr_lib.b2_of; BAL>m/4 is
tr_bal4.bal_b1 with thr=4; NBLK>3 is the same visit count as pcr_lib's NBLK>2
with the threshold at 3.  The read accounting is pcr_lib's INC accounting,
copied literally (opening pass n + 2(n-1) + 2 + n//2, one read per cell taken in
or given up, re-span after a cut, head re-read), so the numbers are comparable
across rules and across the two models.

MODES
-----
  0  BTHM-0    no memory at all; the walk stops at the first repeated quadruple
               (by determinism that is the point of no return: from there it
               repeats forever and dies at the step cap).
  1  BTHM-C    the one-bit oracle; on a repeat the rule trims one cell off the
               window (trim=1 left, trim=2 right), pays 1 read + the head
               re-read, and continues.
  2  LEGACY    pcr_lib's own machine: stamp on (ell, r), stop on a repeat
               (end=2).  Used only to check this engine against tr_lib.
  3  VERIFY    mode 0, then run one more period and report the second repeat,
               so the periodicity claim can be checked rather than asserted.
  4  BTHM-CUT  BTHM-C with the repair as a real trim: on a repeat the rule
               removes the BOUNDARY CELL OF THE TAPE on the side its own action
               table selects for that state (L <- L+1 or R <- R-1), pays the
               re-span exactly as a rule-initiated cut does (abs(y-L) +
               abs(x-R_) plus the head re-read), re-opens the window on the
               narrower tape and continues.  The tape only shrinks, so at most
               n repairs are possible; this is mode 4, the mode probe_cycle_cut
               measures.

END CODES: 0 meet, 1 region empty, 2 legacy (ell,r) repeat, 3 step cap,
           4 quadruple cycle (BTHM-0 death), 5 table full (unresolved),
           6 step cap reached after at least one repair.
"""
import os
import sys

import numpy as np
from numba import njit

HERE = os.path.dirname(os.path.abspath(__file__))
PAPER = os.path.dirname(HERE)
for _d in ('probe_class_rerun', 'probe_e_uniform', 'probe_triple_nat'):
    _p = os.path.join(PAPER, _d)
    if _p not in sys.path:
        sys.path.insert(0, _p)

import pcr_lib as P                       # noqa: E402
import tr_bal4                            # noqa: E402
from pcr_lib import b1_of, b2_of          # noqa: E402

LPLUS, LMINUS, RPLUS, RMINUS, CUTL, CUTR, CUTB, HALT = range(8)
TAB = P.TABLE_NP


@njit(cache=True)
def _b1_variant(variant, rid, a, Pa, Paa, PR, Pd, Pn2, Ppk, Ptr, Ps, cpos, pnz,
                C2, bstart, bpos, psmin, SPs, SNs, LOG2, r2s, r2e, SP2, rcs,
                rce, rcodd, SPC, R, IP, HX, Ppat, d1, d2, ell, r):
    if variant == 1:                                  # BAL>m/4
        return tr_bal4.bal_b1(a, Pa, Ps, psmin, bstart, bpos, 4, ell, r)
    if variant == 2:                                  # NBLK>3 (nb > 3)
        b = Ps[ell] - psmin
        j0 = bstart[b]
        j1 = bstart[b + 1]
        lo = j0
        hi = j1
        while lo < hi:                                # first bpos > ell
            mid = (lo + hi) >> 1
            if bpos[mid] <= ell:
                lo = mid + 1
            else:
                hi = mid
        k0 = lo
        lo = j0
        hi = j1
        while lo < hi:                                # first bpos > r+1
            mid = (lo + hi) >> 1
            if bpos[mid] <= r + 1:
                lo = mid + 1
            else:
                hi = mid
        k1 = lo
        return 1 if (k1 - k0) > 3 else 0
    return b1_of(rid, a, Pa, Paa, PR, Pd, Pn2, Ppk, Ptr, Ps, cpos, pnz, C2,
                 bstart, bpos, psmin, SPs, SNs, LOG2, r2s, r2e, SP2, rcs, rce,
                 rcodd, SPC, R, IP, HX, Ppat, d1, d2, ell, r)


@njit(cache=True)
def walk(variant, rid, conv, a, Pa, Paa, PR, Pd, Pn2, Ppk, Ptr, Ps, cpos, pnz,
         C2, bstart, bpos, psmin, SPs, SNs, LOG2, r2s, r2e, SP2, rcs, rce,
         rcodd, SPC, R, IP, HX, Ppat, d1, d2, n, cap, mode, trim,
         hkeys, hvals, hmask, stamp, genbox, org):
    """One tape.  org[0..11] = steps, reads, cuts, end, cyc_start, cyc_len,
    width, ntrim, second_start, second_step, out, msum; org[12] = repairs that
    removed the left boundary cell, org[13] = the right one, org[14] = hash
    slots used."""
    L = 0
    R_ = n - 1
    ell = L
    r = R_
    y = L
    x = R_
    hy = 0
    hx = max(1, n // 2) - 1
    gen = genbox[0] + 1
    steps = 0
    cuts = 0
    ntrim = 0
    msum = 0
    reads = n + 2 * (n - 1) + 2 + max(1, n // 2)
    out = -1.0
    end = 1
    cstart = -1
    clen = -1
    cwidth = -1
    j0 = -1
    target = -1
    s2start = -1
    s2step = -1
    nused = 0
    nl = 0
    nr = 0
    found = False
    while True:
        if ell >= r:
            out = ell / (n - 1.0)
            end = 0
            break
        if R_ - L + 1 < 3:
            end = 1
            break
        steps += 1
        if steps + ntrim > cap:
            end = 3 if ntrim == 0 else 6
            break
        if mode == 2:
            key2 = ell * (n + 1) + r
            if stamp[key2] == gen:
                end = 2
                break
            stamp[key2] = gen
        m = r - ell + 1
        msum += m
        # ---- the oracle: one bit, "has this quadruple occurred before?"
        if mode == 0 or mode == 1 or mode == 3 or mode == 4:
            key = ((ell * (n + 1) + r) * (n + 1) + L) * (n + 1) + R_
            h = key * 0x9E3779B97F4A7C15
            h = (h ^ (h >> 29)) * 0xBF58476D1CE4E5B9
            idx = (h ^ (h >> 32)) & hmask
            first = -1
            while hkeys[idx] != -1:
                if hkeys[idx] == key:
                    first = hvals[idx]
                    break
                idx = (idx + 1) & hmask
            if first >= 0:
                if found:
                    if steps == target:
                        s2start = first
                        s2step = steps
                        end = 4
                        break
                    if steps > target:
                        s2start = -2          # the period did not close
                        s2step = steps
                        break
                else:
                    cstart = first
                    clen = steps - first
                    cwidth = m
                    if mode == 0:
                        end = 4
                        break
                    if mode == 3:
                        found = True          # keep walking: one full period
                        j0 = steps
                        target = steps + clen
                        # falls through to the test, no trim
                    elif mode == 4:
                        # ---- THE REPAIR, AS A REAL TRIM ----------------
                        # The triple repeats; the model's repair removes the
                        # boundary cell of the TAPE on the side its own action
                        # table selects for this state, and continues.  The
                        # side is read off the same table entry the step block
                        # below would have used, with the c1a rules inverted,
                        # exactly as the trim=5 variant does.  Charge: the
                        # re-span a rule-initiated cut pays (the cells taken in
                        # / given up plus the head re-read), no oracle read.
                        ntrim += 1
                        v1 = _b1_variant(variant, rid, a, Pa, Paa, PR, Pd,
                                         Pn2, Ppk, Ptr, Ps, cpos, pnz, C2,
                                         bstart, bpos, psmin, SPs, SNs,
                                         LOG2, r2s, r2e, SP2, rcs, rce,
                                         rcodd, SPC, R, IP, HX, Ppat, d1,
                                         d2, ell, r)
                        if variant == 3:
                            act = CUTL if v1 == 1 else RMINUS
                        else:
                            v2 = b2_of(a, Pa, ell, r)
                            act = TAB[2 * v2 + v1]
                        side = 1 if (act == LPLUS or act == LMINUS
                                     or act == CUTL) else 2
                        if conv == 1:
                            side = 3 - side
                        if side == 1:
                            L += 1
                            nl += 1
                        else:
                            R_ -= 1
                            nr += 1
                        if L > R_:
                            end = 1
                            break
                        ell = L
                        r = R_
                        gen += 1
                        reads += abs(y - L) + abs(x - R_)
                        y = L
                        x = R_
                        k = (x - y + 1) >> 1
                        if k == 0:
                            k = 1
                        ty = y
                        tx = y + k - 1
                        reads += abs(hy - ty) + abs(hx - tx)
                        hy = ty
                        hx = tx
                        continue
                    else:
                        # mode 1: one trim per repeat, then continue.  The side is
                        # the rule's own cut side for that state: the direction of
                        # the action the rule takes here (a cut action names its
                        # side directly), with the c1a rules inverting it -- the
                        # same parity pcr_lib uses when it turns a blocked move
                        # into a cut.
                        ntrim += 1
                        if trim == 5:
                            v1 = _b1_variant(variant, rid, a, Pa, Paa, PR, Pd,
                                             Pn2, Ppk, Ptr, Ps, cpos, pnz, C2,
                                             bstart, bpos, psmin, SPs, SNs,
                                             LOG2, r2s, r2e, SP2, rcs, rce,
                                             rcodd, SPC, R, IP, HX, Ppat, d1,
                                             d2, ell, r)
                            if variant == 3:
                                act = CUTL if v1 == 1 else RMINUS
                            else:
                                v2 = b2_of(a, Pa, ell, r)
                                act = TAB[2 * v2 + v1]
                            side = 1 if (act == LPLUS or act == LMINUS
                                         or act == CUTL) else 2
                            if conv == 1:
                                side = 3 - side
                        elif trim >= 3:
                            # a cut instead of a window trim (kept for reference)
                            if trim == 3:
                                L += 1
                            else:
                                R_ -= 1
                            if L > R_:
                                end = 1
                                break
                            ell = L
                            r = R_
                            gen += 1
                            reads += abs(y - L) + abs(x - R_)
                            y = L
                            x = R_
                            k = (x - y + 1) >> 1
                            if k == 0:
                                k = 1
                            ty = y
                            tx = y + k - 1
                            reads += abs(hy - ty) + abs(hx - tx)
                            hy = ty
                            hx = tx
                            continue
                        else:
                            side = trim
                        if side == 1:
                            ell += 1
                        else:
                            r -= 1
                        y = ell
                        x = r
                        reads += 1
                        k = (x - y + 1) >> 1
                        if k == 0:
                            k = 1
                        ty = y
                        tx = y + k - 1
                        reads += abs(hy - ty) + abs(hx - tx)
                        hy = ty
                        hx = tx
                        continue
            else:
                hkeys[idx] = key
                hvals[idx] = steps
                nused += 1
                if nused > (hmask + 1) * 3 // 4:
                    end = 5
                    break
        b1 = _b1_variant(variant, rid, a, Pa, Paa, PR, Pd, Pn2, Ppk, Ptr, Ps,
                         cpos, pnz, C2, bstart, bpos, psmin, SPs, SNs, LOG2,
                         r2s, r2e, SP2, rcs, rce, rcodd, SPC, R, IP, HX, Ppat,
                         d1, d2, ell, r)
        if variant == 3:
            act = CUTL if b1 == 1 else RMINUS     # hashmod's own action rule
        else:
            b2 = b2_of(a, Pa, ell, r)
            act = TAB[2 * b2 + b1]
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
            ell = L
            r = R_
            gen += 1
            reads += abs(y - L) + abs(x - R_)
            y = L
            x = R_
            k = (x - y + 1) >> 1
            if k == 0:
                k = 1
            ty = y
            tx = y + k - 1
            reads += abs(hy - ty) + abs(hx - tx)
            hy = ty
            hx = tx
            continue
        dy = 1 if act == LPLUS else (-1 if act == LMINUS else 0)
        dx = 1 if act == RPLUS else (-1 if act == RMINUS else 0)
        if ell + dy < L or r + dx > R_:
            cuts += 1
            if conv == 1:
                if ell + dy < L:
                    R_ -= 1
                else:
                    L += 1
            else:
                if ell == L and r == R_:
                    L += 1
                    R_ -= 1
                    y += 1
                    x -= 1
                    reads += 2
                elif ell + dy < L:
                    L += 1
                else:
                    R_ -= 1
            if L > R_:
                end = 1
                break
            ell = L
            r = R_
            gen += 1
            reads += abs(y - L) + abs(x - R_)
            y = L
            x = R_
            k = (x - y + 1) >> 1
            if k == 0:
                k = 1
            ty = y
            tx = y + k - 1
            reads += abs(hy - ty) + abs(hx - tx)
            hy = ty
            hx = tx
            continue
        if ell + dy >= r + dx:
            out = (ell + dy) / (n - 1.0)
            end = 0
            break
        ell += dy
        r += dx
        y = ell
        x = r
        reads += 1
        k = (x - y + 1) >> 1
        if k == 0:
            k = 1
        ty = y
        tx = y + k - 1
        reads += abs(hy - ty) + abs(hx - tx)
        hy = ty
        hx = tx
    genbox[0] = gen
    org[0] = steps
    org[1] = reads
    org[2] = cuts
    org[3] = end
    org[4] = cstart
    org[5] = clen
    org[6] = cwidth
    org[7] = ntrim
    org[8] = s2start
    org[9] = s2step
    org[10] = out
    org[11] = msum
    org[12] = nl
    org[13] = nr
    org[14] = nused
    return


# ------------------------------------------------------------------- runner
HLOG = 23


class Cyc(object):
    """One tape's structures plus the oracle table, reused across tapes."""

    def __init__(self, X, hlog=HLOG):
        self.S = P.St(X)
        self.n = self.S.n
        self.hkeys = np.full(1 << hlog, -1, np.int64)
        self.hvals = np.zeros(1 << hlog, np.int32)
        self.hmask = np.int64((1 << hlog) - 1)
        self.stamp = np.zeros((self.n + 1) * (self.n + 1), np.int32)
        self.genbox = np.zeros(1, np.int64)
        self.org = np.zeros(16, np.float64)

    def run(self, rid, conv=0, variant=0, mode=0, trim=1, cap=None):
        S = self.S
        n = self.n
        if cap is None:
            cap = 40 * n * n + 100
        self.hkeys.fill(-1)
        self.hvals.fill(0)
        walk(variant, rid, conv, S.a, S.Pa, S.Paa, S.PR, S.Pd, S.Pn2, S.Ppk,
             S.Ptr, S.Ps, S.cpos, S.pnz, S.C2, S.bstart, S.bpos, S.psmin,
             S.SPs, S.SNs, S.LOG2, S.r2s, S.r2e, S.SP2, S.rcs, S.rce, S.rcodd,
             S.SPC, S.R, S.IP, S.HX, S.Ppat, S.d1, S.d2, n, cap, mode, trim,
             self.hkeys, self.hvals, self.hmask, self.stamp, self.genbox,
             self.org)
        o = self.org
        return dict(steps=int(o[0]), reads=int(o[1]), cuts=int(o[2]),
                    end=int(o[3]), cyc_start=int(o[4]), cyc_len=int(o[5]),
                    width=int(o[6]), ntrim=int(o[7]), s2start=int(o[8]),
                    s2step=int(o[9]), delta=float(o[10]) if o[10] != -1.0 else None,
                    msum=int(o[11]), nl=int(o[12]), nr=int(o[13]),
                    nused=int(o[14]))


def mono_run(a):
    """lin_engine.mono_walk: ell += 1 every step, no cut, reads 2n-1."""
    n = len(a)
    return dict(steps=n - 1, reads=2 * n - 1, cuts=0, end=0, cyc_start=-1,
                cyc_len=-1, width=-1, ntrim=0, s2start=-1, s2step=-1,
                delta=0.5, msum=0)
