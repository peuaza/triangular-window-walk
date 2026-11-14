# -*- coding: utf-8 -*-
"""col_engine -- BTHM-C_k: the oracle answers the visit count of the state
modulo k, and the rule may read the colour.  k = 4 here (BTHM-C_4).

The walk is cut_engine.walk (mode 4, the cut repair) with the oracle changed in
ONE place: instead of one bit, the oracle table stores the VISIT COUNT of the
quadruple (ell, r, L, R), capped at k.  Colour c = (visit count - 1) mod k, so
colour 0 is the first visit, 1 the second, 2 the third, 3 the fourth; on the
(k+1)-th visit the colour returns to 0 and the REPAIR fires -- the same repair
BTHM-C uses (remove the boundary cell of the TAPE on the side the rule's
colour-0 table selects for that state), with the same charge (the re-span a
rule-initiated cut pays, abs(y-L)+abs(r... x-R_) plus the head re-read, no
oracle read).  Visits 1..k are actually PLAYED by the rule, which reads the
colour.

THE FOUR COLOUR-VARIANTS (design; written out in full in COLOUR.md)
  At a played step on a state with colour c the rule computes its base action
  a = TAB[2*b2 + b1] and maps it:
    c odd      : a <- mirror(a):  LPLUS<->RMINUS, LMINUS<->RPLUS;
                 the three trims (CUTL, CUTR, CUTB) and HALT are fixed;
    c in {2,3} : the side of a cut is swapped: CUTL<->CUTR (CUTB, HALT and the
                 four moves fixed; the tables of the measured rules carry no
                 CUTL/CUTR entry, so this clause acts on the boundary cut
                 below).
  The mapped action is then executed exactly as the model executes any action.
  A move that would leave the region is replaced by a cut, and that side is
  part of the RULE (pcr_lib: the four-move rules replace the move by a cut on
  the side it would have left).  The variant inherits this: for c in {0,1} the
  cut falls on the side the (mapped) move would leave, for c in {2,3} on the
  OPPOSITE side.  The double cut the model uses when the window fills the
  region (ell == L and r == R_) is kept for every colour.
  The repair itself (the (k+1)-th visit) reads the side off the colour-0
  variant, i.e. the base table -- identical to cut_engine mode 4.

kmax: the repair fires at the visit whose count reaches kmax (= k).  kmax = 1
with the colour pinned to 0 reproduces BTHM-C (cut_engine mode 4) tape for
tape; that run is the control.

Modes of this file: 0 = BTHM-0 (no memory, die at the first repeated
quadruple, end=4), 5 = the colour model.  The BTHM-0 mode is a copy of
cut_engine's mode 0 and is checked against cut_engine tape for tape by
col_check.py.

END CODES (unchanged): 0 meet, 1 region empty, 3 cap, 4 quadruple (BTHM-0),
5 table full, 6 cap after at least one repair.
"""
import os
import sys

import numpy as np
from numba import njit

HERE = os.path.dirname(os.path.abspath(__file__))
PAPER = os.path.dirname(HERE)
for _d in ('probe_trim', 'probe_class_rerun', 'probe_e_uniform',
           'probe_triple_nat'):
    _p = os.path.join(PAPER, _d)
    if _p not in sys.path:
        sys.path.insert(0, _p)

import pcr_lib as P                                  # noqa: E402
from pcr_lib import b1_of, b2_of                     # noqa: E402

LPLUS, LMINUS, RPLUS, RMINUS, CUTL, CUTR, CUTB, HALT = range(8)
TAB = P.TABLE_NP
MIR = (3, 2, 1, 0, 4, 5, 6, 7)          # LPLUS<->RMINUS, LMINUS<->RPLUS
SWAPCUT = (0, 1, 2, 3, 5, 4, 6, 7)      # CUTL<->CUTR

HLOG = 24


@njit(cache=True)
def walkv(rid, conv, a, Pa, Paa, PR, Pd, Pn2, Ppk, Ptr, Ps, cpos, pnz,
          C2, bstart, bpos, psmin, SPs, SNs, LOG2, r2s, r2e, SP2, rcs, rce,
          rcodd, SPC, R, IP, HX, Ppat, d1, d2, n, cap, mode, kmax,
          hkeys, hvals, hmask, genbox, org):
    """One tape, modes 0 (BTHM-0) and 5 (colour model, kmax visits played).

    org[0..14] as cut_engine.walk; 15..18 = steps played on colour 0/1/2/3;
    org[19] = repairs whose stored count was exactly kmax (must equal ntrim);
    org[20] = step of the first repeat seen in this run (0 if none).
    """
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
    nl = 0
    nr = 0
    nused = 0
    kc0 = 0
    kc1 = 0
    kc2 = 0
    kc3 = 0
    repok = 0
    fstep = 0
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
        m = r - ell + 1
        msum += m
        # ---- the oracle: the visit count of this quadruple, capped at kmax
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
        colour = 0
        if first >= 0:
            if mode == 0:
                cstart = first
                clen = steps - first
                cwidth = m
                end = 4
                break
            if fstep == 0:
                fstep = steps
            if first >= kmax:
                # ---- THE REPAIR, on the (kmax+1)-th visit (colour 0 again)
                # identical to cut_engine mode 4: side read off the colour-0
                # variant (= the base table), tape loses one boundary cell
                ntrim += 1
                if first == kmax:
                    repok += 1
                v1 = b1_of(rid, a, Pa, Paa, PR, Pd, Pn2, Ppk, Ptr, Ps, cpos,
                           pnz, C2, bstart, bpos, psmin, SPs, SNs, LOG2, r2s,
                           r2e, SP2, rcs, rce, rcodd, SPC, R, IP, HX, Ppat,
                           d1, d2, ell, r)
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
            colour = first               # 0..kmax-1, the colour the rule reads
            hvals[idx] = first + 1
        else:
            hkeys[idx] = key
            if mode == 0:
                hvals[idx] = steps        # mode 0 stores the step index
            else:
                hvals[idx] = 1            # the colour model stores a count
            nused += 1
            if nused > (hmask + 1) * 3 // 4:
                end = 5
                break
            colour = 0
        # ---- the played step: the rule reads the colour
        b1 = b1_of(rid, a, Pa, Paa, PR, Pd, Pn2, Ppk, Ptr, Ps, cpos, pnz, C2,
                   bstart, bpos, psmin, SPs, SNs, LOG2, r2s, r2e, SP2, rcs,
                   rce, rcodd, SPC, R, IP, HX, Ppat, d1, d2, ell, r)
        b2 = b2_of(a, Pa, ell, r)
        act = TAB[2 * b2 + b1]
        if colour == 0:
            kc0 += 1
        elif colour == 1:
            kc1 += 1
        elif colour == 2:
            kc2 += 1
        else:
            kc3 += 1
        if colour & 1:
            act = MIR[act]
        if colour >= 2:
            act = SWAPCUT[act]
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
                else:
                    left = ell + dy < L
                    if (colour < 2) == left:
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
    org[8] = -1
    org[9] = -1
    org[10] = out
    org[11] = msum
    org[12] = nl
    org[13] = nr
    org[14] = nused
    org[15] = kc0
    org[16] = kc1
    org[17] = kc2
    org[18] = kc3
    org[19] = repok
    org[20] = fstep
    return


class Col(object):
    """One tape's structures plus the count table, reused across runs."""

    def __init__(self, X, hlog=HLOG):
        self.S = P.St(X)
        self.n = self.S.n
        self.hkeys = np.full(1 << hlog, -1, np.int64)
        self.hvals = np.zeros(1 << hlog, np.int32)
        self.hmask = np.int64((1 << hlog) - 1)
        self.genbox = np.zeros(1, np.int64)
        self.org = np.zeros(24, np.float64)

    def run(self, rid, conv=0, mode=5, kmax=4, cap=None):
        S = self.S
        if cap is None:
            cap = 40 * self.n * self.n + 100
        self.hkeys.fill(-1)
        self.hvals.fill(0)
        self.genbox[0] = 0
        walkv(rid, conv, S.a, S.Pa, S.Paa, S.PR, S.Pd, S.Pn2, S.Ppk,
              S.Ptr, S.Ps, S.cpos, S.pnz, S.C2, S.bstart, S.bpos, S.psmin,
              S.SPs, S.SNs, S.LOG2, S.r2s, S.r2e, S.SP2, S.rcs, S.rce,
              S.rcodd, S.SPC, S.R, S.IP, S.HX, S.Ppat, S.d1, S.d2, self.n,
              cap, mode, kmax, self.hkeys, self.hvals, self.hmask,
              self.genbox, self.org)
        o = self.org
        return dict(steps=int(o[0]), reads=int(o[1]), cuts=int(o[2]),
                    end=int(o[3]), cyc_start=int(o[4]), cyc_len=int(o[5]),
                    width=int(o[6]), ntrim=int(o[7]),
                    delta=float(o[10]) if o[10] != -1.0 else None,
                    msum=int(o[11]), nl=int(o[12]), nr=int(o[13]),
                    nused=int(o[14]), kc0=int(o[15]), kc1=int(o[16]),
                    kc2=int(o[17]), kc3=int(o[18]), repok=int(o[19]),
                    fstep=int(o[20]))
