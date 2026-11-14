# -*- coding: utf-8 -*-
"""k_engine -- ONE colour model, k a real parameter: BTHM-C_k.

THE ORACLE, once, cleanly
-------------------------
The oracle answers the VISIT COUNT of the current quadruple (ell, r, L, R)
EXACTLY: on arrival it reports how many times that quadruple has already been
visited in this run (0 on the first visit, 1 on the second, ...).

THE RULE
--------
The rule's action table reads that count MODULO k: colour c = count mod k, so
a base rule has k colour-variants (c = 0..k-1).  Visits 1..k of a quadruple
are PLAYED by the rule with colours 0,1,...,k-1.

THE REPAIR
----------
The repair fires when the count is positive and the count mod k == 0, i.e. on
the (k+1)-th visit, so that k intermediate visits are played by the rule
before any trim.  The repair is BTHM-C's cut repair (cut_engine mode 4): the
boundary cell of the tape goes (L <- L+1 or R <- R-1) on the side the
COLOUR-0 (base) table selects for that state, with the re-span charge; the
key (L, R) then changes and never returns (L only rises, R only falls), so
each quadruple fires at most one repair.

k = 1  is the colour-blind model (colour pinned to 0, repair on the 2nd
       visit): it must reproduce probe_trim/cut_cells.json bit for bit.
k = 4  reproduces probe_colour/col_engine.py (BTHM-C_4 of COLOUR.md).

THE VARIANT DESIGN (fixed across k)
-----------------------------------
At a played step on a state with colour c the rule computes its base action
a = TAB[2*b2 + b1] and maps it:
    c odd        : a <- mirror(a):  LPLUS<->RMINUS, LMINUS<->RPLUS;
                   the three trims (CUTL, CUTR, CUTB) and HALT are fixed;
    c is 2 or 3  : the side of a cut is swapped: CUTL<->CUTR (CUTB, HALT
                   and the four moves fixed).
For k = 2: only c = 0,1 occur -- only the mirror acts.
For k = 3: c = 0,1,2 -- the mirror at c = 1, the swap at c = 2.
For k = 4: c = 0,1,2,3 -- mirror at 1 and 3, swap at 2 and 3 (COLOUR.md).
A move that would leave the region is replaced by a cut; the side is part of
the rule: the cut falls on the side the (mapped) move would leave for
c in {0,1}, on the OPPOSITE side for c in {2,3}.  The double cut used when
the window fills the region (ell == L and r == R_) is kept for every colour.
The repair itself reads the side off the colour-0 variant (= base table).

Modes of this file: 0 = BTHM-0 (no memory, die at the first repeated
quadruple, end=4), 5 = the colour model.  The BTHM-0 mode is a copy of
col_engine's mode 0 and is checked against the stored cut files tape for tape.

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
def walkk(rid, conv, a, Pa, Paa, PR, Pd, Pn2, Ppk, Ptr, Ps, cpos, pnz,
          C2, bstart, bpos, psmin, SPs, SNs, LOG2, r2s, r2e, SP2, rcs, rce,
          rcodd, SPC, R, IP, HX, Ppat, d1, d2, n, cap, mode, k,
          hkeys, hvals, hmask, genbox, org):
    """One tape, modes 0 (BTHM-0) and 5 (colour model, k colour-variants).

    org[0..14] as col_engine.walkv; 15..18 = steps played on colour 0/1/2/3;
    org[19] = repairs whose stored count was exactly k (must equal ntrim);
    org[20] = step of the first repeat seen in this run (0 if none);
    org[21] = repairs whose stored count was > 0 and == 0 mod k.
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
    rmod = 0
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
        # ---- the oracle: the EXACT visit count of this quadruple
        key = ((ell * (n + 1) + r) * (n + 1) + L) * (n + 1) + R_
        h = key * 0x9E3779B97F4A7C15
        h = (h ^ (h >> 29)) * 0xBF58476D1CE4E5B9
        idx = (h ^ (h >> 32)) & hmask
        cnt = -1
        while hkeys[idx] != -1:
            if hkeys[idx] == key:
                cnt = hvals[idx]
                break
            idx = (idx + 1) & hmask
        colour = 0
        if cnt >= 0:
            if mode == 0:
                cstart = cnt
                clen = steps - cnt
                cwidth = m
                end = 4
                break
            if fstep == 0:
                fstep = steps
            if cnt > 0 and cnt % k == 0:
                # ---- THE REPAIR, on the (k+1)-th visit (colour 0 again):
                # identical to cut_engine mode 4: the side is read off the
                # colour-0 variant (= the base table), the tape loses one
                # boundary cell.  The visit count is kept exact.
                ntrim += 1
                if cnt == k:
                    repok += 1
                rmod += 1
                hvals[idx] = cnt + 1
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
                kk = (x - y + 1) >> 1
                if kk == 0:
                    kk = 1
                ty = y
                tx = y + kk - 1
                reads += abs(hy - ty) + abs(hx - tx)
                hy = ty
                hx = tx
                continue
            colour = cnt % k              # the colour the rule reads
            hvals[idx] = cnt + 1
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
        if colour == 2 or colour == 3:
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
            kk = (x - y + 1) >> 1
            if kk == 0:
                kk = 1
            ty = y
            tx = y + kk - 1
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
            kk = (x - y + 1) >> 1
            if kk == 0:
                kk = 1
            ty = y
            tx = y + kk - 1
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
        kk = (x - y + 1) >> 1
        if kk == 0:
            kk = 1
        ty = y
        tx = y + kk - 1
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
    org[21] = rmod
    return


class Col(object):
    """One tape's structures plus the visit-count table, reused across runs."""

    def __init__(self, X, hlog=HLOG):
        self.S = P.St(X)
        self.n = self.S.n
        self.hkeys = np.full(1 << hlog, -1, np.int64)
        self.hvals = np.zeros(1 << hlog, np.int32)
        self.hmask = np.int64((1 << hlog) - 1)
        self.genbox = np.zeros(1, np.int64)
        self.org = np.zeros(24, np.float64)

    def run(self, rid, conv=0, mode=5, k=4, cap=None):
        S = self.S
        if cap is None:
            cap = 40 * self.n * self.n + 100
        self.hkeys.fill(-1)
        self.hvals.fill(0)
        self.genbox[0] = 0
        walkk(rid, conv, S.a, S.Pa, S.Paa, S.PR, S.Pd, S.Pn2, S.Ppk,
              S.Ptr, S.Ps, S.cpos, S.pnz, S.C2, S.bstart, S.bpos, S.psmin,
              S.SPs, S.SNs, S.LOG2, S.r2s, S.r2e, S.SP2, S.rcs, S.rce,
              S.rcodd, S.SPC, S.R, S.IP, S.HX, S.Ppat, S.d1, S.d2, self.n,
              cap, mode, k, self.hkeys, self.hvals, self.hmask,
              self.genbox, self.org)
        o = self.org
        return dict(steps=int(o[0]), reads=int(o[1]), cuts=int(o[2]),
                    end=int(o[3]), cyc_start=int(o[4]), cyc_len=int(o[5]),
                    width=int(o[6]), ntrim=int(o[7]),
                    delta=float(o[10]) if o[10] != -1.0 else None,
                    msum=int(o[11]), nl=int(o[12]), nr=int(o[13]),
                    nused=int(o[14]), kc0=int(o[15]), kc1=int(o[16]),
                    kc2=int(o[17]), kc3=int(o[18]), repok=int(o[19]),
                    fstep=int(o[20]), rmod=int(o[21]))
