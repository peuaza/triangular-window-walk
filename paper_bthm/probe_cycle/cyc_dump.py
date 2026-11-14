# -*- coding: utf-8 -*-
"""cyc_dump -- record every state (ell, r, L, R) the numba walk visits.

A copy of cyc_engine.walk's BTHM-0 path with recording, used only to locate the
first state at which two engines disagree.
"""
import numpy as np
from numba import njit

from pcr_lib import b1_of, b2_of
import pcr_lib as P

LPLUS, LMINUS, RPLUS, RMINUS, CUTL, CUTR, CUTB, HALT = range(8)
TAB = P.TABLE_NP


@njit(cache=True)
def dump(variant, rid, conv, a, Pa, Paa, PR, Pd, Pn2, Ppk, Ptr, Ps, cpos, pnz,
         C2, bstart, bpos, psmin, SPs, SNs, LOG2, r2s, r2e, SP2, rcs, rce,
         rcodd, SPC, R, IP, HX, Ppat, d1, d2, n, cap, traj, want_b1):
    """Fill traj[k] = (ell, r, L, R) and want_b1[k] = b1 at that state."""
    L = 0
    R_ = n - 1
    ell = L
    r = R_
    y = L
    x = R_
    hy = 0
    hx = max(1, n // 2) - 1
    steps = 0
    nrec = 0
    end = 1
    out = -1.0
    while True:
        if ell >= r:
            end = 0
            out = ell / (n - 1.0)
            break
        if R_ - L + 1 < 3:
            end = 1
            break
        steps += 1
        if steps > cap:
            end = 3
            break
        if nrec < traj.shape[0]:
            traj[nrec, 0] = ell
            traj[nrec, 1] = r
            traj[nrec, 2] = L
            traj[nrec, 3] = R_
        b1 = b1_of(rid, a, Pa, Paa, PR, Pd, Pn2, Ppk, Ptr, Ps, cpos, pnz, C2,
                   bstart, bpos, psmin, SPs, SNs, LOG2, r2s, r2e, SP2, rcs, rce,
                   rcodd, SPC, R, IP, HX, Ppat, d1, d2, ell, r)
        if nrec < want_b1.shape[0]:
            want_b1[nrec] = b1
        nrec += 1
        if variant == 3:
            act = CUTL if b1 == 1 else RMINUS
        else:
            b2 = b2_of(a, Pa, ell, r)
            act = TAB[2 * b2 + b1]
        if act >= CUTL:
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
            y = L
            x = R_
            k = (x - y + 1) >> 1
            if k == 0:
                k = 1
            hy = y
            hx = y + k - 1
            continue
        dy = 1 if act == LPLUS else (-1 if act == LMINUS else 0)
        dx = 1 if act == RPLUS else (-1 if act == RMINUS else 0)
        if ell + dy < L or r + dx > R_:
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
                elif ell + dy < L:
                    L += 1
                else:
                    R_ -= 1
            if L > R_:
                end = 1
                break
            ell = L
            r = R_
            y = L
            x = R_
            k = (x - y + 1) >> 1
            if k == 0:
                k = 1
            hy = y
            hx = y + k - 1
            continue
        if ell + dy >= r + dx:
            out = (ell + dy) / (n - 1.0)
            end = 0
            break
        ell += dy
        r += dx
        y = ell
        x = r
        k = (x - y + 1) >> 1
        if k == 0:
            k = 1
        hy = y
        hx = y + k - 1
    return steps, end, out, nrec
