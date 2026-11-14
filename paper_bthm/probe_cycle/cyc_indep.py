# -*- coding: utf-8 -*-
"""cyc_indep -- a second opinion on the cycle and trim logic, in plain Python.

The numba walk is validated against tr_lib for (delta, end); this script
re-derives the trajectory with a dict-based seen-set and its own trim loop, and
checks the four numbers the table quotes: end0, steps0, cyc_len, width, and
(for BTHM-C) endC, ntrim, readsC.  Slow, so a handful of tapes only.
"""
import sys

import numpy as np

import cyc_engine as E
import cyc_verify as V
import tr_lib as T
import pcr_lib as P

TAB = P.TABLE_NP
LPLUS, LMINUS, RPLUS, RMINUS, CUTL, CUTR, CUTB, HALT = range(8)


def b1(rule, S, ell, r):
    if rule == 'BAL>m/4':
        import tr_bal4
        return int(tr_bal4.bal_b1(S.a, S.Pa, S.Ps, S.psmin, S.bstart, S.bpos, 4, ell, r))
    if rule == 'NBLK>3':
        Ps, psmin = S.Ps, S.psmin
        b = Ps[ell] - psmin
        j0, j1 = S.bstart[b], S.bstart[b + 1]
        idx = [j for j in range(j0, j1) if ell + 1 <= S.bpos[j] <= r + 1]
        return 1 if len(idx) > 3 else 0
    rid = V.rrid(rule)
    return int(P.b1_of(rid, S.a, S.Pa, S.Paa, S.PR, S.Pd, S.Pn2, S.Ppk, S.Ptr,
                       S.Ps, S.cpos, S.pnz, S.C2, S.bstart, S.bpos, S.psmin,
                       S.SPs, S.SNs, S.LOG2, S.r2s, S.r2e, S.SP2, S.rcs, S.rce,
                       S.rcodd, S.SPC, S.R, S.IP, S.HX, S.Ppat, S.d1, S.d2,
                       ell, r))


def trace(rule, S, n, cap, trim=None, maxsteps=4_000_000):
    """mode 0 (trim=None): stop at the first repeated quadruple.
    mode 1 (trim=1/2): oracle + trim, stop at a terminal state."""
    L, R = 0, n - 1
    ell, r = L, R
    y, x = L, R
    hy, hx = 0, max(1, n // 2) - 1
    steps = ntrim = cuts = 0
    reads = n + 2 * (n - 1) + 2 + max(1, n // 2)
    seen = {}
    cyc = None
    while True:
        if ell >= r:
            return dict(end=0, steps=steps, reads=reads, ntrim=ntrim,
                        cyc_len=(cyc[0] if cyc else None),
                        width=(cyc[1] if cyc else None))
        if R - L + 1 < 3:
            return dict(end=1, steps=steps, reads=reads, ntrim=ntrim,
                        cyc_len=(cyc[0] if cyc else None),
                        width=(cyc[1] if cyc else None))
        steps += 1
        if steps + ntrim > cap or steps > maxsteps:
            return dict(end=3, steps=steps, reads=reads, ntrim=ntrim,
                        cyc_len=(cyc[0] if cyc else None),
                        width=(cyc[1] if cyc else None))
        key = (ell, r, L, R)
        m = r - ell + 1
        if key in seen:
            if trim is None:
                return dict(end=4, steps=steps, reads=reads, ntrim=0,
                            cyc_len=steps - seen[key], width=m)
            if cyc is None:
                cyc = (steps - seen[key], m)
            ntrim += 1
            if trim == 1:
                ell += 1
            else:
                r -= 1
            y, x = ell, r
            reads += 1
            k = max(1, (x - y + 1) >> 1)
            ty, tx = y, y + k - 1
            reads += abs(hy - ty) + abs(hx - tx)
            hy, hx = ty, tx
            continue
        seen[key] = steps
        v1 = b1(rule, S, ell, r)
        if rule == 'hashmod':
            act = CUTL if v1 == 1 else RMINUS
        else:
            v2 = int(P.b2_of(S.a, S.Pa, ell, r))
            act = int(TAB[2 * v2 + v1])
        if act >= CUTL:
            cuts += 1
            if act == CUTL:
                L += 1
            elif act == CUTR:
                R -= 1
            else:
                L += 1
                R -= 1
            if L > R:
                return dict(end=1, steps=steps, reads=reads, ntrim=ntrim,
                            cyc_len=(cyc[0] if cyc else None),
                            width=(cyc[1] if cyc else None))
            ell, r = L, R
            reads += abs(y - L) + abs(x - R)
            y, x = L, R
            k = max(1, (x - y + 1) >> 1)
            ty, tx = y, y + k - 1
            reads += abs(hy - ty) + abs(hx - tx)
            hy, hx = ty, tx
            continue
        dy = 1 if act == LPLUS else (-1 if act == LMINUS else 0)
        dx = 1 if act == RPLUS else (-1 if act == RMINUS else 0)
        if ell + dy < L or r + dx > R:
            cuts += 1
            if V.conv(rule) == 1:
                if ell + dy < L:
                    R -= 1
                else:
                    L += 1
            else:
                if ell == L and r == R:
                    L += 1
                    R -= 1
                    y += 1
                    x -= 1
                    reads += 2
                elif ell + dy < L:
                    L += 1
                else:
                    R -= 1
            if L > R:
                return dict(end=1, steps=steps, reads=reads, ntrim=ntrim,
                            cyc_len=(cyc[0] if cyc else None),
                            width=(cyc[1] if cyc else None))
            ell, r = L, R
            reads += abs(y - L) + abs(x - R)
            y, x = L, R
            k = max(1, (x - y + 1) >> 1)
            ty, tx = y, y + k - 1
            reads += abs(hy - ty) + abs(hx - tx)
            hy, hx = ty, tx
            continue
        if ell + dy >= r + dx:
            return dict(end=0, steps=steps, reads=reads, ntrim=ntrim,
                        cyc_len=(cyc[0] if cyc else None),
                        width=(cyc[1] if cyc else None))
        ell += dy
        r += dx
        y, x = ell, r
        reads += 1
        k = max(1, (x - y + 1) >> 1)
        ty, tx = y, y + k - 1
        reads += abs(hy - ty) + abs(hx - tx)
        hy, hx = ty, tx


def main():
    n = 512
    rules = ['canonical', 'dense2', 'pattern', 'BAL>m/2', 'NBLK>3', 'AMAX',
             'longrun', 'pal', 'hashmod']
    bad = 0
    tot = 0
    for rule in rules:
        for i, a, s, X in T.tapes(n, 4):
            S = P.St(X)
            ref0 = trace(rule, S, n, 40 * n * n + 100, trim=None)
            C = E.Cyc(X)
            got0 = C.run(V.rrid(rule), V.conv(rule), V.variant(rule), mode=0)
            refC = trace(rule, S, n, 40 * n * n + 100, trim=1)
            gotC = C.run(V.rrid(rule), V.conv(rule), V.variant(rule), mode=1, trim=1)
            tot += 1
            checks = [
                ('end0', ref0['end'], got0['end']),
                ('steps0', ref0['steps'], got0['steps']),
                ('cyc_len', ref0['cyc_len'], (got0['cyc_len'] if got0['end'] == 4 else None)),
                ('width', ref0['width'], (got0['width'] if got0['end'] == 4 else None)),
                ('endC', refC['end'], gotC['end']),
                ('ntrim', refC['ntrim'], gotC['ntrim']),
                ('readsC', refC['reads'], gotC['reads']),
            ]
            for name, a_, b_ in checks:
                if a_ != b_:
                    bad += 1
                    print('DIFF %-9s tape %d %-8s python %s numba %s'
                          % (rule, i, name, a_, b_), flush=True)
    print('indep check n=%d: %d/%d tape-runs agree on all 7 fields'
          % (n, tot - bad, tot))
    return 0 if bad == 0 else 1


if __name__ == '__main__':
    sys.exit(main())
