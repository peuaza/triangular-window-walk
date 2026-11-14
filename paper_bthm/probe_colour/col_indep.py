# -*- coding: utf-8 -*-
"""col_indep -- a second implementation of the colour walk, in plain Python.

Written from the design text (COLOUR.md) alone: a dict holds the visit count of
every quadruple (ell, r, L, R); visits 0..kmax-1 are played with the colour
c = count, the count kmax fires the repair (cut_engine mode 4's repair: the
boundary cell of the tape goes, on the side the colour-0 table selects).
Predicates are pcr_lib.b1_of / b2_of -- not under test.

Compared with col_engine.Col on the same tapes: kmax = 1 (BTHM-C) and kmax = 4
(BTHM-C_4), fields (end, steps, reads, ntrim, nl, nr).

usage: col_indep.py [ntapes]        writes col_indep.json
"""
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import col_engine as E                                          # noqa: E402
import pcr_lib as P                                             # noqa: E402
import tr_lib as T                                              # noqa: E402

MIR = (3, 2, 1, 0, 4, 5, 6, 7)
SWAP = (0, 1, 2, 3, 5, 4, 6, 7)
TABv = [int(x) for x in P.TABLE_NP]          # [LMINUS, RPLUS, RMINUS, LPLUS]
RULES = ('pattern', 'runs', 'peak')


def walk_indep(S, rid, kmax, cap=None, conv=0):
    n = S.n
    a = S.a
    if cap is None:
        cap = 40 * n * n + 100
    seen = {}
    L, R_ = 0, n - 1
    ell, r = L, R_
    y, x = L, R_
    hy, hx = 0, max(1, n // 2) - 1
    steps = ntrim = cuts = nl = nr = 0
    msum = 0
    reads = n + 2 * (n - 1) + 2 + max(1, n // 2)
    out = None
    end = 1

    def b1b2(ell, r):
        v1 = int(P.b1_of(rid, a, S.Pa, S.Paa, S.PR, S.Pd, S.Pn2, S.Ppk,
                         S.Ptr, S.Ps, S.cpos, S.pnz, S.C2, S.bstart, S.bpos,
                         S.psmin, S.SPs, S.SNs, S.LOG2, S.r2s, S.r2e, S.SP2,
                         S.rcs, S.rce, S.rcodd, S.SPC, S.R, S.IP, S.HX,
                         S.Ppat, S.d1, S.d2, ell, r))
        v2 = int(P.b2_of(a, S.Pa, ell, r))
        return v1, v2

    def respan(L, R_, y, x, hy, hx):
        rd = abs(y - L) + abs(x - R_)
        y, x = L, R_
        k = (x - y + 1) >> 1
        if k == 0:
            k = 1
        ty, tx = y, y + k - 1
        rd += abs(hy - ty) + abs(hx - tx)
        return rd, y, x, ty, tx

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
        key = (ell, r, L, R_)
        cnt = seen.get(key, 0)
        if cnt >= kmax:
            # ---- the repair, on the (kmax+1)-th visit: the boundary cell of
            # the tape goes, on the side the colour-0 (base) table selects
            ntrim += 1
            v1, v2 = b1b2(ell, r)
            act = TABv[2 * v2 + v1]
            side = 1 if act in (0, 1, 4) else 2
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
            ell, r = L, R_
            rd, y, x, hy, hx = respan(L, R_, y, x, hy, hx)
            reads += rd
            continue
        colour = cnt
        seen[key] = cnt + 1
        v1, v2 = b1b2(ell, r)
        act = TABv[2 * v2 + v1]
        if colour & 1:
            act = MIR[act]
        if colour >= 2:
            act = SWAP[act]
        if act >= 4:
            cuts += 1
            if act == 4:
                L += 1
            elif act == 5:
                R_ -= 1
            else:
                L += 1
                R_ -= 1
            if L > R_:
                end = 1
                break
            ell, r = L, R_
            rd, y, x, hy, hx = respan(L, R_, y, x, hy, hx)
            reads += rd
            continue
        dy = 1 if act == 0 else (-1 if act == 1 else 0)
        dx = 1 if act == 2 else (-1 if act == 3 else 0)
        if ell + dy < L or r + dx > R_:
            cuts += 1
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
            ell, r = L, R_
            rd, y, x, hy, hx = respan(L, R_, y, x, hy, hx)
            reads += rd
            continue
        if ell + dy >= r + dx:
            out = (ell + dy) / (n - 1.0)
            end = 0
            break
        ell += dy
        r += dx
        y, x = ell, r
        reads += 1
        k = (x - y + 1) >> 1
        if k == 0:
            k = 1
        ty, tx = y, y + k - 1
        reads += abs(hy - ty) + abs(hx - tx)
        hy, hx = ty, tx
    return dict(steps=steps, reads=reads, cuts=cuts, end=end, ntrim=ntrim,
                nl=nl, nr=nr, out=out, msum=msum)


def main():
    nt = int(sys.argv[1]) if len(sys.argv) > 1 else 8
    bad = []
    nchk = 0
    for rule in RULES:
        rid = P.RID[rule]
        for i, a, s, X in T.tapes(512, nt):
            S = P.St(X)
            C = E.Col(X)
            for kmax in (1, 4):
                r1 = walk_indep(S, rid, kmax)
                r2 = C.run(rid, 0, mode=5, kmax=kmax)
                nchk += 1
                f = [k for k in ('end', 'steps', 'reads', 'ntrim', 'nl', 'nr')
                     if int(r1[k]) != int(r2[k])]
                if f:
                    bad.append((rule, i, kmax, f,
                                [r1[k] for k in f], [r2[k] for k in f]))
            print('%-8s tape %-3d checked (kmax 1 and 4)' % (rule, i),
                  flush=True)
    out = dict(n=512, tapes=nt, rules=list(RULES), runs=nchk, bad=bad)
    json.dump(out, open(os.path.join(HERE, 'col_indep.json'), 'w'), indent=1)
    print('col_indep: %d runs, %d mismatches' % (nchk, len(bad)), flush=True)
    for b in bad[:10]:
        print('MISMATCH', b)
    return 0 if not bad else 1


if __name__ == '__main__':
    sys.exit(main())
