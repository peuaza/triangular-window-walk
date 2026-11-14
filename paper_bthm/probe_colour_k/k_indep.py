# -*- coding: utf-8 -*-
"""k_indep -- a second implementation of the BTHM-C_k walk, in plain Python.

Written from the model text alone (k_engine's docstring), NOT from k_engine's
code: a dict holds the EXACT visit count of every quadruple (ell, r, L, R);
on a visit with count cnt the rule reads the colour c = cnt mod k and plays
the variant; the repair fires when cnt > 0 and cnt mod k == 0 (the (k+1)-th
visit) and removes one boundary cell of the tape on the side the colour-0
(base) table selects, with the re-span charge.  Predicates are pcr_lib's
b1_of / b2_of (not under test); the charge is pcr_lib's INC accounting,
copied literally, as in k_engine.

usage: k_indep.py [ntapes512] [ntapes1024]   writes k_indep.json
"""
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import k_engine as E                                          # noqa: E402
import pcr_lib as P                                           # noqa: E402
import tr_lib as T                                            # noqa: E402

MIR = (3, 2, 1, 0, 4, 5, 6, 7)          # LPLUS<->RMINUS, LMINUS<->RPLUS
SWAP = (0, 1, 2, 3, 5, 4, 6, 7)         # CUTL<->CUTR
TABv = [int(x) for x in P.TABLE_NP]
RULES = ('pattern', 'runs', 'peak')


def walk_indep(S, rid, k, cap=None, conv=0):
    """The model: exact visit count; colour = count mod k; repair on the
    (k+1)-th visit (count positive and count mod k == 0)."""
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
        q = (x - y + 1) >> 1
        if q == 0:
            q = 1
        ty, tx = y, y + q - 1
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
        if cnt > 0 and cnt % k == 0:
            # ---- the repair, on the (k+1)-th visit: the boundary cell of
            # the tape goes, on the side the colour-0 (base) table selects
            ntrim += 1
            seen[key] = cnt + 1
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
        colour = cnt % k                    # the colour the rule reads
        seen[key] = cnt + 1
        if cnt == 0:
            pass                            # first visit: colour 0
        v1, v2 = b1b2(ell, r)
        act = TABv[2 * v2 + v1]
        if colour & 1:
            act = MIR[act]
        if colour == 2 or colour == 3:
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
        q = (x - y + 1) >> 1
        if q == 0:
            q = 1
        ty, tx = y, y + q - 1
        reads += abs(hy - ty) + abs(hx - tx)
        hy, hx = ty, tx
    return dict(steps=steps, reads=reads, cuts=cuts, end=end, ntrim=ntrim,
                nl=nl, nr=nr, out=out, msum=msum)


def main():
    nt1 = int(sys.argv[1]) if len(sys.argv) > 1 else 12
    nt2 = int(sys.argv[2]) if len(sys.argv) > 2 else 4
    bad = []
    nchk = 0
    runs = []
    for n, nt in ((512, nt1), (1024, nt2)):
        for rule in RULES:
            rid = P.RID[rule]
            for i, a, s, X in T.tapes(n, nt):
                S = P.St(X)
                C = E.Col(X)
                for k in (1, 2, 3, 4):
                    r1 = walk_indep(S, rid, k)
                    r2 = C.run(rid, 0, mode=5, k=k)
                    nchk += 1
                    f = [key for key in ('end', 'steps', 'reads', 'ntrim',
                                         'nl', 'nr')
                         if int(r1[key]) != int(r2[key])]
                    runs.append(dict(n=n, rule=rule, i=i, k=k, fields=f))
                    if f:
                        bad.append((rule, n, i, k, f,
                                    [r1[key] for key in f],
                                    [r2[key] for key in f]))
            print('%-8s n=%-5d tapes %d checked, k 1-4' % (rule, n, nt),
                  flush=True)
    out = dict(ns=[512, 1024], tapes512=nt1, tapes1024=nt2,
               runs=nchk, bad=bad)
    json.dump(out, open(os.path.join(HERE, 'k_indep.json'), 'w'), indent=1)
    print('k_indep: %d runs, %d mismatches' % (nchk, len(bad)), flush=True)
    for b in bad[:10]:
        print('MISMATCH', b)
    return 0 if not bad else 1


if __name__ == '__main__':
    sys.exit(main())
