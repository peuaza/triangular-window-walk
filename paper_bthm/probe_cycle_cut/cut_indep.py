# -*- coding: utf-8 -*-
"""cut_indep -- a second opinion on the CUT repair, in plain Python.

The numba branch (cut_engine.py mode 4) is checked here against a dict-based
trace with its own repair loop: on a repeated quadruple it reads the rule's
action from TABLE[2*b2+b1], takes the side from it (c1a inverted), removes that
boundary cell of the tape, re-spans the window and continues.  Compares
end, steps, reads, ntrim, nl, nr.

Slow, so a handful of tapes per rule, n = 512.  Writes cut_indep.json.
"""
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import cut_engine as E                                           # noqa: E402
import cut_verify as V                                           # noqa: E402
import tr_lib as T                                               # noqa: E402
import pcr_lib as P                                              # noqa: E402

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
    return int(P.b1_of(rid, S.a, S.Pa, S.Paa, S.PR, S.Pd, S.Pn2, S.Ppk,
                       S.Ptr, S.Ps, S.cpos, S.pnz, S.C2, S.bstart, S.bpos, S.psmin,
                       S.SPs, S.SNs, S.LOG2, S.r2s, S.r2e, S.SP2, S.rcs, S.rce,
                       S.rcodd, S.SPC, S.R, S.IP, S.HX, S.Ppat, S.d1, S.d2,
                       ell, r))


def act_of(rule, S, ell, r):
    v1 = b1(rule, S, ell, r)
    if V.variant(rule) == 3:
        return CUTL if v1 == 1 else RMINUS
    v2 = int(P.b2_of(S.a, S.Pa, ell, r))
    return int(TAB[2 * v2 + v1])


def trace(rule, S, n, cap, repair=False, maxsteps=40 * 512 * 512 + 100):
    """repair=False: mode 0 (stop at the first repeated quadruple).
    repair=True:  mode 4 (cut the boundary cell of the tape on the rule's side)."""
    L, R = 0, n - 1
    ell, r = L, R
    y, x = L, R
    hy, hx = 0, max(1, n // 2) - 1
    steps = ntrim = cuts = nl = nr = 0
    reads = n + 2 * (n - 1) + 2 + max(1, n // 2)
    seen = {}
    cyc = None
    while True:
        if ell >= r:
            return dict(end=0, steps=steps, reads=reads, ntrim=ntrim, nl=nl, nr=nr)
        if R - L + 1 < 3:
            return dict(end=1, steps=steps, reads=reads, ntrim=ntrim, nl=nl, nr=nr)
        steps += 1
        if steps + ntrim > cap or steps > maxsteps:
            return dict(end=(3 if ntrim == 0 else 6), steps=steps, reads=reads,
                        ntrim=ntrim, nl=nl, nr=nr)
        key = (ell, r, L, R)
        m = r - ell + 1
        if key in seen:
            if not repair:
                return dict(end=4, steps=steps, reads=reads, ntrim=0, nl=0, nr=0)
            if cyc is None:
                cyc = (steps - seen[key], m)
            ntrim += 1
            act = act_of(rule, S, ell, r)
            side = 1 if (act == LPLUS or act == LMINUS or act == CUTL) else 2
            if V.conv(rule) == 1:
                side = 3 - side
            if side == 1:
                L += 1
                nl += 1
            else:
                R -= 1
                nr += 1
            if L > R:
                return dict(end=1, steps=steps, reads=reads, ntrim=ntrim,
                            nl=nl, nr=nr)
            ell, r = L, R
            reads += abs(y - L) + abs(x - R)
            y, x = L, R
            k = max(1, (x - y + 1) >> 1)
            ty, tx = y, y + k - 1
            reads += abs(hy - ty) + abs(hx - tx)
            hy, hx = ty, tx
            continue
        seen[key] = steps
        act = act_of(rule, S, ell, r)
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
                            nl=nl, nr=nr)
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
                            nl=nl, nr=nr)
            ell, r = L, R
            reads += abs(y - L) + abs(x - R)
            y, x = L, R
            k = max(1, (x - y + 1) >> 1)
            ty, tx = y, y + k - 1
            reads += abs(hy - ty) + abs(hx - tx)
            hy, hx = ty, tx
            continue
        if ell + dy >= r + dx:
            return dict(end=0, steps=steps, reads=reads, ntrim=ntrim, nl=nl, nr=nr)
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
    nt = int(sys.argv[1]) if len(sys.argv) > 1 else 3
    rules = ['canonical', 'dense2', 'pattern', 'BAL>m/2', 'BAL>m/4', 'NBLK>3',
             'AMAX', 'longrun', 'pal', 'hashmod']
    bad = 0
    tot = 0
    cap = 40 * n * n + 100
    for rule in rules:
        for i, a, s, X in T.tapes(n, nt):
            S = P.St(X)
            ref = trace(rule, S, n, cap, repair=True)
            C = E.Cyc(X)
            got = C.run(V.rrid(rule), V.conv(rule), V.variant(rule), mode=4)
            tot += 1
            for name, a_, b_ in (('end', ref['end'], got['end']),
                                 ('steps', ref['steps'], got['steps']),
                                 ('reads', ref['reads'], got['reads']),
                                 ('ntrim', ref['ntrim'], got['ntrim']),
                                 ('nl', ref['nl'], got['nl']),
                                 ('nr', ref['nr'], got['nr'])):
                if a_ != b_:
                    bad += 1
                    print('DIFF %-9s tape %d %-6s python %s numba %s'
                          % (rule, i, name, a_, b_), flush=True)
        print('%-9s %d tapes checked' % (rule, nt), flush=True)
    print('cut_indep n=%d: %d/%d tape-runs agree on all 6 fields'
          % (n, tot - bad, tot))
    json.dump(dict(n=n, tapes=tot, bad=bad),
              open(os.path.join(HERE, 'cut_indep.json'), 'w'), indent=1)
    return 0 if bad == 0 else 1


if __name__ == '__main__':
    sys.exit(main())
