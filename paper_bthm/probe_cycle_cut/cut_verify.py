# -*- coding: utf-8 -*-
"""cyc_verify -- is the instrumented walk the record's walk?

Two independent checks, both before any cycle number is quoted.

A  gw_cells.json (probe_g_witness), n = 512, nine witness tapes, the eleven
   rules whose tr_lib machine is charged with pcr_lib's INC accounting
   (canonical, dense2, peak, trip, oddrun, pattern, LAS, MAXP, KCHG, BAL>m/2,
   BAL>m/4): my reads, delta and end must equal the recorded cell exactly.
   monotone must give 2n-1 reads.  This is the authoritative cell set for the
   model, so a match here is a match on reads, not only on the walk.

B  tr_lib itself, n = 512, six random tapes, all 21 catalogue rules, LEGACY
   mode (the same (ell, r) stamp pcr_lib carries): delta and end must agree.
   Steps are not compared (tr_lib counts the endgame test differently), reads
   are not compared (each machine has its own charge).
"""
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import cut_engine as E                                           # noqa: E402

PAPER = E.PAPER
sys.path.insert(0, os.path.join(PAPER, 'probe_g_witness'))
sys.path.insert(0, os.path.join(PAPER, 'probe_triple_nat'))
import tr_lib as T                                               # noqa: E402
import pcr_lib as P                                              # noqa: E402
import gw                                                       # noqa: E402

INC = ('canonical', 'dense2', 'peak', 'trip', 'oddrun', 'pattern', 'LAS',
       'MAXP', 'KCHG', 'BAL>m/2', 'BAL>m/4')
RUL = [r for r, _, _, _, _ in T.SPEC]


def rid(rule):
    return P.RID[rule]


def variant(rule):
    return {'BAL>m/4': 1, 'NBLK>3': 2, 'hashmod': 3}.get(rule, 0)


def rrid(rule):
    if rule == 'BAL>m/4':
        return P.RID['BAL>m/2']
    if rule == 'NBLK>3':
        return P.RID['NBLK>2']
    if rule == 'hashmod':
        return P.R_HASHMOD
    if rule == 'monotone':
        return -1
    return P.RID[rule]


def conv(rule):
    return 1 if rule == 'canonical' else 0


def endname(e):
    if isinstance(e, str):
        return e
    return {0: 'meet', 1: 'empty', 2: 'repeat', 3: 'cap', None: 'none'}.get(e, str(e))


def main():
    n = 512
    bad = []
    # ---------------- A: gw_cells, INC rules, nine witnesses
    cells = json.load(open(os.path.join(PAPER, 'probe_g_witness', 'gw_cells.json')))[str(n)]
    na = nb = 0
    for spec in gw.SPECS:
        a, s, X = gw._cols(n, spec)
        C = E.Cyc(X)
        for rule in INC:
            r = C.run(rrid(rule), conv(rule), variant(rule), mode=2)
            key = '%s|%s' % (rule, spec)
            want = cells[key]
            got_d = r['delta']
            wd = want['delta']
            if wd == -1.0:
                wd = None
            okd = (wd is None and got_d is None) or (
                wd is not None and got_d is not None
                and abs(wd - got_d) < 1e-12)
            oke = endname(want['end']) == endname(r['end'])
            okr = int(want['reads']) == r['reads']
            na += 1
            if not (okd and oke and okr):
                bad.append(('A', key, want['reads'], r['reads'], want['delta'], got_d,
                            endname(want['end']), endname(r['end'])))
        # monotone: 2n-1 by theorem
        mo = E.mono_run(a)
        nb += 1
        if mo['reads'] != 2 * n - 1:
            bad.append(('A', 'monotone|%s' % spec, 2 * n - 1, mo['reads'], 0, 0, 0, 0))
    print('A  gw_cells n=%d: %d/%d INC cells exact, %d monotone cells' % (n, na - len([b for b in bad if b[0] == 'A']), na, nb), flush=True)

    # ---------------- B: tr_lib, all rules, legacy mode
    nc = 0
    for i, a, s, X in T.tapes(n, 6):
        C = E.Cyc(X)
        for rule in RUL:
            if rule in ('monotone',):
                continue
            rt = T.run(rule, a, s)
            rm = C.run(rrid(rule), conv(rule), variant(rule), mode=2)
            dt = rt.get('delta')
            dm = rm['delta']
            if dt == -1.0:
                dt = None
            okd = (dt is None and dm is None) or (dt is not None and dm is not None
                                                  and abs(float(dt) - dm) < 1e-12)
            et = endname(rt.get('end'))
            em = endname(rm['end'])
            if rule in ('pal', 'longrun'):      # those machines: end=1 == no delta
                oke = (em == 'meet') if et == 'meet' else (em != 'meet')
            else:
                oke = et == em
            nc += 1
            if not (okd and oke):
                bad.append(('B', rule, i, rt.get('delta'), dm, et, em,
                            rt.get('steps'), rm['steps']))
    print('B  tr_lib n=%d: %d/%d (delta, end) agree' % (n, nc - len([b for b in bad if b[0] == 'B']), nc), flush=True)

    for b in bad:
        print('MISMATCH', b, flush=True)
    out = dict(n=n, A_cells=na, B_cells=nc, bad=bad)
    json.dump(out, open(os.path.join(HERE, 'cut_verify.json'), 'w'), indent=1)
    print('bad:', len(bad))
    return 0 if not bad else 1


if __name__ == '__main__':
    sys.exit(main())
