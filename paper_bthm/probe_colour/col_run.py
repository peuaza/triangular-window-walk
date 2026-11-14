# -*- coding: utf-8 -*-
"""col_run -- one cell of the BTHM-C_k measurement.

usage: col_run.py <n> <rule> <kmax> [ntapes]

rule in {pattern, runs, peak} (the three variant-0, conv-0 rules of CUT.md);
kmax = 4 is the colour model (BTHM-C_4); kmax = 1 is the control, which must
reproduce probe_trim/cut_<rule>_<n>.json tape for tape (same oracle use: one
trim per repeat; colour pinned to 0).

For every tape of tr_lib.tapes(n):
  mode 0  BTHM-0 (col_engine's copy of cut_engine mode 0): the walk stops at
          the first repeated quadruple.  Its seven fields (end0, steps0,
          reads0, cyc_start, cyc_len, width, cuts0) are compared tape by tape
          with the stored probe_trim/cut_<rule>_<n>.json rows -- the same tapes
          and the same no-memory walk as the colour-blind BTHM-C run, so the
          comparison in COLOUR.md is on identical tapes.
  mode 5  the colour model at kmax: visits 1..kmax of a state are played (the
          rule reads the colour), the (kmax+1)-th fires the cut repair.
          The oracle bit (the colour) is free and takes no step; a repair
          pays the re-span a rule-initiated cut pays, abs(y-L)+abs(x-R_)
          plus the head re-read -- unchanged from cut_engine mode 4.

writes col_<rule-safe>_<n>.json (kmax=4) or col_k<n>_<rule-safe>_<n>.json.
"""
import json
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import col_engine as E                                          # noqa: E402
import tr_lib as T                                              # noqa: E402
import pcr_lib as P                                             # noqa: E402

TRIM = os.path.join(E.PAPER, 'probe_trim')
K1 = os.path.join(E.PAPER, 'probe_k300')
RULES = ('pattern', 'runs', 'peak')
B0F = ('end0', 'steps0', 'reads0', 'cyc_start', 'cyc_len', 'width', 'cuts0')


def safe(name):
    return (name.replace('>', 'g').replace('=', 'e').replace('/', 'd')
            .replace('.', 'p'))


def stored_rows(rule, n):
    for d in (TRIM, K1):
        f = os.path.join(d, 'cut_%s_%d.json' % (safe(rule), n))
        if os.path.exists(f):
            return {int(r['i']): r for r in json.load(open(f))['rows']}, f
    return None, None


def one(rule, n, kmax, nt, hlog=E.HLOG):
    assert rule in RULES, rule
    rid = P.RID[rule]
    conv = 0
    cap = 40 * n * n + 100
    old, oldf = stored_rows(rule, n)
    rows = []
    bad0 = []
    t0 = time.time()
    for i, a, s, X in T.tapes(n, nt):
        C = E.Col(X, hlog=hlog)
        r0 = C.run(rid, conv, mode=0)
        r5 = C.run(rid, conv, mode=5, kmax=kmax)
        row = dict(i=i, end0=r0['end'], steps0=r0['steps'], reads0=r0['reads'],
                   cyc_start=r0['cyc_start'], cyc_len=r0['cyc_len'],
                   width=r0['width'], cuts0=r0['cuts'], endC=r5['end'],
                   stepsC=r5['steps'], readsC=r5['reads'], ntrim=r5['ntrim'],
                   cutsC=r5['cuts'], nl=r5['nl'], nr=r5['nr'],
                   nused=r5['nused'], kc0=r5['kc0'], kc1=r5['kc1'],
                   kc2=r5['kc2'], kc3=r5['kc3'], repok=r5['repok'],
                   fstep=r5['fstep'])
        if old is not None and i in old:
            o = old[i]
            d = [f for f in B0F if int(o[f]) != int(row[f])]
            if d:
                bad0.append((i, d, [o[f] for f in d], [row[f] for f in d]))
                for f in d:
                    row[f + '_stored'] = o[f]
        rows.append(row)
    e0 = [x['end0'] for x in rows]
    eC = [x['endC'] for x in rows]
    out = dict(n=n, rule=rule, kmax=kmax, tapes=nt, hlog=hlog, cap=cap,
               rrid=rid, conv=conv,
               model=('BTHM-C_4: oracle = visit count mod 4; colour 0 first '
                      'visit, 1 second, 2 third, 3 fourth; repair at the 5th'
                      if kmax == 4 else
                      'BTHM-C control: kmax=1, colour pinned to 0, repair at '
                      'the 2nd visit'),
               repair='cut_engine mode 4: on the (kmax+1)-th visit remove the '
                      'boundary cell of the tape (L+=1 or R-=1) on the side '
                      'the colour-0 (base) action table selects',
               rows=rows, bad0=bad0,
               b0_stored=oldf, secs=time.time() - t0)
    f = ('col_%s_%d.json' % (safe(rule), n) if kmax == 4 else
         'col_k%d_%s_%d.json' % (kmax, safe(rule), n))
    with open(os.path.join(HERE, f), 'w') as fh:
        json.dump(out, fh)
    print('%-8s n=%-5d kmax=%d meet0=%3d cyc0=%3d cap0=%3d | col: meet=%3d '
          'empty=%3d cap=%3d full=%3d rep=%3d trims med=%-6.0f max=%-5d '
          'steps=%d/%d/%d/%d repok_bad=%d b0bad=%d [%.0fs]'
          % (rule, n, kmax, e0.count(0), e0.count(4), e0.count(3),
             eC.count(0), eC.count(1), eC.count(3) + eC.count(6),
             eC.count(5), sum(x['ntrim'] for x in rows),
             np.median([x['ntrim'] for x in rows]) if rows else -1,
             max([x['ntrim'] for x in rows], default=-1),
             sum(x['kc0'] for x in rows), sum(x['kc1'] for x in rows),
             sum(x['kc2'] for x in rows), sum(x['kc3'] for x in rows),
             sum(1 for x in rows if x['repok'] != x['ntrim']),
             len(bad0), out['secs']), flush=True)
    return f, out


def main():
    n = int(sys.argv[1])
    rule = sys.argv[2]
    kmax = int(sys.argv[3]) if len(sys.argv) > 3 else 4
    nt = int(sys.argv[4]) if len(sys.argv) > 4 else T.NT
    one(rule, n, kmax, nt)


if __name__ == '__main__':
    main()
