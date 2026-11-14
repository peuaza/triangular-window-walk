# -*- coding: utf-8 -*-
"""k_run -- one cell of the BTHM-C_k sweep.

usage: k_run.py <n> <rule> <k> [ntapes]

rule in {pattern, runs, peak} (the three variant-0, conv-0 rules of CUT.md);
k is the model parameter: colour = (visit count) mod k, repair on the
(k+1)-th visit (count positive and count mod k == 0).  k = 1 is the
colour-blind model (BTHM-C) and must reproduce probe_trim/cut_cells.json;
k = 4 must reproduce the already-stored probe_colour engine runs.

For every tape of tr_lib.tapes(n):
  mode 0  BTHM-0 (k_engine's copy of cut_engine mode 0): the walk stops at
          the first repeated quadruple.  Its seven fields (end0, steps0,
          reads0, cyc_start, cyc_len, width, cuts0) are compared tape by tape
          with the stored probe_trim/cut_<rule>_<n>.json rows.
  mode 5  the colour model at this k.

Checks written into the cell file:
  b0bad  -- mode-0 fields vs the stored cut rows (same tapes, same no-memory
            walk), must be 0 for every cell;
  k1bad  -- (only k=1) the colour run vs the stored cut rows on
            (endC, stepsC, readsC, ntrim, cutsC, nl, nr, nused): the k=1 gate;
  k4bad  -- (only k=4) the colour run vs the stored probe_colour/col_*.json
            rows on the same eight fields plus kc0..kc3: continuity with the
            BTHM-C_4 measurement.

writes k<k>_<rule-safe>_<n>.json
"""
import json
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import k_engine as E                                            # noqa: E402
import tr_lib as T                                              # noqa: E402
import pcr_lib as P                                             # noqa: E402

TRIM = os.path.join(E.PAPER, 'probe_trim')
COLOUR4 = os.path.join(E.PAPER, 'probe_colour')
RULES = ('pattern', 'runs', 'peak')
B0F = ('end0', 'steps0', 'reads0', 'cyc_start', 'cyc_len', 'width', 'cuts0')
C1F = ('endC', 'stepsC', 'readsC', 'ntrim', 'cutsC', 'nl', 'nr', 'nused')
C4F = C1F + ('kc0', 'kc1', 'kc2', 'kc3')


def safe(name):
    return (name.replace('>', 'g').replace('=', 'e').replace('/', 'd')
            .replace('.', 'p'))


def load_rows(path):
    if not os.path.exists(path):
        return None
    return {int(r['i']): r for r in json.load(open(path))['rows']}


def one(rule, n, k, nt, hlog=E.HLOG):
    assert rule in RULES, rule
    assert k >= 1
    rid = P.RID[rule]
    conv = 0
    cap = 40 * n * n + 100
    old = load_rows(os.path.join(TRIM, 'cut_%s_%d.json' % (safe(rule), n)))
    old4 = (load_rows(os.path.join(COLOUR4, 'col_%s_%d.json'
                                   % (safe(rule), n))) if k == 4 else None)
    rows = []
    b0bad = []
    c1bad = []
    c4bad = []
    t0 = time.time()
    for i, a, s, X in T.tapes(n, nt):
        C = E.Col(X, hlog=hlog)
        r0 = C.run(rid, conv, mode=0)
        r5 = C.run(rid, conv, mode=5, k=k)
        row = dict(i=i, end0=r0['end'], steps0=r0['steps'], reads0=r0['reads'],
                   cyc_start=r0['cyc_start'], cyc_len=r0['cyc_len'],
                   width=r0['width'], cuts0=r0['cuts'], endC=r5['end'],
                   stepsC=r5['steps'], readsC=r5['reads'], ntrim=r5['ntrim'],
                   cutsC=r5['cuts'], nl=r5['nl'], nr=r5['nr'],
                   nused=r5['nused'], kc0=r5['kc0'], kc1=r5['kc1'],
                   kc2=r5['kc2'], kc3=r5['kc3'], repok=r5['repok'],
                   fstep=r5['fstep'], rmod=r5['rmod'])
        if old is not None and i in old:
            o = old[i]
            d = [f for f in B0F if int(o[f]) != int(row[f])]
            if d:
                b0bad.append((i, d, [o[f] for f in d], [row[f] for f in d]))
        if k == 1 and old is not None and i in old:
            o = old[i]
            d = [f for f in C1F if int(o[f]) != int(row[f])]
            if d:
                c1bad.append((i, d, [o[f] for f in d], [row[f] for f in d]))
        if k == 4 and old4 is not None and i in old4:
            o = old4[i]
            d = [f for f in C4F if int(o[f]) != int(row[f])]
            if d:
                c4bad.append((i, d, [o[f] for f in d], [row[f] for f in d]))
        rows.append(row)
    e0 = [x['end0'] for x in rows]
    eC = [x['endC'] for x in rows]
    out = dict(
        n=n, rule=rule, k=k, tapes=nt, hlog=hlog, cap=cap, rrid=rid, conv=conv,
        model='BTHM-C_k: oracle = EXACT visit count of (ell, r, L, R); the '
              'rule reads count mod k (k colour-variants); the repair fires '
              'when the count is positive and count mod k == 0, i.e. on the '
              '(k+1)-th visit, so k visits are played before any trim.  '
              'k=%d: colours 0..%d, repair on visit %d.' % (k, k - 1, k + 1),
        design='fixed across k: mirror the move when c is odd '
               '(LPLUS<->RMINUS, LMINUS<->RPLUS; trims and halt fixed); swap '
               'the cut side when c is 2 or 3 (CUTL<->CUTR).  Boundary cut '
               'side: the side the (mapped) move would leave for c in {0,1}, '
               'the opposite side for c in {2,3}; double cut kept for every '
               'colour.  Repair reads the side off the colour-0 (base) table.',
        repair='cut_engine mode 4: on the (k+1)-th visit remove the boundary '
               'cell of the tape (L+=1 or R-=1) on the side the colour-0 '
               '(base) action table selects for that state',
        rows=rows, b0bad=b0bad, b0_stored=('probe_trim/cut_%s_%d.json'
                                           % (safe(rule), n)),
        c1bad=c1bad, c4bad=c4bad,
        c1_stored=('probe_trim/cut_%s_%d.json' % (safe(rule), n)),
        c4_stored=('probe_colour/col_%s_%d.json' % (safe(rule), n)),
        secs=time.time() - t0)
    f = 'k%d_%s_%d.json' % (k, safe(rule), n)
    with open(os.path.join(HERE, f), 'w') as fh:
        json.dump(out, fh)
    print('%-8s n=%-5d k=%d meet0=%3d cyc0=%3d cap0=%3d | k-run: meet=%3d '
          'empty=%3d cap=%3d full=%3d rep=%3d trims med=%-6.0f max=%-5d '
          'repinv_bad=%d rt_bad=%d b0bad=%d k1bad=%s k4bad=%s [%.0fs]'
          % (rule, n, k, e0.count(0), e0.count(4), e0.count(3),
             eC.count(0), eC.count(1), eC.count(3) + eC.count(6),
             eC.count(5), sum(x['ntrim'] for x in rows),
             np.median([x['ntrim'] for x in rows]) if rows else -1,
             max([x['ntrim'] for x in rows], default=-1),
             sum(1 for x in rows if x['repok'] != x['ntrim']),
             sum(1 for x in rows if x['rmod'] != x['ntrim']),
             len(b0bad), len(c1bad) if k == 1 else 'n/a',
             len(c4bad) if k == 4 else 'n/a', out['secs']), flush=True)
    return f, out


def main():
    n = int(sys.argv[1])
    rule = sys.argv[2]
    k = int(sys.argv[3])
    nt = int(sys.argv[4]) if len(sys.argv) > 4 else T.NT
    one(rule, n, k, nt)


if __name__ == '__main__':
    main()
