# -*- coding: utf-8 -*-
"""par_run -- one cell of the BTHM-C_2 measurement.

usage: par_run.py <n> <rule> [ntapes]

For every tape of tr_lib.tapes(n) the same pair of runs as probe_trim/cut_run.py,
with the flag replaced by the parity colour:

  BTHM-0    no memory; stops at the first repeated quadruple (mode 0).  Its
            numbers must equal probe_trim's stored cut_<rule>_<n>.json rows
            bit for bit -- the no-memory walk is untouched by the new model.
  BTHM-C_2  the parity colour (mode 5): one bit per window state (ell, r, L, R),
            toggled on every visit -- black iff an odd number of visits so far.
            The rule's action tables are UNCHANGED (variant (a), colour-blind
            rules); the colour is used only by the repair, which fires when the
            window READ is black: visit 1 reads white, visit 2 black (trim),
            visit 3 white again.  The repair is exactly mode 4's: remove the
            boundary cell of the tape on the side the action table selects,
            re-open the window, continue; the colour table is never cleared.

reads: pcr_lib's INC accounting in both models; a repair pays the re-span a
rule-initiated cut pays, never an oracle read.

writes par_<rule-safe>_<n>.json
"""
import json
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import par_engine as E                                           # noqa: E402
import par_verify as V                                           # noqa: E402
import tr_lib as T                                               # noqa: E402

HLOG = 24          # same table as probe_trim/cut_run.py, so the two models
                   # face the same instrument


def safe(name):
    return (name.replace('>', 'g').replace('=', 'e').replace('/', 'd')
            .replace('.', 'p'))


def one(rule, n, nt, hlog=HLOG):
    t0 = time.time()
    cap = 40 * n * n + 100
    rows = []
    for i, a, s, X in T.tapes(n, nt):
        C = E.Cyc(X, hlog=hlog)
        r0 = C.run(V.rrid(rule), V.conv(rule), V.variant(rule), mode=0)
        r5 = C.run(V.rrid(rule), V.conv(rule), V.variant(rule), mode=5)
        rows.append(dict(i=i, end0=r0['end'], steps0=r0['steps'],
                         reads0=r0['reads'], cyc_start=r0['cyc_start'],
                         cyc_len=r0['cyc_len'], width=r0['width'],
                         cuts0=r0['cuts'], endP=r5['end'], stepsP=r5['steps'],
                         readsP=r5['reads'], ntrim=r5['ntrim'],
                         cutsP=r5['cuts'], nl=r5['nl'], nr=r5['nr'],
                         nused=r5['nused']))
    out = dict(n=n, rule=rule, tapes=nt,
               repair='parity: on a BLACK window (odd number of visits so far) '
                      'remove the boundary cell of the tape (L+=1 or R-=1) on '
                      'the side the rule action table selects for that state; '
                      'the colour bit toggles on every visit and the table is '
                      'never cleared',
               mode=5, cap=cap, hlog=hlog,
               rrid=V.rrid(rule), variant=V.variant(rule), conv=V.conv(rule),
               rows=rows, secs=time.time() - t0)
    f = 'par_%s_%d.json' % (safe(rule), n)
    with open(os.path.join(HERE, f), 'w') as fh:
        json.dump(out, fh)
    return f, out


def main():
    n = int(sys.argv[1])
    rule = sys.argv[2]
    nt = int(sys.argv[3]) if len(sys.argv) > 3 else T.NT
    f, out = one(rule, n, nt)
    r = out['rows']
    e0 = [x['end0'] for x in r]
    eP = [x['endP'] for x in r]
    ntr = [x['ntrim'] for x in r]
    print('%-10s n=%-5d meet0=%3d cyc0=%3d empty0=%3d cap0=%3d | par: meet=%3d '
          'empty=%3d cap=%3d full=%3d trims med=%-6.0f max=%-5d [%.0fs]'
          % (rule, n, e0.count(0), e0.count(4), e0.count(1), e0.count(3),
             eP.count(0), eP.count(1), eP.count(3) + eP.count(6),
             eP.count(5),
             np.median([x['ntrim'] for x in r]) if r else -1,
             max([x['ntrim'] for x in r], default=-1), out['secs']), flush=True)


if __name__ == '__main__':
    main()
