# -*- coding: utf-8 -*-
"""cyc_run -- one cell: 300 tapes, BTHM-0 and BTHM-C on the same tape.

usage: cyc_run.py <n> <rule> [ntapes] [trim]
writes cyc_<rule-safe>_<n>.json  (and ..._trim2.json when trim=2)

Per tape:
  BTHM-0  no memory.  end0=0 meet, 1 region empty, 3 step cap, 4 quadruple
          cycle (the walk closed a period; by determinism it then repeats for
          ever, so under BTHM-0 it dies at the cap), 5 oracle table full.
          steps0 = steps to the terminal / to the first repeat.
  BTHM-C  the one-bit oracle; on a repeat trim one cell off the window and
          continue.  endC, stepsC, readsC, ntrim.
reads: pcr_lib's INC accounting in both models (opening pass
       n + 2(n-1) + 2 + n//2, 1 per cell taken in or given up, re-span after a
       cut, head re-read).  The oracle bit is free; a trim pays the cell it
       gives up plus the head re-read.
"""
import json
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import cyc_engine as E                                           # noqa: E402
import cyc_verify as V                                           # noqa: E402
import tr_lib as T                                               # noqa: E402


def safe(name):
    return (name.replace('>', 'g').replace('=', 'e').replace('/', 'd')
            .replace('.', 'p'))


def one(rule, n, nt, trim):
    t0 = time.time()
    rows = []
    for i, a, s, X in T.tapes(n, nt):
        C = E.Cyc(X)
        r0 = C.run(V.rrid(rule), V.conv(rule), V.variant(rule), mode=0)
        rC = C.run(V.rrid(rule), V.conv(rule), V.variant(rule), mode=1,
                   trim=trim)
        rows.append(dict(i=i, end0=r0['end'], steps0=r0['steps'],
                         reads0=r0['reads'], cyc_start=r0['cyc_start'],
                         cyc_len=r0['cyc_len'], width=r0['width'],
                         cuts0=r0['cuts'], endC=rC['end'], stepsC=rC['steps'],
                         readsC=rC['reads'], ntrim=rC['ntrim'],
                         cutsC=rC['cuts']))
    out = dict(n=n, rule=rule, tapes=nt, trim=trim,
               rrid=V.rrid(rule), variant=V.variant(rule), conv=V.conv(rule),
               rows=rows, secs=time.time() - t0)
    f = 'cyc_%s_%d%s.json' % (safe(rule), n, '' if trim == 1 else '_trim%d' % trim)
    with open(os.path.join(HERE, f), 'w') as fh:
        json.dump(out, fh)
    return f, out


def main():
    n = int(sys.argv[1])
    rule = sys.argv[2]
    nt = int(sys.argv[3]) if len(sys.argv) > 3 else T.NT
    trim = int(sys.argv[4]) if len(sys.argv) > 4 else 1
    f, out = one(rule, n, nt, trim)
    r = out['rows']
    e0 = [x['end0'] for x in r]
    print('%-10s n=%-5d trim=%d meet0=%3d cyc=%3d empty=%3d cap=%3d meetC=%3d '
          'med_cyc=%6.0f med_w=%5.0f [%.0fs]'
          % (rule, n, trim, e0.count(0), e0.count(4), e0.count(1), e0.count(3),
             [x['endC'] for x in r].count(0),
             np.median([x['cyc_len'] for x in r if x['end0'] == 4] or [0]),
             np.median([x['width'] for x in r if x['end0'] == 4] or [0]),
             out['secs']), flush=True)


if __name__ == '__main__':
    main()
