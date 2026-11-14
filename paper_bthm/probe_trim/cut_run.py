# -*- coding: utf-8 -*-
"""cut_run -- one cell of the CUT-repair measurement.

usage: cut_run.py <n> <rule> [ntapes]

For every tape of tr_lib.tapes(n) the same two runs as probe_cycle/cyc_run.py,
except that the repair is the companion paper's real trim:

  BTHM-0    no memory; stops at the first repeated quadruple (mode 0).  Its
            numbers must equal probe_cycle's stored rows exactly -- the
            no-memory walk is untouched by the repair.
  BTHM-CUT  the one-bit oracle; on a repeat the rule removes the BOUNDARY CELL
            OF THE TAPE on the side its own action table selects for that state
            (L <- L+1 or R <- R-1), re-opens the window on the narrower tape and
            continues (mode 4).  The tape only shrinks, so at most n repairs are
            possible; the step cap and the charge are unchanged.

reads: pcr_lib's INC accounting in both models (opening pass n + 2(n-1) + 2 +
n//2, 1 per cell taken in or given up, re-span after a cut, head re-read).  The
oracle bit is free; a repair pays the re-span a rule-initiated cut pays
(abs(y-L) + abs(x-R_) plus the head re-read), never an oracle read.

writes cut_<rule-safe>_<n>.json
"""
import json
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import cut_engine as E                                           # noqa: E402
import cut_verify as V                                           # noqa: E402
import tr_lib as T                                               # noqa: E402

HLOG = 24          # 16.7M slots: a cut run visits far more states than BTHM-0
                   # (many repairs per tape); the table-full guard is 3/4 of it.


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
        r4 = C.run(V.rrid(rule), V.conv(rule), V.variant(rule), mode=4)
        rows.append(dict(i=i, end0=r0['end'], steps0=r0['steps'],
                         reads0=r0['reads'], cyc_start=r0['cyc_start'],
                         cyc_len=r0['cyc_len'], width=r0['width'],
                         cuts0=r0['cuts'], endC=r4['end'], stepsC=r4['steps'],
                         readsC=r4['reads'], ntrim=r4['ntrim'],
                         cutsC=r4['cuts'], nl=r4['nl'], nr=r4['nr'],
                         nused=r4['nused']))
    out = dict(n=n, rule=rule, tapes=nt,
               repair='cut: on a repeat remove the boundary cell of the tape '
                      '(L+=1 or R-=1) on the side the rule action table '
                      'selects for that state',
               mode=4, cap=cap, hlog=hlog,
               rrid=V.rrid(rule), variant=V.variant(rule), conv=V.conv(rule),
               rows=rows, secs=time.time() - t0)
    f = 'cut_%s_%d.json' % (safe(rule), n)
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
    eC = [x['endC'] for x in r]
    print('%-10s n=%-5d meet0=%3d cyc0=%3d empty0=%3d cap0=%3d | cut: meet=%3d '
          'empty=%3d cap=%3d full=%3d trims med=%-6.0f max=%-5d [%.0fs]'
          % (rule, n, e0.count(0), e0.count(4), e0.count(1), e0.count(3),
             eC.count(0), eC.count(1), eC.count(3) + eC.count(6),
             eC.count(5),
             np.median([x['ntrim'] for x in r]) if r else -1,
             max([x['ntrim'] for x in r], default=-1), out['secs']), flush=True)


if __name__ == '__main__':
    main()
