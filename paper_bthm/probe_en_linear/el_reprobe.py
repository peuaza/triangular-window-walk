# -*- coding: utf-8 -*-
"""el_reprobe -- peak on a different segment of the tape family.

The headline claim is that `peak` is exact on 700/700 meeting tapes at W2 = 128.
That is one sample.  This re-runs the same machine on tapes i = OFF, OFF+1, ...
(a disjoint segment of the same family) at n = 1024 and 2048, and prints the
exactness and e/n at each W2.

Usage: python el_reprobe.py RULE N OFF NMEET
"""
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import el_cell as C
import el_lib as E

rule = sys.argv[1]
n = int(sys.argv[2])
off = int(sys.argv[3])
nmeet = int(sys.argv[4]) if len(sys.argv) > 4 else 200
mode = 'sliding' if E.EVAL[rule] == 'sliding' else 'scan'
i = off
meet = []
t0 = __import__('time').time()
while len(meet) < nmeet and i < off + 200000:
    X = E.tape(n, i)
    M = E.Mach(X, n)
    t = M.trace(rule)
    if t['end'] == 0:
        rec = dict(i=i, mach={})
        for W in (32, 64, 128, 256):
            d, en, rd, epc, cu, egs, nsc = C.jm_walk(rule, M.S, n, W, mode)
            ok = (int(en) == int(t['end'])) and (
                (d is None and t['delta'] < 0) or
                (d is not None and t['delta'] >= 0 and abs(d - t['delta']) < 1e-12))
            rec['mach'][str(W)] = (int(ok), rd)
        meet.append(rec)
    i += 1
print('%s n=%d offset=%d meetings=%d draws=%d secs=%.0f' % (
    rule, n, off, len(meet), i - off, __import__('time').time() - t0))
for W in (32, 64, 128, 256):
    ex = sum(r['mach'][str(W)][0] for r in meet)
    e = np.mean([r['mach'][str(W)][1] for r in meet]) / n
    print('   W2=%3d  exact %3d/%3d   e/n %.2f' % (W, ex, len(meet), e))
json.dump([{k: v for k, v in r.items()} for r in meet],
          open(os.path.join(HERE, 'el_reprobe_%s_%d_%d.json'
                            % (C.safe(rule), n, off)), 'w'))
