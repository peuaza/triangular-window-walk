# -*- coding: utf-8 -*-
"""el_fail -- inspect the meeting tapes the machine gets wrong in a reprobe."""
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
mode = 'sliding' if E.EVAL[rule] == 'sliding' else 'scan'
recs = json.load(open(os.path.join(
    HERE, 'el_reprobe_%s_%d_%d.json' % (C.safe(rule), n, off))))
bad = [r for r in recs if not r['mach']['256'][0]]
print('%s n=%d offset=%d  meetings %d, wrong at W2=256: %d'
      % (rule, n, off, len(recs), len(bad)))
for r in bad:
    X = E.tape(n, r['i'])
    M = E.Mach(X, n)
    t = M.trace(rule)
    meet = [b for b in []]
    big = [b for b in []]
    mm = ''
    worst = [b for b in t and []]
    print('  tape i=%d  truth delta %.6f end %d  nep %d  termM0 %d'
          % (r['i'], t['delta'], t['end'], t['nep'], t['epM'][-1]))
    for W in (32, 64, 128, 256, 512):
        d, en, rd, epc, cu, egs, nsc = C.jm_walk(rule, M.S, n, W, mode)
        ok = (int(en) == int(t['end'])) and (
            (d is None and t['delta'] < 0) or
            (d is not None and t['delta'] >= 0 and abs(d - t['delta']) < 1e-12))
        print('     W2=%3d  mach delta %s end %d  exact %d'
              % (W, 'None' if d is None else '%.6f' % d, en, ok))
    bads = [(int(b[0]), int(b[1]), int(b[2]), int(b[5])) for b in []]
    bl = []
    for k in range(t['nep']):
        if t['epb2'][k] < 0 or t['eplen'][k] == 0:
            continue
        pred = C.pred_of(rule, int(t['epb2'][k]), int(t['epb1'][k]))
        if int(t['epout'][k]) != pred:
            bl.append((int(t['epM'][k]), pred, int(t['epout'][k]),
                       int(t['eplen'][k])))
    print('     mispredicted epochs (M0, pred, out, len):', bl[:10],
          '... total', len(bl))
