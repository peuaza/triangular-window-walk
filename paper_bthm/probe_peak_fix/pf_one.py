# -*- coding: utf-8 -*-
"""pf_one -- one tape, the trace, the EN_LINEAR machine and its first divergence.

Usage: python pf_one.py N I [RULE]
"""
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
PAPER = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(PAPER, 'probe_en_linear'))

import el_cell as C
import el_lib as E

rule = sys.argv[3] if len(sys.argv) > 3 else 'peak'
n = int(sys.argv[1])
i = int(sys.argv[2])

X = E.tape(n, i)
M = E.Mach(X, n)
t = M.trace(rule)
print('tape n=%d i=%d  rule=%s' % (n, i, rule))
print('  truth: end=%d delta=%.12f steps=%d cuts=%d nep=%d termM0=%d'
      % (t['end'], t['delta'], t['steps'], t['cuts'], t['nep'], t['epM'][-1]))

# first mispredicted epoch of the truth
bad = []
for k in range(t['nep']):
    if t['epb2'][k] < 0 or t['eplen'][k] == 0:
        continue
    pred = C.pred_of(rule, int(t['epb2'][k]), int(t['epb1'][k]))
    if int(t['epout'][k]) != pred:
        bad.append((k, int(t['epM'][k]), pred, int(t['epout'][k]),
                    int(t['eplen'][k])))
print('  mispredicted truth epochs: %d' % len(bad))
for b in bad[:6]:
    print('    epoch %d  M0=%d  pred=%d  out=%d  len=%d' % b)

for W in (32, 64, 128, 256, 512):
    d, en, rd, epc, cu, eg, nsc = C.jm_walk(rule, M.S, n, W, 'sliding')
    ok = (int(en) == int(t['end'])) and (
        (d is None and t['delta'] < 0) or
        (d is not None and t['delta'] >= 0 and abs(d - t['delta']) < 1e-12))
    print('  W2=%3d  exact=%d  end=%d  delta=%s  reads=%d (e/n %.3f)  epochs=%d eg=%d'
          % (W, ok, en, 'None' if d is None else '%.9f' % d, rd, rd / n, epc, eg))

# epoch where the first epoch ends: M0 of epoch 1 and its length
if t['nep'] > 1:
    print('  epoch0: M0=%d len=%d out=%d ; epoch1: M0=%d len=%d out=%d'
          % (t['epM'][0], t['eplen'][0], t['epout'][0],
             t['epM'][1], t['eplen'][1], t['epout'][1]))
print('  cumulative steps at epoch 443:', int(np.sum(t['eplen'][:443])))
