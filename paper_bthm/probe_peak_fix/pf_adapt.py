# -*- coding: utf-8 -*-
"""pf_adapt -- the adaptive machine on one segment of the peak tape family.

Usage: python pf_adapt.py N OFF NMEET OUT.jsonl

For every meeting tape: the trace, then `pf_var.jm_adapt` at W2 = 128 and 256.
The adaptive machine runs epochs literally until a literal epoch's cut agrees
with pred_of(b2,b1), then hands over to the jump machine.  The first epoch is
charged n (the "at most n, once" convention); later literal epochs are charged
their measured length.  Exactness is against the trace.
"""
import json
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import pf_lib as F
import pf_var as V
import el_lib as E

rule = 'peak'
n = int(sys.argv[1])
off = int(sys.argv[2])
nmeet = int(sys.argv[3])
outp = sys.argv[4]

i = off
recs = []
t0 = time.time()
while len(recs) < nmeet and i < off + 400000:
    X = E.tape(n, i)
    M = E.Mach(X, n)
    t = M.trace(rule)
    i += 1
    if t['end'] != 0:
        continue
    rec = dict(i=i - 1, n=n, end=int(t['end']), delta=float(t['delta']),
               steps=int(t['steps']), nep=int(t['nep']))
    for W in (128, 256):
        m = V.jm_adapt(rule, M.S, n, W, first_charge=n)
        rec['a%d' % W] = dict(exact=int(F.exact_of(m, t, n)),
                              reads=int(m['reads']), end=int(m['end']),
                              delta=(-1.0 if m['delta'] is None else float(m['delta'])),
                              lit_steps=int(m['lit_steps']),
                              lit_epochs=int(m['lit_epochs']))
    recs.append(rec)

with open(outp, 'w') as fh:
    for r in recs:
        fh.write(json.dumps(r) + '\n')
print('n=%d off=%d meetings=%d draws=%d secs=%.0f -> %s'
      % (n, off, len(recs), i - off, time.time() - t0, outp))
for W in (128, 256):
    ex = sum(r['a%d' % W]['exact'] for r in recs)
    e = np.mean([r['a%d' % W]['reads'] for r in recs]) / n
    lits = np.array([r['a%d' % W]['lit_steps'] for r in recs]) / n
    print('  adapt W2=%3d exact %3d/%3d  e/n %.3f  lit/n mean %.3f p90 %.3f max %.2f'
          % (W, ex, len(recs), e, lits.mean(), np.percentile(lits, 90), lits.max()))
