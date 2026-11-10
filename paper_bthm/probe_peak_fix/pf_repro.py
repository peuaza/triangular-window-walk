# -*- coding: utf-8 -*-
"""pf_repro -- one segment of the peak tape family, all machines, one JSONL row.

Usage: python pf_repro.py N OFF NMEET OUT.jsonl

For every meeting tape drawn from i = OFF, OFF+1, ...:
  * the trace (el_lib.trace_walk, gated against pcr_lib.Runner by el_gate.py);
  * the EN_LINEAR machine at W2 = 128 and 256 (el_cell.jm_walk semantics, with
    an epoch audit);
  * the fixed machine (literal first epoch, charged n, then the same jump
    machine) at W2 = 128 and 256.

The row carries the truth, each machine's (exact, reads, end, delta), the first
divergence (epoch index, M0, predicted out, true out), and the first 48 epochs
of the trace (M0, length, outcome, b2, b1) so the top-of-tape structure can be
read off the file rather than re-traced.
"""
import json
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import pf_lib as F
import el_cell as C
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
               steps=int(t['steps']), cuts=int(t['cuts']), nep=int(t['nep']),
               termM0=int(t['epM'][-1]))
    # truth epochs at the top
    top = []
    for k in range(min(48, t['nep'])):
        top.append([int(t['epM'][k]), int(t['eplen'][k]), int(t['epout'][k]),
                    int(t['epb2'][k]), int(t['epb1'][k])])
    rec['top'] = top
    bad = []
    for k in range(t['nep']):
        if t['epb2'][k] < 0 or t['eplen'][k] == 0:
            continue
        pred = C.pred_of(rule, int(t['epb2'][k]), int(t['epb1'][k]))
        if int(t['epout'][k]) != pred:
            bad.append([k, int(t['epM'][k]), pred, int(t['epout'][k]),
                        int(t['eplen'][k])])
    rec['nbad'] = len(bad)
    rec['bad'] = bad[:64]
    rec['steps_at_first_bad'] = int(np.sum(t['eplen'][:bad[0][0]])) if bad else -1

    for W in (128, 256):
        b = F.jm_base_audit(rule, M.S, n, W)
        fd = F.first_divergence(b, t)
        rec['base%d' % W] = dict(exact=int(F.exact_of(b, t, n)),
                                 reads=int(b['reads']), end=int(b['end']),
                                 delta=(-1.0 if b['delta'] is None else float(b['delta'])),
                                 epochs=int(b['epochs']), eg=int(b['eg']),
                                 fd=(-1 if fd is None else int(fd[0])),
                                 fdM0=(-1 if fd is None else int(fd[1])),
                                 fdpred=(-1 if fd is None else int(fd[2])),
                                 fdout=(-1 if fd is None else int(fd[3])))
        f = F.jm_fixed(rule, M.S, n, W)
        rec['fix%d' % W] = dict(exact=int(F.exact_of(f, t, n)),
                                reads=int(f['reads']), end=int(f['end']),
                                delta=(-1.0 if f['delta'] is None else float(f['delta'])),
                                epochs=int(f['epochs']), eg=int(f['eg']),
                                first_len=int(f['first_len']))
    recs.append(rec)

with open(outp, 'w') as fh:
    for r in recs:
        fh.write(json.dumps(r) + '\n')
print('n=%d off=%d meetings=%d draws=%d secs=%.0f -> %s'
      % (n, off, len(recs), i - off, time.time() - t0, outp))
for W in (128, 256):
    for pfx in ('base', 'fix'):
        ex = sum(r['%s%d' % (pfx, W)]['exact'] for r in recs)
        e = np.mean([r['%s%d' % (pfx, W)]['reads'] for r in recs]) / n
        print('  %s W2=%3d exact %3d/%3d  e/n %.3f'
              % (pfx, W, ex, len(recs), e))
