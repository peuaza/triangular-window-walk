"""Per-tape costs of the four fully-linear rules at n=4096, for the difficulty-density figure.

Writes probe_cost_dist/fl_costs.json:  {rule: [reads per tape, ...]}  (300 tapes, meeting tapes only
for the rules whose machine can report a halt; `peak` uses the W2=128 machine of the paper).
"""
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
PAPER = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(PAPER, 'probe_triple_nat'))
sys.path.insert(0, os.path.join(PAPER, 'probe_peak_fix'))
sys.path.insert(0, os.path.join(PAPER, 'probe_en_linear'))

import tr_lib as T                      # noqa: E402
import pf_lib as F                      # noqa: E402
import el_lib as E                      # noqa: E402

N = 4096
out = {}
for rule in ('monotone', 'pattern', 'runs'):
    vals = []
    for i, a, s, X in T.tapes(N, nt=100):
        r = T.run(rule, a, s)
        if r.get('end', 1) == 0:        # 0 = halted with an output, 1 = divergence
            vals.append(float(r['reads']))
    out[rule] = vals
    print(rule, len(vals), np.mean(vals) if vals else 0)

vals = []
for i, a, s, X in T.tapes(N, nt=100):
    Xc = E.tape(N, i)
    M = E.Mach(Xc, N)
    m = F.jm_fixed('peak', M.S, N, 128)
    if m.get('end', 1) == 0:
        vals.append(float(m['reads']))
out['peak'] = vals
print('peak', len(vals), np.mean(vals) if vals else 0)

os.makedirs(HERE, exist_ok=True)
json.dump(out, open(os.path.join(HERE, 'fl_costs.json'), 'w'), indent=1)
print('wrote', os.path.join(HERE, 'fl_costs.json'))
