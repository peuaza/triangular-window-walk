"""Per-tape costs for the figure, restricted to machines whose numbers the paper already quotes.

monotone: the counter machine of Theorem thm:collapse (2n-1 reads on every tape).
peak:     the W2=128 jump machine of Section sec:cost (2.03 n reads on the meeting tapes).
The band rules come from probe_k300 (the e-best band machine).  pattern and runs are deliberately
absent: the paper's 2.57 n and 2.77 n come from the CES variant used in probe_ces_exact, whose
per-tape reads are not in any stored array, and a figure that silently used the W=16 variant would
contradict the text.
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

N, NT = 4096, 100
out = {}

vals = []
for i, a, s, X in T.tapes(N, nt=NT):
    r = T.run('monotone', a, s)
    if r.get('end', 1) == 0:
        vals.append(float(r['reads']))
out['monotone'] = vals
print('monotone', len(vals), np.mean(vals) / N if vals else 0)

vals = []
for i, a, s, X in T.tapes(N, nt=NT):
    Xc = E.tape(N, i)
    M = E.Mach(Xc, N)
    m = F.jm_fixed('peak', M.S, N, 128)
    if m.get('end', 1) == 0:
        vals.append(float(m['reads']))
out['peak'] = vals
print('peak', len(vals), np.mean(vals) / N if vals else 0)

json.dump(out, open(os.path.join(HERE, 'fl_costs.json'), 'w'), indent=1)
print('wrote fl_costs.json')
