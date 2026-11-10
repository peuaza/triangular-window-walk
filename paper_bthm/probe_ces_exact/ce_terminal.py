# -*- coding: utf-8 -*-
"""ce_terminal -- terminal-epoch M0 by outcome class (meet vs halt).

For each tape, run the exact walk and record the LAST epoch's M0, its outcome
('meet' -> the walk answers; 'cycle'/'narrow' -> it halts with no output) and
its length.  This is the quantity that decides whether a W2 threshold can be
exact: meeting tapes must be caught to get delta right, halting tapes must be
caught to avoid inventing a delta.

usage: ce_terminal.py <rule> <n> [nt_each]
Writes ce_terminal_<rule>_<n>.json
"""
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
for _p in (HERE, os.path.join(ROOT, 'probe_e_uniform'),
           os.path.join(ROOT, 'probe_kappa_pred'),
           os.path.join(ROOT, 'probe_ologn2')):
    if _p not in sys.path:
        sys.path.insert(0, _p)

import ce_lib as CL                                              # noqa: E402
import eu_lib as EL                                              # noqa: E402


def q(v):
    v = np.array(sorted(v), float)
    if len(v) == 0:
        return {}
    return dict(n=int(len(v)), min=float(v.min()), med=float(np.median(v)),
                p90=float(np.percentile(v, 90)), p99=float(np.percentile(v, 99)),
                max=float(v.max()))


def main():
    rule = sys.argv[1]; n = int(sys.argv[2])
    nt = int(sys.argv[3]) if len(sys.argv) > 3 else 100
    out = dict(rule=rule, n=n, nt=nt,
               meet=[], halt=[], meet_ln=[], halt_ln=[], halt_type={})
    i = 0
    while (len(out['meet']) < nt or len(out['halt']) < nt) and i < 40 * nt:
        a, s, X = EL.make_tape(n, 5000 + n + 100003 * i)
        a = np.asarray(a, np.int64)
        d, why, eps = CL.trace(a, rule)
        M0, b2, b1, oc, ln, L, R, ee = eps[-1]
        if why == 'meet':
            out['meet'].append(int(M0)); out['meet_ln'].append(int(ln))
        else:
            out['halt'].append(int(M0)); out['halt_ln'].append(int(ln))
            out['halt_type'][oc] = out['halt_type'].get(oc, 0) + 1
        i += 1
    out['scanned'] = i
    with open(os.path.join(HERE, 'ce_terminal_%s_%d.json' % (rule, n)), 'w') as f:
        json.dump(out, f, indent=1)
    print('%s n=%d scanned=%d' % (rule, n, i))
    print('  meet  M0 %s' % q(out['meet']))
    print('  meet  ln %s' % q(out['meet_ln']))
    print('  halt  M0 %s  types %s' % (q(out['halt']), out['halt_type']))
    print('  halt  ln %s' % q(out['halt_ln']))


if __name__ == '__main__':
    main()
