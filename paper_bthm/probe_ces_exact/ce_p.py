# -*- coding: utf-8 -*-
"""ce_p -- p and e on the SAME meeting tapes, same arm, same accounting.

Reads ce_patch_<rule>_<n>.json (per-tape rows, corrected opening pass) and
compares with probe_table_p's pt_<rule>_<n>.json.

p protocol (pt_cell.py): population = meeting tapes; axis = the arm's mean
per-tape spend; a non-finite prediction is scored at the population mean; p is
the smallest spend whose population R^2 is >= 0.99.  For an arm of fixed cost
that is simply its mean spend, once its R^2 clears the threshold.

usage: ce_p.py
"""
import glob
import json
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
PT = os.path.join(os.path.dirname(HERE), 'probe_table_p')


def r2(y, p):
    y = np.asarray(y, float); p = np.asarray(p, float)
    return 1.0 - ((y - p) ** 2).sum() / ((y - y.mean()) ** 2).sum()


def cell(fn):
    J = json.load(open(fn))
    n = J['n']; rule = J['rule']
    rec = None
    ptf = os.path.join(PT, 'pt_%s_%d.json' % (rule, n))
    if os.path.exists(ptf):
        rec = json.load(open(ptf))
    print('\n%s n=%d  (record: p=%.1f = %.2f n, e_meet=%.0f = %.2f n)'
          % (rule, n, rec['p'] if rec else -1,
             (rec['p'] / n) if rec else -1,
             rec['e_meet'] if rec else -1,
             (rec['e_meet'] / n) if rec else -1))
    print('  arm      e_meet/n   R2(meet)  exact_meet  p/n   exact?')
    arms = []
    for name in ('ces', 'w2_32', 'w2_64', 'w2_128', 'w2_256', 'w2_512',
                 'w2_q', 'w2_half', 'w2_3q', 'w2_n', 'oracle'):
        v = J['machines'].get(name)
        if not v:
            continue
        rows = [r for r in v['rows'] if r['tag'] == 'meet']
        y = np.array([r['truth'] for r in rows], float)
        g = np.array([np.nan if r['got'] is None else r['got'] for r in rows], float)
        mu = y.mean()
        pred = np.where(np.isfinite(g), g, mu)
        rr = r2(y, pred)
        spend = np.mean([r['reads'] for r in rows]) / n
        arms.append((name, spend, rr, v['exact_meet'], v['n_meet']))
    legal = [a for a in arms if a[0] != 'oracle']
    p_min = min([a[1] for a in legal if a[2] >= 0.99], default=None)
    e_min = min([a[1] for a in legal if a[3] == a[4]], default=None)
    for (name, spend, rr, em, nm) in arms:
        print('  %-8s %8.2f   %8.5f  %4d/%-4d  %s  %s'
              % (name, spend, rr, em, nm,
                 ('%.2f' % spend) if rr >= 0.99 else '-',
                 'yes' if em == nm else ''))
    print('  -> p/n = %s, e/n = %s, p<=e %s'
          % (('%.2f' % p_min) if p_min else '-', ('%.2f' % e_min) if e_min else '-',
             ('OK' if (p_min is not None and e_min is not None and p_min <= e_min + 1e-9)
              else 'VIOLATION')))


def main():
    for fn in sorted(glob.glob(os.path.join(HERE, 'ce_patch_*.json'))):
        cell(fn)


if __name__ == '__main__':
    main()
