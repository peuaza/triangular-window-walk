# -*- coding: utf-8 -*-
"""ce_summary -- tables for CES_EXACT.md, straight from the JSONs."""
import glob
import json
import os
from collections import Counter

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))


def key(fn):
    p = os.path.basename(fn).replace('.json', '').split('_')
    return (p[-2], int(p[-1]))


def t1():
    print('## per-epoch failure of the CES outcome rule (true trajectory)\n')
    print('| rule | n | tapes | meet | epochs | cut-mism | cut rate | early-term >W | term rate | rate*n | tapes w/ term | max M0 | med ln |')
    print('|---|---|---|---|---|---|---|---|---|---|---|---|---|')
    for fn in sorted(glob.glob(os.path.join(HERE, 'ce_fail_*.json')), key=key):
        J = json.load(open(fn))
        n = J['n']; nt = J['nt']
        vals = []
        for m, c in J['term_M0']:
            vals += [int(m)] * c
        vals = sorted(vals)
        # median terminal length is not stored; use mean
        print('| %s | %d | %d | %d | %d | %d | %.1e | %d | %.2e | %.2f | %d (%.0f%%) | %d | %.0f |'
              % (J['rule'], n, nt, J['meet_tapes'], J['tot_epochs'],
                 J['n_cut_mismatch'],
                 J['n_cut_mismatch'] / max(1, J['tot_epochs']),
                 J['n_term'] - (J['n_term'] - len(vals)) if False else len(vals),
                 len(vals) / max(1, J['tot_epochs']),
                 len(vals) / max(1, J['tot_epochs']) * n,
                 len(vals), 100 * len(vals) / nt,
                 max(vals) if vals else 0, J['term_ln_mean'] or 0))
    print()


def t2():
    print('## first failing epoch: position and size\n')
    print('| rule | n | tapes failing | median j/n | p90 j/n | median M0 | M0<=14 | 15-64 | 65-256 | >256 |')
    print('|---|---|---|---|---|---|---|---|---|---|')
    for fn in sorted(glob.glob(os.path.join(HERE, 'ce_fail_*.json')), key=key):
        J = json.load(open(fn))
        fi = np.array(J['first_idx_frac']); fm = np.array(J['first_M0'])
        if len(fi) == 0:
            continue
        b = [int((fm <= 14).sum()), int(((fm > 14) & (fm <= 64)).sum()),
             int(((fm > 64) & (fm <= 256)).sum()), int((fm > 256).sum())]
        print('| %s | %d | %d/%d | %.2f | %.2f | %d | %d | %d | %d | %d |'
              % (J['rule'], J['n'], len(fi), J['nt'], np.median(fi),
                 np.percentile(fi, 90), np.median(fm), b[0], b[1], b[2], b[3]))
    print()


def t3():
    print('## machines: the three counts that matter (200 meeting + 100 halting tapes)\n')
    print('(a) true meeting tapes with the right delta; (b) true halting tapes the machine also')
    print('calls halting; (c) true meeting tapes the machine wrongly calls halting (declines);')
    print('(d) true halting tapes the machine gives a delta for ("invented"); (e) meeting tapes')
    print('with a wrong delta.  e needs a = n_meet and c = 0.\n')
    print('| rule | n | machine | e/n meet | (a) exact meet | (b) halt correct | (c) meet halted | (d) delta on halt | (e) wrong delta |')
    print('|---|---|---|---|---|---|---|---|---|')
    order = ['ces', 'k4', 'w2_32', 'w2_64', 'w2_128', 'w2_256', 'w2_512',
             'w2_q', 'w2_half', 'w2_3q', 'w2_n', 'oracle']
    for fn in sorted(glob.glob(os.path.join(HERE, 'ce_patch_*.json')), key=key):
        J = json.load(open(fn))
        M = J['machines']
        for name in order:
            if name not in M:
                continue
            v = M[name]
            cm = v['counts_meet']; ch = v['counts_halt']
            print('| %s | %d | %s | %.2f | %d/%d | %d/%d | %d | %d | %d |'
                  % (J['rule'], J['n'], name,
                     v['reads_meet_mean'] / J['n'],
                     cm.get('exact', 0), v['n_meet'],
                     ch.get('both_nan', 0), sum(ch.values()),
                     cm.get('declined', 0),
                     ch.get('invented', 0),
                     cm.get('wrong', 0)))
    print()


def t4():
    print('## CES output classes on all tapes (random T_n)\n')
    print('| rule | n | exact | both_nan | invented | declined | wrong |')
    print('|---|---|---|---|---|---|---|')
    for fn in sorted(glob.glob(os.path.join(HERE, 'ce_align_*.json')), key=key):
        J = json.load(open(fn))
        c = J['counts']['cheap']
        n = sum(c.values())
        print('| %s | %d | %d | %d | %d | %d | %d |'
              % (J['rule'], J['n'], c.get('exact', 0), c.get('both_nan', 0),
                 c.get('invented', 0), c.get('declined', 0), c.get('wrong_delta', 0)))
    print()


if __name__ == '__main__':
    t1(); t2(); t3(); t4()
