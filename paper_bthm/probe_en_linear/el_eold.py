# -*- coding: utf-8 -*-
"""el_eold -- the record's own e (mean reads over meeting tapes / n), from the
gated five-size cells of probe_table_eg.  One number per (rule, n), for the
comparison column of EN_LINEAR.md.

Also writes the e/n over ALL draws, which is what the pcr table's means are.
"""
import json
import os

import numpy as np

PAPER = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
HERE = os.path.dirname(os.path.abspath(__file__))
TN = os.path.join(PAPER, 'probe_triple_nat')

SAFE = {'BAL>m/2': 'BALgmd2', 'BAL>m/4': 'BALgmd4', 'NBLK=1': 'NBLKe1',
        'NBLK>2': 'NBLKg2', 'NBLK>3': 'NBLKg3'}
RULES = ['dense2', 'peak', 'trip', 'oddrun', 'canonical', 'MAXP', 'LAS', 'KCHG',
         'BAL>m/2', 'BAL>m/4', 'NBLK=1', 'NBLK>2', 'NBLK>3', 'AMAX', 'AMIN']
NS = (512, 1024, 2048, 4096)


def main():
    out = {}
    for r in RULES:
        s = SAFE.get(r, r)
        for n in NS:
            p = os.path.join(PAPER, 'probe_table_eg', 'teg_cell_%s_%d.json' % (s, n))
            src = 'teg'
            if not os.path.exists(p):
                p = os.path.join(TN, 'tr_%s_%d.json' % (s, n))
                src = 'tr'
            if not os.path.exists(p):
                continue
            d = json.load(open(p))
            rd = np.array(d['reads'], float)
            e = np.array(d['end'])
            meet = (e == 0)
            out['%s|%d' % (r, n)] = dict(
                src=src, ndraw=len(rd), nmeet=int(meet.sum()),
                e_meet=float(rd[meet].mean() / n) if meet.any() else None,
                e_all=float(rd.mean() / n),
                machine=d.get('machine'))
    with open(os.path.join(HERE, 'el_eold.json'), 'w') as fh:
        json.dump(out, fh, indent=1)
    print('%-9s %5s %6s %7s %9s %9s  %s' % ('rule', 'n', 'src', 'nmeet',
                                            'e_meet/n', 'e_all/n', 'machine'))
    for r in RULES:
        for n in NS:
            k = '%s|%d' % (r, n)
            if k not in out:
                continue
            v = out[k]
            print('%-9s %5d %6s %7d %9s %9.1f  %s' % (
                r, n, v['src'], v['nmeet'],
                '-' if v['e_meet'] is None else '%.1f' % v['e_meet'],
                v['e_all'], v['machine']))


if __name__ == '__main__':
    main()
