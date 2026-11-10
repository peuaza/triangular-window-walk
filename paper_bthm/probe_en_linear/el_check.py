# -*- coding: utf-8 -*-
"""el_check -- verify the numbers quoted in EN_LINEAR.md against the cells."""
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import el_cell as C
import el_lib as E

NS = (512, 1024, 2048, 4096)


def main():
    print('%-9s %-26s %-26s %6s %6s %6s %6s' % (
        'rule', 'e/n @W2=256', 'exact @W2=256', 'badM', 'badH', 'tapM', 'tapH'))
    for rule in E.RULES:
        en, ex = [], []
        bm = bh = tm = th = 0
        for n in NS:
            p = os.path.join(HERE, 'el_%s_%d.json' % (C.safe(rule), n))
            if not os.path.exists(p):
                en.append('--')
                ex.append('--')
                continue
            d = json.load(open(p))
            meet = [t for t in d['taps'] if t['end'] == 0]
            oth = [t for t in d['taps'] if t['end'] != 0]
            en.append('%.1f' % (np.mean([t['mach']['256']['reads'] for t in meet]) / n))
            ex.append('%d/%d' % (sum(t['mach']['256']['exact'] for t in meet),
                                 len(meet)))
            bm += sum(len(t['bad']) for t in meet)
            bh += sum(len(t['bad']) for t in oth)
            tm += sum(1 for t in meet if t['bad'])
            th += sum(1 for t in oth if t['bad'])
        print('%-9s %-26s %-26s %6d %6d %6d %6d' % (
            rule, ' '.join(en), ' '.join(ex), bm, bh, tm, th))


if __name__ == '__main__':
    main()
