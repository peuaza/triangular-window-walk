# -*- coding: utf-8 -*-
"""el_grid -- per rule, the W2 x n grid: exactness on meeting tapes and e/n.

Usage: python el_grid.py [RULE ...]
Each cell prints  <exact>/<nmeet> <e/n>.
"""
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import el_cell as C
import el_lib as E

NS = (512, 1024, 2048, 4096)
W2S = (8, 16, 32, 64, 128, 256)


def cell(rule, n):
    p = os.path.join(HERE, 'el_%s_%d.json' % (C.safe(rule), n))
    if not os.path.exists(p):
        return None
    return json.load(open(p))


def main():
    rules = sys.argv[1:] or E.RULES
    for rule in rules:
        ds = {n: cell(rule, n) for n in NS}
        have = [n for n in NS if ds[n]]
        if not have:
            continue
        print('== %-9s mode %s   e_old/n: %s' % (
            rule, ds[have[0]]['mode'],
            '  '.join('%d:%.0f' % (n, np.mean([t['inc'] for t in ds[n]['taps']
                                               if t['end'] == 0]) / n)
                      for n in have)))
        hdr = '   %5s %6s' % ('n', 'meet')
        for W in W2S:
            hdr += ' %16s' % ('W2=%d' % W)
        print(hdr)
        for n in have:
            d = ds[n]
            meet = [t for t in d['taps'] if t['end'] == 0]
            line = '   %5d %6d' % (n, len(meet))
            for W in W2S:
                w = str(W)
                em = sum(t['mach'][w]['exact'] for t in meet)
                e = np.mean([t['mach'][w]['reads'] for t in meet]) / n
                line += ' %9d/%-3d%5.1f' % (em, len(meet), e)
            print(line)
        # halting bonus
        line = '   %5s %6s' % ('halt', '')
        for W in W2S:
            line += ' %16s' % ''
        print('   exact on halting tapes: ' + '  '.join(
            'W%d:%s' % (W, '/'.join(str(sum(t['mach'][str(W)]['exact']
                                            for t in ds[n]['taps'] if t['end']))
                                    for n in have))
            for W in (32, 64, 256)))


if __name__ == '__main__':
    main()
