# -*- coding: utf-8 -*-
"""el_split -- mispredicted epochs below and above the buffer, per rule.

The distinction that decides the machine: a failure at M0 <= W2 is absorbed by
the endgame, one above is not.
"""
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import el_cell as C
import el_lib as E

NS = (512, 1024, 2048, 4096)
W = 256


def main():
    print('%-9s %8s %8s %8s %6s %6s' % ('rule', 'below', 'above', 'total',
                                        'below%', 'tapes>W'))
    for rule in E.RULES:
        lo = hi = tapes = tot = 0
        for n in NS:
            p = os.path.join(HERE, 'el_%s_%d.json' % (C.safe(rule), n))
            if not os.path.exists(p):
                continue
            d = json.load(open(p))
            for t in d['taps']:
                if t['end'] != 0:
                    continue
                for b in t['bad']:
                    if b[0] > W:
                        hi += 1
                    else:
                        lo += 1
                if any(b[0] > W for b in t['bad']):
                    tapes += 1
        tot = lo + hi
        print('%-9s %8d %8d %8d %5.1f%% %6d' % (
            rule, lo, hi, tot, 100.0 * lo / max(1, tot), tapes))


if __name__ == '__main__':
    main()
