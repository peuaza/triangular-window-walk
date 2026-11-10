# -*- coding: utf-8 -*-
"""el_lawm -- the epoch law by region width: does a W2 buffer repair it?

For bucket (b2,b1) = (1,0) and (1,1) tabulate the outcome accuracy against M0 in
three bands.  If the errors live at small M0 the buffer repairs them; if they are
spread over the whole tape it cannot.

Usage: python el_lawm.py [NT] [N]
Writes el_lawm.json.
"""
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import el_cell as C
import el_lib as E

NT = int(sys.argv[1]) if len(sys.argv) > 1 else 40
N = int(sys.argv[2]) if len(sys.argv) > 2 else 1024
BANDS = ((3, 64), (64, 256), (256, 1 << 30))


def main():
    out = {}
    print('rule       bucket  band(M0)      n   acc')
    for rule in E.RULES:
        tab = {}
        for i in range(NT):
            X = E.tape(N, i)
            M = E.Mach(X, N)
            t = M.trace(rule)
            for k in range(t['nep'] - 1):
                if t['epb2'][k] < 0 or t['eplen'][k] == 0:
                    continue
                key = (int(t['epb2'][k]), int(t['epb1'][k]))
                m0 = int(t['epM'][k])
                o = int(t['epout'][k])
                for j, (lo, hi) in enumerate(BANDS):
                    if lo <= m0 < hi:
                        tab.setdefault((key, j), {}).setdefault(o, 0)
                        tab[(key, j)][o] += 1
                        break
        row = {}
        for key in ((1, 0), (1, 1), (0, 0)):
            for j, (lo, hi) in enumerate(BANDS):
                v = tab.get((key, j))
                if not v:
                    continue
                tot = sum(v.values())
                pred = 2 if key[1] else 1
                acc = v.get(pred, 0) / tot
                row['%d%d_%d' % (key[0], key[1], j)] = dict(
                    n=tot, acc=acc, out=v)
                print('%-9s %d%d      [%d,%d)  %6d   %.4f'
                      % (rule, key[0], key[1], lo, hi, tot, acc))
        out[rule] = row
    with open(os.path.join(HERE, 'el_lawm.json'), 'w') as fh:
        json.dump(out, fh, indent=1)
    print('wrote el_lawm.json  (n=%d, %d tapes)' % (N, NT))


if __name__ == '__main__':
    main()
