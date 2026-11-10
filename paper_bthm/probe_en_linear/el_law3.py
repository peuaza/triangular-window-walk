# -*- coding: utf-8 -*-
"""el_law3 -- the epoch law on MEETING tapes vs HALTING tapes, by M0 band.

The load-bearing measurement of EN_LINEAR.md §3.  For every interior epoch:
  bucket = (b2,b1) at the epoch start
  ok     = the epoch's outcome equals the machine's prediction
  band   = M0 < 256 or >= 256
  kind   = the tape meets or halts
The claim is that for the level/extremum rules the law is exact on meeting tapes
above the band where the buffer sits, and breaks on halting tapes.

Usage: python el_law3.py [NT] [N]
Writes el_law3.json.
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


def main():
    out = {}
    print('%-9s %-4s %-5s %-8s %8s %8s' % ('rule', 'bkt', 'kind', 'band', 'n',
                                           'acc'))
    for rule in E.RULES:
        tab = {}
        for i in range(NT):
            X = E.tape(N, i)
            M = E.Mach(X, N)
            t = M.trace(rule)
            kind = 'meet' if t['end'] == 0 else 'halt'
            for k in range(t['nep'] - 1):
                if t['epb2'][k] < 0 or t['eplen'][k] == 0:
                    continue
                key = (int(t['epb2'][k]), int(t['epb1'][k]))
                m0 = int(t['epM'][k])
                band = 'lo' if m0 < 256 else 'hi'
                ok = 1 if int(t['epout'][k]) == C.pred_of(rule, key[0], key[1]) else 0
                d = tab.setdefault((key, kind, band), [0, 0])
                d[0] += 1
                d[1] += ok
        row = {}
        for key in sorted(tab):
            kk = '%d%d_%s_%s' % (key[0][0], key[0][1], key[1], key[2])
            n, ok = tab[key]
            row[kk] = dict(n=n, ok=ok, acc=ok / n)
            print('%-9s %d%d   %-5s %-8s %8d %8.4f'
                  % (rule, key[0][0], key[0][1], key[1], key[2], n, ok / n))
        out[rule] = row
    with open(os.path.join(HERE, 'el_law3.json'), 'w') as fh:
        json.dump(out, fh, indent=1)
    print('wrote el_law3.json (n=%d, %d tapes)' % (N, NT))


if __name__ == '__main__':
    main()
