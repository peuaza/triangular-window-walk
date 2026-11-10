# -*- coding: utf-8 -*-
"""el_law -- the epoch law measured directly: outcome by (b2,b1) at the epoch start.

No machine, no W2: for every interior epoch of every tape, tabulate the epoch's
outcome against the two bits at its start.  This is the evidence behind §3 of
EN_LINEAR.md and it does not depend on any buffer.

Usage: python el_law.py [NT]
Writes el_law.json and prints the table.
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
NS = (512, 1024)
OUTS = (1, 2, 3)


def main():
    out = {}
    for rule in E.RULES:
        out[rule] = {}
        for n in NS:
            tab = {}
            for i in range(NT):
                X = E.tape(n, i)
                M = E.Mach(X, n)
                t = M.trace(rule)
                for k in range(t['nep'] - 1):
                    if t['epb2'][k] < 0 or t['eplen'][k] == 0:
                        continue
                    key = (int(t['epb2'][k]), int(t['epb1'][k]))
                    o = int(t['epout'][k])
                    tab.setdefault(key, {}).setdefault(o, 0)
                    tab[key][o] += 1
            out[rule][str(n)] = {('%d%d' % k): v for k, v in tab.items()}
        print('== %s' % rule)
        for n in NS:
            row = out[rule][str(n)]
            parts = []
            for key in ('00', '01', '10', '11'):
                if key not in row:
                    continue
                v = row[key]
                tot = sum(v.values())
                best = max(v.items(), key=lambda kv: kv[1])
                parts.append('%s n=%d acc=%.4f out=%d' % (
                    key, tot, best[1] / tot, int(best[0])))
            print('   n=%4d  %s' % (n, '   '.join(parts)))
    with open(os.path.join(HERE, 'el_law.json'), 'w') as fh:
        json.dump(out, fh, indent=1)
    print('wrote el_law.json')


if __name__ == '__main__':
    main()
