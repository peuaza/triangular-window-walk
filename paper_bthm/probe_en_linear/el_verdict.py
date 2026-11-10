# -*- coding: utf-8 -*-
"""el_verdict -- the verdict table of EN_LINEAR.md.

For each rule: the device class, the smallest W2 that is exact on every meeting
tape at every size that has been run, the exact counts, e/n, the slope of e/n
over n, and how many meeting tapes had at least one mispredicted epoch above
that W2 (the difference between exact-by-construction and exact-by-luck).

Usage: python el_verdict.py
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


def main():
    print('%-9s %-7s %5s %-28s %-26s %6s %6s %s' % (
        'rule', 'device', 'W2', 'exact / nmeet by n', 'e/n by n', 'slope',
        'tapes_bad', 'kind'))
    rows = []
    for rule in E.RULES:
        cells = {}
        for n in NS:
            p = os.path.join(HERE, 'el_%s_%d.json' % (C.safe(rule), n))
            if os.path.exists(p):
                cells[n] = json.load(open(p))
        if not cells:
            continue
        ns = sorted(cells)
        mode = cells[ns[0]]['mode']
        best = None
        for W in W2S:
            w = str(W)
            if all(sum(t['mach'][w]['exact'] for t in cells[n]['taps']
                       if t['end'] == 0) ==
                   sum(1 for t in cells[n]['taps'] if t['end'] == 0) for n in ns):
                best = W
                break
        exact_all = best is not None
        if best is None:
            tot = {W: sum(sum(t['mach'][str(W)]['exact'] for t in cells[n]['taps']
                              if t['end'] == 0) for n in ns) for W in W2S}
            best = max(tot, key=tot.get)
        w = str(best)
        ex = []
        en = []
        tb = 0
        struct = True
        for n in ns:
            meet = [t for t in cells[n]['taps'] if t['end'] == 0]
            ex.append('%d/%d' % (sum(t['mach'][w]['exact'] for t in meet), len(meet)))
            en.append(np.mean([t['mach'][w]['reads'] for t in meet]) / n)
            b = sum(1 for t in meet if any(x[0] > best for x in t['bad']))
            tb += b
            struct = struct and (b == 0)
        slope = float(np.polyfit(np.log(np.array(ns, float)),
                                 np.log(np.array(en, float)), 1)[0])
        kind = ('struct' if struct else ('output' if exact_all else 'inexact'))
        print('%-9s %-7s %5d %-28s %-26s %+6.2f %6d %s' % (
            rule, mode, best, ' '.join(ex),
            ' '.join('%.1f' % x for x in en), slope, tb, kind))
        rows.append(dict(rule=rule, mode=mode, w2=best, exact=ex, e=en,
                         slope=slope, tapes_bad=tb, kind=kind))
    with open(os.path.join(HERE, 'el_verdict.json'), 'w') as fh:
        json.dump(rows, fh, indent=1)


if __name__ == '__main__':
    main()
