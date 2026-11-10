# -*- coding: utf-8 -*-
"""pf_gate -- the three checks PEAK_FIX.md rests on.

1. pf_lib.jm_base_audit == el_cell.jm_walk on (end, reads, epochs, cuts, eg),
   6 W2 values, 15 meeting tapes at n = 512 and 1024.
2. The fixed machine is a no-op on (end, delta, exact) over every segment.
3. The 1 500-tape exactness totals at W2 = 128 and 256.

Usage: python pf_gate.py
"""
import glob
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(os.path.dirname(HERE), 'probe_en_linear'))

import numpy as np

import pf_lib as F
import el_cell as C
import el_lib as E


def gate1():
    bad = cells = 0
    for n in (512, 1024):
        i = m = 0
        while m < 15 and i < 20000:
            X = E.tape(n, i)
            M = E.Mach(X, n)
            t = M.trace('peak')
            i += 1
            if t['end'] != 0:
                continue
            m += 1
            for W in (8, 16, 32, 64, 128, 256):
                d, en, rd, epc, cu, eg, nsc = C.jm_walk('peak', M.S, n, W, 'sliding')
                b = F.jm_base_audit('peak', M.S, n, W)
                cells += 1
                if (int(en) != int(b['end']) or rd != b['reads']
                        or epc != b['epochs'] or cu != b['cuts'] or eg != b['eg']):
                    bad += 1
    print('gate 1: %d cells, %d mismatches' % (cells, bad))
    return bad == 0


def gate2():
    ok = True
    for W in (128, 256):
        ne = nd = nx = tot = 0
        for f in sorted(glob.glob(os.path.join(HERE, 'seg_*.jsonl'))):
            for line in open(f):
                r = json.loads(line)
                b, x = r['base%d' % W], r['fix%d' % W]
                tot += 1
                ne += (b['end'] != x['end'])
                nd += (abs(b['delta'] - x['delta']) > 1e-12)
                nx += (b['exact'] != x['exact'])
        print('gate 2 W2=%d: %d tapes, end diffs %d, delta diffs %d, exact diffs %d'
              % (W, tot, ne, nd, nx))
        ok = ok and (ne == nd == nx == 0)
    return ok


def gate3():
    tot = {p: 0 for p in ('base128', 'fix128', 'base256', 'fix256')}
    k = 0
    for f in sorted(glob.glob(os.path.join(HERE, 'seg_*.jsonl'))):
        if '4096_900000' in f:
            continue
        for line in open(f):
            r = json.loads(line)
            k += 1
            for p in tot:
                tot[p] += r[p]['exact']
    print('gate 3: %d tapes  base128 %d  fix128 %d  base256 %d  fix256 %d'
          % (k, tot['base128'], tot['fix128'], tot['base256'], tot['fix256']))
    return k == 1500


if __name__ == '__main__':
    a = gate1()
    b = gate2()
    c = gate3()
    print('GATE: %s' % ('all pass' if (a and b and c) else 'FAIL'))
