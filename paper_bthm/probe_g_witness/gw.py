"""A witness panel for g: structured tapes on which the walk degenerates.

The ensemble's tapes are pairs of balanced columns X0, X1 in {0,1}^n (each with
exactly n/2 ones); sigma = X0 + X1 in {0,1,2}.  Every witness here is built as a
pair of balanced columns, so it lies in the ensemble's support (probability > 0).

For each rule and each witness we record the machine's reads, whether the tape
meets (delta defined) and delta, then compare with the sample maximum g_K from
MASTER_FINAL / probe_k300.

Usage:  python gw.py [n ...]
"""
import json
import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
PAPER = os.path.dirname(HERE)
TAN = os.path.join(PAPER, 'probe_triple_nat')
if TAN not in sys.path:
    sys.path.insert(0, TAN)

import tr_lib as T                                              # noqa: E402


def _cols(n, spec):
    """Build the pair of balanced columns from a spec string.

    'zero'    : 0^n                       (not in the support; a control)
    'one'     : all 1
    'two'     : all 2
    'h2'      : 2^{n/2} 0^{n/2}          (the longrun witness shape)
    'h0'      : 0^{n/2} 2^{n/2}
    'alt01'   : 0,1,0,1...
    'alt02'   : 0,2,0,2...
    'blk'     : blocks of b: 2^b 0^b     (b = n//8)
    'half2a'  : 1^{n/2} 2^{n/4} 0^{n/4}  (a three-block split)
    'odds2'   : X0 = X1 = 1 on the odd positions
    """
    h = n // 2
    X0 = np.zeros(n, np.int8)
    X1 = np.zeros(n, np.int8)
    if spec == 'zero':
        pass
    elif spec == 'one':
        X0[:h] = 1
        X1[h:] = 1
    elif spec == 'two':
        X0[:] = 1
        X1[:] = 1
    elif spec == 'h2':
        X0[:h] = 1
        X1[:h] = 1
    elif spec == 'h0':
        X0[h:] = 1
        X1[h:] = 1
    elif spec == 'alt01':
        X0[1::2] = 1
    elif spec == 'alt02':
        X0[1::2] = 1
        X1[1::2] = 1
    elif spec == 'blk':
        b = max(2, n // 8)
        for st in range(0, n, 2 * b):
            X0[st:st + b] = 1
            X1[st:st + b] = 1
    elif spec == 'half2a':
        q = n // 4
        X0[:h] = 1
        X1[:q] = 1
    elif spec == 'odds2':
        X0[1::2] = 1
        X1[1::2] = 1
    else:
        raise KeyError(spec)
    # keep every column balanced: n/2 ones
    for X in (X0, X1):
        ones = np.flatnonzero(X)
        need = h - ones.size
        if need > 0:
            free = np.flatnonzero(X == 0)
            X[free[:need]] = 1
        elif need < 0:
            X[ones[:(-need)]] = 0
    a = (X0 + X1).astype(np.int64)
    s = (X0 - X1).astype(np.int64)
    return a, s, np.stack([X0, X1], 1)


SPECS = ('one', 'two', 'h2', 'h0', 'alt01', 'alt02', 'blk', 'half2a', 'odds2')

RULES = ('monotone', 'longrun', 'pattern', 'runs', 'dense2', 'peak', 'trip', 'oddrun',
         'canonical', 'pal', 'hashmod', 'MAXP', 'LAS', 'KCHG', 'BAL>m/2', 'BAL>m/4',
         'NBLK=1', 'NBLK>2', 'NBLK>3', 'AMAX', 'AMIN')


def main(ns):
    path = os.path.join(HERE, 'gw_cells.json')
    out = json.load(open(path)) if os.path.exists(path) else {}
    for n in ns:
        rows = {}
        for spec in SPECS:
            a, s, X = _cols(n, spec)
            for rule in RULES:
                try:
                    r = T.run(rule, a, s)
                except Exception as e:                          # noqa: BLE001
                    rows[f'{rule}|{spec}'] = {'err': repr(e)[:80]}
                    continue
                rows[f'{rule}|{spec}'] = {
                    'reads': int(r.get('reads') or 0),
                    'delta': (None if r.get('delta') is None or r.get('delta') != r.get('delta')
                              else float(r['delta'])),
                    'end': r.get('end'),
                }
        out[str(n)] = rows
        with open(os.path.join(HERE, 'gw_cells.json'), 'w') as f:
            json.dump(out, f)
        print(f'n={n} done', flush=True)
    return out


if __name__ == '__main__':
    ns = [int(x) for x in sys.argv[1:]] or [512, 1024, 2048]
    main(ns)
