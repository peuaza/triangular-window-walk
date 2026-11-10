# -*- coding: utf-8 -*-
"""el_peak3 -- a third independent segment for peak at n = 2048.

The headline number is peak's exactness.  Main grid: 200/200 at n = 2048.
Replay at i = 700 000: 199/200.  This runs a third segment (i = 300 000) to see
whether the failure rate is ~0.5 % or the second sample was unlucky.

Usage: python el_peak3.py [N] [OFF] [NMEET]
"""
import sys
import time

import numpy as np

sys.path.insert(0, '.')
import el_cell as C
import el_lib as E

n = int(sys.argv[1]) if len(sys.argv) > 1 else 2048
off = int(sys.argv[2]) if len(sys.argv) > 2 else 300000
nmeet = int(sys.argv[3]) if len(sys.argv) > 3 else 200

i = off
recs = []
t0 = time.time()
while len(recs) < nmeet and i < off + 300000:
    X = E.tape(n, i)
    M = E.Mach(X, n)
    t = M.trace('peak')
    if t['end'] == 0:
        out = {}
        for W in (64, 128, 256):
            d, en, rd, epc, cu, eg, nsc = C.jm_walk('peak', M.S, n, W, 'sliding')
            ok = (int(en) == int(t['end'])) and (
                (d is None and t['delta'] < 0) or
                (d is not None and t['delta'] >= 0 and abs(d - t['delta']) < 1e-12))
            out[str(W)] = (int(ok), int(rd))
        recs.append((i, out))
    i += 1
print('peak n=%d offset=%d meetings=%d draws=%d secs=%.0f'
      % (n, off, len(recs), i - off, time.time() - t0))
for W in (64, 128, 256):
    ex = sum(r[1][str(W)][0] for r in recs)
    e = np.mean([r[1][str(W)][1] for r in recs]) / n
    print('   W2=%3d exact %3d/%3d  e/n %.2f' % (W, ex, len(recs), e))
