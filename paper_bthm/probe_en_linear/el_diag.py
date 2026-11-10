# -*- coding: utf-8 -*-
"""el_diag -- where the prediction fails: M0 of the bad epochs, per rule."""
import json
import sys

import numpy as np

for f in sys.argv[1:]:
    d = json.load(open(f))
    n = d['n']
    meet = [t for t in d['taps'] if t['end'] == 0]
    bad = [b for t in meet for b in t['bad']]
    print('==', d['rule'], n, 'meet', len(meet), 'bad', len(bad),
          'tapes with >=1 bad', sum(1 for t in meet if t['bad']))
    if not bad:
        continue
    m0 = np.array([b[0] for b in bad])
    ln = np.array([b[5] for b in bad])
    kind = {}
    for b in bad:
        k = 'terminal' if b[2] in (0, 4, 5) else 'cut'
        kind[k] = kind.get(k, 0) + 1
    print('   M0 of bad: med %d p10 %d p90 %d max %d ; %s'
          % (np.median(m0), np.percentile(m0, 10), np.percentile(m0, 90),
             m0.max(), kind))
    print('   epoch length of bad: med %d p10 %d p90 %d  (<=4: %d, <=16: %d)'
          % (np.median(ln), np.percentile(ln, 10), np.percentile(ln, 90),
             int((ln <= 4).sum()), int((ln <= 16).sum())))
    for W in d['w2s']:
        above = [b for b in bad if b[0] > W]
        nT = sum(1 for t in meet if any(b[0] > W for b in t['bad']))
        print('     W2=%4d  bad above %5d  tapes with bad above %3d/%d'
              % (W, len(above), nT, len(meet)))
