# -*- coding: utf-8 -*-
"""sl2_calib -- what pcr_lib's two accountings actually charge, per tape."""
import os, sys, statistics
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
PAPER = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(PAPER, 'probe_class_rerun'))
import pcr_lib as P

n = int(sys.argv[1]) if len(sys.argv) > 1 else 512
nt = int(sys.argv[2]) if len(sys.argv) > 2 else 20
X = P.make_tapes(n, nt, 7001 + n)
for rule in ('canonical', 'runs', 'pattern'):
    rows = []
    for t in range(nt):
        R = P.Runner(X[t])
        out, steps, msum, inc, cuts, end, snaps = R.run(rule)
        rows.append((out, steps, msum, inc, cuts, end))
    meet = [r for r in rows if r[5] == 0]
    print('rule=%-9s n=%d  meet %d/%d' % (rule, n, len(meet), nt))
    if meet:
        for k, name in ((1, 'steps'), (2, 'msum'), (3, 'inc'), (4, 'cuts')):
            v = [r[k] for r in meet]
            print('   %-6s mean %10.1f  med %10.1f  max %10.1f   /n %8.2f'
                  % (name, np.mean(v), np.median(v), np.max(v), np.mean(v)/n))
        print('   ends:', {e: sum(1 for r in rows if r[5]==e) for e in set(r[5] for r in rows)})
