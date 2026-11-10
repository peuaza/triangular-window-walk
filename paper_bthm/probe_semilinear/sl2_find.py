# -*- coding: utf-8 -*-
"""sl2_find -- locate defect tapes (hash % n == 0) and run the defect branch on them."""
import os, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); PAPER = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(PAPER, 'probe_class_rerun')); sys.path.insert(0, HERE)
import pcr_lib as P, sl2_engine as E

n = int(sys.argv[1]) if len(sys.argv) > 1 else 128
want = int(sys.argv[2]) if len(sys.argv) > 2 else 2
batch, found, tried = 4000, 0, 0
seed = 7
while found < want and tried < 400000:
    X = P.make_tapes(n, batch, seed + tried // batch)
    Z = E.tape_hashes(X, n) % n
    for t in np.nonzero(Z == 0)[0]:
        S = P.St(X[t]); out, steps, msum, inc, cuts, end = E.Runner2(S).run('ballast')
        print('DEFECT n=%d tried=%d  end=%d out=%.4f steps=%d msum=%d inc=%d (%.3f n^2) cuts=%d'
              % (n, tried+t, end, out, steps, msum, inc, inc/float(n*n), cuts), flush=True)
        found += 1
        if found >= want: break
    tried += batch
print('defect fraction so far: %d/%d = %.5f   (1/n = %.5f)' % (found, tried, found/float(tried), 1.0/n))
