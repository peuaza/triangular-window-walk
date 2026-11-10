# -*- coding: utf-8 -*-
"""sl2_smoke -- does the ballast rule behave, and does the defect branch fire?"""
import os, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); PAPER = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(PAPER, 'probe_class_rerun')); sys.path.insert(0, HERE)
import pcr_lib as P, sl2_engine as E

n = int(sys.argv[1]) if len(sys.argv) > 1 else 512
nt = int(sys.argv[2]) if len(sys.argv) > 2 else 20
X = P.make_tapes(n, nt, 424242 + n)
Z = E.tape_hashes(X, n) % n
print('n=%d tapes=%d  defect (hash %% n == 0): %d  -> q=%.4f  (1/n=%.4f)'
      % (n, nt, int((Z == 0).sum()), float((Z == 0).mean()), 1.0/n))
for t in range(nt):
    S = P.St(X[t]); R = E.Runner2(S)
    out, steps, msum, inc, cuts, end = R.run('ballast')
    tag = 'DEFECT' if Z[t] == 0 else 'normal'
    print('  t=%-3d %-6s end=%d  out=%7.4f steps=%9d msum=%11d inc=%9d cuts=%5d  inc/n=%7.2f inc/n2=%7.4f'
          % (t, tag, end, out, steps, msum, inc, cuts, inc/n, inc/n**2 if n else 0))
