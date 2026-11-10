# -*- coding: utf-8 -*-
"""sl2_gate -- the new walk must be pcr_lib's walk, tape by tape.

The canonical rule exists in both engines.  For every tape and both cut conventions
the two must agree EXACTLY on: output, step count, m-sum, INC reads, cut count,
ending.  Any mismatch is printed and counted; a non-zero count voids the probe.
"""
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
PAPER = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(PAPER, 'probe_class_rerun'))
sys.path.insert(0, HERE)

import pcr_lib as P                                              # noqa: E402
import sl2_engine as E                                           # noqa: E402

FIELDS = ('out', 'steps', 'msum', 'inc', 'cuts', 'end')


def main():
    ns = [int(v) for v in (sys.argv[1].split(',') if len(sys.argv) > 1
                           else ['128', '512', '1024'])]
    nt = int(sys.argv[2]) if len(sys.argv) > 2 else 40
    bad = 0
    checked = 0
    for n in ns:
        X = P.make_tapes(n, nt, 90001 + n)
        for t in range(nt):
            S = P.St(X[t])
            R1 = P.Runner(X[t])
            for conv in (0, 1):
                ref = R1.run('canonical', conv=conv)[:6]
                got = E.Runner2(S).run('canonical', conv=conv)
                checked += 1
                for k in range(6):
                    a, b = ref[k], got[k]
                    if isinstance(a, float) or isinstance(b, float):
                        ok = (a == b) or (abs(a - b) < 1e-12)
                    else:
                        ok = int(a) == int(b)
                    if not ok:
                        bad += 1
                        print('MISMATCH n=%d tape=%d conv=%d %s: %r vs %r'
                              % (n, t, conv, FIELDS[k], a, b), flush=True)
        print('n=%-5d %d tapes x 2 conv   mismatch so far %d/%d'
              % (n, nt, bad, checked), flush=True)
    print('GATE %s  (%d runs checked, %d field mismatches)'
          % ('PASS' if bad == 0 else 'FAIL', checked, bad))
    return 0 if bad == 0 else 1


if __name__ == '__main__':
    sys.exit(main())
