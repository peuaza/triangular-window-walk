# -*- coding: utf-8 -*-
"""ce_check -- validate the machine against the reference walk and against
kp_lib.ces.  Three gates:

  A  machine(W=3, W2=n)  == kp_lib.exact        (delta and nan)
  B  machine(W=W0)       == kp_lib.ces(W0)      (delta and nan)
  C  machine(W=W0, K=4)  == kp_lib.ces(W0)      where the epoch is a cut
     (prefix check must not change a cut outcome)
"""
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
for _p in (HERE, os.path.join(ROOT, 'probe_e_uniform'),
           os.path.join(ROOT, 'probe_kappa_pred'),
           os.path.join(ROOT, 'probe_ologn2')):
    if _p not in sys.path:
        sys.path.insert(0, _p)

import ce_lib as CL                                              # noqa: E402
import ce_machine as CM                                          # noqa: E402
import eu_lib as EL                                              # noqa: E402
import kp_lib as KP                                              # noqa: E402


def same(x, y):
    xn = x is None or x != x
    yn = y is None or y != y
    if xn or yn:
        return xn and yn
    return abs(x - y) < 1e-12


def main():
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 128
    nt = int(sys.argv[2]) if len(sys.argv) > 2 else 20
    bad = dict(A=0, B=0, Bval=0, C=0)
    for i in range(nt):
        a, s, X = EL.make_tape(n, 5000 + n + 100003 * i)
        a = np.asarray(a, np.int64)
        for rule in ('pattern', 'runs'):
            W0 = 14 if rule == 'pattern' else max(16, 4 * int(round(np.log2(n))))
            d, why, eps = CL.trace(a, rule)
            ref = None if (d is None or d != d) else float(d)
            # A: full exact machine
            rA = CM.machine(a, rule, W=3, W2=n)[0]
            if not same(ref, rA):
                bad['A'] += 1
                print('A mismatch i=%d %s ref=%s got=%s why=%s' % (i, rule, ref, rA, why))
            # B: pure CES
            dB, _, _ = KP.ces(a, rule, W=W0)
            rB, rd, inf = CM.machine(a, rule, W=W0)
            if not same(dB, rB):
                bad['B'] += 1
                print('B mismatch i=%d %s kpces=%s mach=%s end=%s' % (i, rule, dB, rB, inf))
            # C: prefix-4 CES (must equal B whenever no terminal is caught)
            rC, _, infC = CM.machine(a, rule, W=W0, K=4)
            if not same(rB, rC):
                bad['C'] += 1
                print('C differs i=%d %s ces=%s K4=%s inf=%s' % (i, rule, rB, rC, infC))
    print('n=%d nt=%d mismatches: A(fullexact)=%d B(ces)=%d C(K4)=%d'
          % (n, nt, bad['A'], bad['B'], bad['C']))


if __name__ == '__main__':
    main()
