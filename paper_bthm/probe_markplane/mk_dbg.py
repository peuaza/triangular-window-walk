# -*- coding: utf-8 -*-
"""mk_dbg -- why does a `cert` tape diverge?  Compare the oracle answer streams."""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
PAPER = os.path.dirname(HERE)
for d in ('probe_triple_nat', 'probe_cert', 'probe_e_uniform'):
    p = os.path.join(PAPER, d)
    if p not in sys.path:
        sys.path.insert(0, p)

import tr_lib as T                                              # noqa: E402
import palfast as pf                                            # noqa: E402
import paldiag_trace as PT                                      # noqa: E402
import mk_pal                                                   # noqa: E402


class Log:
    def __init__(self, inner):
        self.inner = inner
        self.log = []

    @property
    def reads(self):
        return self.inner.reads

    def __call__(self, ell, r):
        v = self.inner(ell, r)
        self.log.append((ell, r, v))
        return v


def main(n=512, i=3, pol='cert'):
    tp = list(T.tapes(n, max(i + 1, 6)))
    _, a, s, X = tp[i]
    o1 = PT.PalTrace(a, K=1)
    L1 = Log(o1)
    d1, t1, st1 = pf.walk(a, L1)
    o2 = mk_pal.PalMarked(a, policy=pol)
    L2 = Log(o2)
    d2, t2, st2 = pf.walk(a, L2)
    print('base  reads=%-8d delta=%s end?=%s steps=%d' % (t1, d1, st1, len(L1.log)))
    print('%s reads=%-8d delta=%s steps=%d consults=%d'
          % (pol, t2, d2, len(L2.log), o2.consults))
    for k in range(min(len(L1.log), len(L2.log))):
        if L1.log[k][2] != L2.log[k][2]:
            print('first differing query at', k, L1.log[k], L2.log[k])
            prev = L1.log[max(0, k - 6):k]
            print('previous queries:', prev)
            break


if __name__ == '__main__':
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 512,
         int(sys.argv[2]) if len(sys.argv) > 2 else 3,
         sys.argv[3] if len(sys.argv) > 3 else 'cert')
