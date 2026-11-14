# -*- coding: utf-8 -*-
"""mk_probe2 -- where does pal's 35% repeat-comparison rate come from, and how
long are the scans a tape-resident certificate would have to beat?"""
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
PAPER = os.path.dirname(HERE)
for d in ('probe_triple_nat', 'probe_cert', 'probe_e_uniform'):
    p = os.path.join(PAPER, d)
    if p not in sys.path:
        sys.path.insert(0, p)
sys.path.insert(0, HERE)

import tr_lib as T                                              # noqa: E402
import palfast as pf                                            # noqa: E402
import paldiag_trace as PT                                      # noqa: E402


class PalDeep(PT.PalTrace):
    """adds scan-length and lag bookkeeping"""

    def __init__(self, a, K=1):
        PT.PalTrace.__init__(self, a, K=K)
        self.scanlens = []
        self.laststep = {}
        self.lags = []
        self.step = 0
        self.hitcnt = 0

    def _cmp(self, s, i):
        PT.PalTrace._cmp(self, s, i)
        k = (s, i)
        if k in self.laststep:
            self.lags.append(self.step - self.laststep[k])
        self.laststep[k] = self.step

    def _scan(self, s, i0, i1):
        p, c = PT.PalTrace._scan(self, s, i0, i1)
        self.scanlens.append((i1 - i0 + 1, c, p < 0))
        return p, c

    def __call__(self, ell, r):
        self.step += 1
        return PT.PalTrace.__call__(self, ell, r)


def one(a, K):
    o = PalDeep(a, K=K)
    d, tot, steps = pf.walk(a, o)
    lens = np.array([x[0] for x in o.scanlens], float)
    used = np.array([x[1] for x in o.scanlens], float)
    clean = np.array([x[2] for x in o.scanlens])
    lags = np.array(o.lags, float)
    return dict(total=tot, b1=o.reads, scans=len(o.scanlens),
                meanlen=lens.mean() if len(lens) else 0,
                maxlen=int(lens.max()) if len(lens) else 0,
                meanused=used.mean() if len(used) else 0,
                clean=float(clean.mean()) if len(clean) else 0,
                cmp=o.cmp, repeats=o.repeats,
                lag_med=float(np.median(lags)) if len(lags) else -1,
                lag_lo=float(np.percentile(lags, 10)) if len(lags) else -1,
                lag_hi=float(np.percentile(lags, 90)) if len(lags) else -1,
                longcmp=int((used >= 8).sum()),
                longreads=float(used[used >= 8].sum()),
                steps=steps)


if __name__ == '__main__':
    for n in (512, 1024):
        for K in (1, 4, 10 ** 9):
            acc = []
            for i, a, s, X in T.tapes(n, 10):
                acc.append(one(a, K))
            m = {k: float(np.mean([r[k] for r in acc])) for k in acc[0]}
            print('n=%-5d K=%-4s total=%-9.0f b1=%-8.0f scans=%-7.0f '
                  'len=%-5.1f max=%-5.0f used=%-4.1f clean=%-4.2f cmp=%-7.0f '
                  'rep=%-6.0f lag50=%-5.0f lag10=%-5.0f lag90=%-6.0f '
                  'longscans=%-5.0f longreads=%-7.0f'
                  % (n, ('inf' if K > 10 ** 8 else K), m['total'], m['b1'],
                     m['scans'], m['meanlen'], m['maxlen'], m['meanused'],
                     m['clean'], m['cmp'], m['repeats'], m['lag_med'],
                     m['lag_lo'], m['lag_hi'], m['longcmp'], m['longreads']))
