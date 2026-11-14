# -*- coding: utf-8 -*-
"""mk_probe -- diagnostics on the two baselines before any plane is added.

Questions:
  pal      how many reads are b1 (comparisons) and how many b2 (the S/Q carry)?
           how much does certificate memory buy?  K = 1 (the machine), K = 4,
           K = inf (a dict -- out of class, the ceiling any tape-resident
           certificate store could reach).  And: how many comparisons repeat a
           pair (s, i) already compared?  A repeat is the only thing a per-cell
           memo can ever save.
  pattern  decomposition of `inc` into opening / per-move / head / cut, and the
           fraction of charged cell reads that touch a cell already read.
"""
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
PAPER = os.path.dirname(HERE)
for d in ('probe_e_uniform', 'probe_class_rerun', 'probe_cert', 'probe_triple_nat'):
    p = os.path.join(PAPER, d)
    if p not in sys.path:
        sys.path.insert(0, p)
sys.path.insert(0, HERE)

import tr_lib as T                                              # noqa: E402
import paldiag_trace as PT                                      # noqa: E402


def pal_split(a, K=1):
    """b1 reads and (total - b1) = b2 reads for palmach2 at slot count K."""
    o = PT.PalTrace(a, K=K)
    d, tot, steps = PT.walk(a, o)
    return dict(total=tot, b1=o.reads, b2=tot - o.reads, steps=steps,
                cmp=o.cmp, pairs=len(o.pairs), repeats=o.repeats,
                cellreads=o.cellreads, cellrep=o.cellrep, seed=o.seedreads)


def pattern_split(n, nt=3):
    import pcr_lib as P
    import tr_lib as T
    rows = []
    for i, a, s, X in T.tapes(n, nt):
        r = P.Runner(X).run('pattern')
        out, steps, msum, inc, cuts, end, _ = r
        rows.append(dict(steps=steps, inc=inc, cuts=cuts, end=end,
                         base=n + 2 * (n - 1) + 2 + max(1, n // 2)))
    return rows


if __name__ == '__main__':
    print('pal: total / b1 (comparisons) / b2 / steps / repeat-comparisons, K=1')
    for n in (512, 1024):
        acc = []
        for i, a, s, X in T.tapes(n, 12):
            acc.append(pal_split(a, K=1))
        m = {k: float(np.mean([r[k] for r in acc])) for k in acc[0]}
        print('  n=%-5d total=%-9.1f b1=%-9.1f b2=%-8.1f steps=%-8.1f '
              'cmp=%-8.1f repeats=%-8.1f cellreads=%-8.1f cellrep=%-6.1f seed=%.1f'
              % (n, m['total'], m['b1'], m['b2'], m['steps'], m['cmp'],
                 m['repeats'], m['cellreads'], m['cellrep'], m['seed']))
    print('pal: b1 reads under K = 1 / 4 / inf (certificate memory ceiling)')
    for n in (512, 1024):
        for K in (1, 4, 10 ** 9):
            acc, tot = [], []
            for i, a, s, X in T.tapes(n, 8):
                o = PT.PalTrace(a, K=K)
                acc.append(o)
                tot.append(PT.walk(a, o)[1])
            b1 = [o.reads for o in acc]
            print('  n=%-5d K=%-10s b1=%-9.1f total=%-9.1f'
                  % (n, ('inf' if K > 10 ** 8 else K), np.mean(b1), np.mean(tot)))
    print('pattern: inc decomposition')
    for n in (512, 1024):
        rows = pattern_split(n, 3)
        for r in rows[:3]:
            print('  n=%-5d steps=%-8d inc=%-9d base=%-6d rest=%-9d cuts=%-5d end=%d'
                  % (n, r['steps'], r['inc'], r['base'], r['inc'] - r['base'],
                     r['cuts'], r['end']))
