# -*- coding: utf-8 -*-
"""q9: first-epoch cost E[cost_0(n)] (the rigorous ingredient) and one big-n value."""
import sys
import numpy as np
from q0_walk import prefix, epoch, walk_b0_fast, sample_tapes


def first_epoch(ns=(64, 128, 256, 512, 1024, 2048, 4096, 8192), ntapes=2000):
    print('%6s %8s %10s %10s %8s' % ('n', 'tapes', 'E[cost_0]', '/n', 'P(meet)'))
    for n in ns:
        cs = []
        meet = 0
        for a in sample_tapes(n, ntapes if n <= 2048 else 600, 77):
            P = prefix(a)
            ms, c = epoch(P, 0, n - 1)
            cs.append(c)
            meet += (ms is None)
        cs = np.array(cs, float)
        print('%6d %8d %10.1f %10.4f %8.4f  se/n %.4f'
              % (n, len(cs), cs.mean(), cs.mean() / n, meet / len(cs),
                 cs.std(ddof=1) / np.sqrt(len(cs)) / n))
        sys.stdout.flush()


def big(n=8192, ntapes=400):
    Ts, ep = [], []
    for a in sample_tapes(n, ntapes, 31337):
        d, T, diag = walk_b0_fast(a, want_diag=True)
        Ts.append(T)
        ep.append(len(diag))
    Ts = np.array(Ts, float)
    print('\nn=%d tapes=%d  E[T]=%.0f  E[T]/n^2=%.4f +- %.4f  median %.4f  '
          'P(T>=0.1n^2)=%.3f  P(T<=n)=%.3f  epochs/n=%.3f'
          % (n, ntapes, Ts.mean(), Ts.mean() / n ** 2,
             Ts.std(ddof=1) / np.sqrt(ntapes) / n ** 2, np.median(Ts) / n ** 2,
             np.mean(Ts >= 0.1 * n ** 2), np.mean(Ts <= n), np.mean(ep) / n))


if __name__ == '__main__':
    first_epoch()
    big()
