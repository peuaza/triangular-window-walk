# -*- coding: utf-8 -*-
"""q11: the half-swap identity and the conditional law given the sign of the drift.

Global drift:  Delta(a) = mean(a[n/2..n-1]) - mean(a[0..n/2-1])   (n even).
The half swap sigma (exchange the two halves) is a permutation of positions, so it
preserves the ensemble, and Delta(sigma a) = -Delta(a).  Hence EXACTLY

    E[T] = E[ (T(a) + T(sigma a)) * 1{Delta(a) >= 0} ]                (*)

so the expectation is carried by the "positive drift" half alone.  This script
verifies (*) and measures E[T/n^2 | Delta >= 0] versus n.
"""
import sys
import numpy as np
from q0_walk import walk_b0_fast, sample_tapes


def main(ns, ntapes):
    for n in ns:
        assert n % 2 == 0
        h = n // 2
        T, D, Tsig = [], [], []
        for a in sample_tapes(n, ntapes, 5150):
            _, t, _ = walk_b0_fast(a)
            s = np.concatenate((a[h:], a[:h]))
            _, ts, _ = walk_b0_fast(s)
            T.append(t)
            Tsig.append(ts)
            D.append(a[h:].mean() - a[:h].mean())
        T = np.array(T, float)
        Tsig = np.array(Tsig, float)
        D = np.array(D, float)
        pos = D >= 0
        rhs = ((T + Tsig) * pos).mean()
        print('n=%5d N=%4d  P(Delta>=0)=%.3f  E[T]/n^2=%.4f   E[T/n^2|D>=0]=%.4f  '
              'E[T/n^2|D<0]=%.4f  identity(*): lhs %.1f vs rhs %.1f (ratio %.4f)'
              % (n, ntapes, pos.mean(), T.mean() / n ** 2, T[pos].mean() / n ** 2,
                 T[~pos].mean() / n ** 2, T.mean(), rhs, T.mean() / rhs))
        # sanity: same for the swapped tape (should reproduce the same numbers)
        print('        T(sigma a): mean/n^2 %.4f  |D>=0 %.4f  |D<0 %.4f   '
              'P(T>=0.01n^2 | D>=0)=%.3f'
              % (Tsig.mean() / n ** 2, Tsig[pos].mean() / n ** 2,
                 Tsig[~pos].mean() / n ** 2, np.mean(T[pos] >= 0.01 * n ** 2)))
        sys.stdout.flush()


if __name__ == '__main__':
    main([int(x) for x in sys.argv[1:]] or [128, 256, 512, 1024, 2048], 600)
