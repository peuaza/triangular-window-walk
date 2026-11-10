# -*- coding: utf-8 -*-
"""q2: exact b1=0 trajectory statistics by exhaustive enumeration, n = 6..13.

Writes q2_exact.json with, per n:
  tapes, mean_T, T/n^2, P(T=n-1), P(T<=2n), P(T<=n^2/10), max_T,
  mean epochs, histogram of T (first 40 values), double-cut fraction.
"""
import json
import sys
import numpy as np
sys.path.insert(0, '../probe_ologn')
from p0_tapes import enum_tapes
from q0_walk import walk_b0_fast


def main(ns):
    out = {}
    for n in ns:
        Ts, Ls, dc, ep, cut = [], [], 0, 0, 0
        hist = {}
        for a in enum_tapes(n):
            d, T, diag = walk_b0_fast(a, want_diag=True)
            Ts.append(T)
            ep += len(diag)
            for (L, R, M, ms, cost) in diag:
                if ms == M:
                    dc += 1
                cut += 1
            if d is not None:
                Ls.append(int(round(d * (n - 1))))
            hist[T] = hist.get(T, 0) + 1
        Ts = np.array(Ts, dtype=float)
        Ls = np.array(Ls)
        rec = dict(tapes=len(Ts), mean_T=float(Ts.mean()), T_over_n2=float(Ts.mean() / n ** 2),
                   P_T_eq_nm1=float(np.mean(Ts == n - 1)),
                   P_T_le_2n=float(np.mean(Ts <= 2 * n)),
                   P_T_le_n2_10=float(np.mean(Ts <= n * n / 10)),
                   max_T=float(Ts.max()), mean_epochs=ep / len(Ts),
                   double_cut_frac=dc / cut,
                   Lstar_mean=float(Ls.mean()) / (n - 1),
                   Lstar_hist=[int(x) for x in np.bincount(Ls, minlength=n)],
                   T_hist={int(k): int(v) for k, v in sorted(hist.items())})
        out[n] = rec
        print('n=%2d tapes=%6d  E[T]=%9.1f  E[T]/n^2=%.4f  P(T=n-1)=%.4f  P(T<=2n)=%.4f  '
              'P(T<=n^2/10)=%.3f  maxT=%d  epochs=%.2f  double-cut=%.3f'
              % (n, rec['tapes'], rec['mean_T'], rec['T_over_n2'], rec['P_T_eq_nm1'],
                 rec['P_T_le_2n'], rec['P_T_le_n2_10'], rec['max_T'], rec['mean_epochs'],
                 rec['double_cut_frac']))
    with open('q2_exact.json', 'w') as f:
        json.dump(out, f, indent=1)


if __name__ == '__main__':
    main([int(x) for x in sys.argv[1:]] or [6, 7, 8, 9, 10, 11, 12, 13])
