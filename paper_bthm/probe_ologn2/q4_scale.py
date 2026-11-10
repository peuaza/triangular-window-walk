# -*- coding: utf-8 -*-
"""q4: scale-up of E[T]/n^2 (does the constant drift down?) + structural profile.

Part A (scale): for each n, sample tapes, report mean T/n^2 with its standard
error, median, P(T>=0.1 n^2), P(T<=n), and the composition of the short tapes
(#double cuts / #epochs).
Part B (profile): for n = 1024, average the 8-block mean profile of the tape
over tapes in T/n^2 bands, to see what structure makes a long trajectory.
"""
import json
import sys
import numpy as np
from q0_walk import walk_b0_fast, sample_tapes


def scale(ns, ntapes, seed=2026):
    out = {}
    for n in ns:
        Ts, dc, ep, Ls, ncut1 = [], [], [], [], []
        for a in sample_tapes(n, ntapes, seed):
            d, T, diag = walk_b0_fast(a, want_diag=True)
            Ts.append(T)
            ep.append(len(diag))
            dc.append(sum(1 for (_, _, M, ms, c) in diag if ms == M))
            ncut1.append(sum(1 for (_, _, M, ms, c) in diag if c == 1))
            Ls.append(int(round(d * (n - 1))) if d is not None else -1)
        Ts = np.array(Ts, dtype=float)
        epa = np.array(ep, dtype=float)
        dca = np.array(dc, dtype=float)
        Lsa = np.array(Ls, dtype=float)
        nca = np.array(ncut1, dtype=float)
        short = Ts <= n
        rec = dict(n=n, ntapes=ntapes, mean=float(Ts.mean()),
                   se=float(Ts.std(ddof=1) / np.sqrt(len(Ts))),
                   mean_over_n2=float(Ts.mean() / n ** 2),
                   se_over_n2=float(Ts.std(ddof=1) / np.sqrt(len(Ts)) / n ** 2),
                   median_over_n2=float(np.median(Ts) / n ** 2),
                   P_ge_010=float(np.mean(Ts >= 0.10 * n ** 2)),
                   P_ge_005=float(np.mean(Ts >= 0.05 * n ** 2)),
                   P_le_n=float(np.mean(short)),
                   mean_epochs=float(epa.mean()), mean_epochs_over_n=float(epa.mean() / n),
                   dc_over_ep=float((dca / np.maximum(epa, 1)).mean()),
                   Lstar_over_n=float((Lsa[Lsa >= 0].mean()) / (n - 1)),
                   short_frac_cost1=float((nca[short] / np.maximum(epa[short], 1)).mean())
                   if short.any() else None,
                   long_frac_cost1=float((nca[~short] / np.maximum(epa[~short], 1)).mean())
                   if (~short).any() else None)
        out[n] = rec
        print('n=%5d N=%4d  E[T]/n^2 = %.4f +- %.4f   med %.4f   P(T>=.1n^2)=%.3f  '
              'P(T<=n)=%.3f  epochs/n=%.3f  dc/ep=%.2f  L*/n=%.2f  cost1frac short/long %.2f/%.2f'
              % (n, ntapes, rec['mean_over_n2'], rec['se_over_n2'], rec['median_over_n2'],
                 rec['P_ge_010'], rec['P_le_n'], rec['mean_epochs_over_n'],
                 rec['dc_over_ep'], rec['Lstar_over_n'],
                 rec['short_frac_cost1'] if rec['short_frac_cost1'] is not None else -1,
                 rec['long_frac_cost1'] if rec['long_frac_cost1'] is not None else -1))
        sys.stdout.flush()
    with open('q4_scale.json', 'w') as f:
        json.dump(out, f, indent=1)


def profile(n=1024, ntapes=600, seed=99):
    bands = [(0, 0.005), (0.005, 0.02), (0.02, 0.06), (0.06, 0.12), (0.12, 2)]
    B = 8
    acc = {b: np.zeros(B) for b in range(len(bands))}
    cnt = {b: 0 for b in range(len(bands))}
    for a in sample_tapes(n, ntapes, seed):
        d, T, _ = walk_b0_fast(a)
        r = T / n ** 2
        prof = a[: (n // B) * B].reshape(B, -1).mean(1)
        for i, (lo, hi) in enumerate(bands):
            if lo <= r < hi:
                acc[i] += prof
                cnt[i] += 1
                break
    print('\nq4B: 8-block mean profile by T/n^2 band (n=%d, %d tapes)' % (n, ntapes))
    for i, (lo, hi) in enumerate(bands):
        if cnt[i]:
            print('  T/n^2 in [%.3f,%.3f): %4d tapes   blocks '
                  % (lo, hi, cnt[i]) + ' '.join('%.2f' % x for x in acc[i] / cnt[i]))
        else:
            print('  T/n^2 in [%.3f,%.3f): 0 tapes' % (lo, hi))


if __name__ == '__main__':
    if len(sys.argv) > 1 and sys.argv[1] == 'prof':
        profile()
    else:
        ns = [int(x) for x in sys.argv[1:]] or [512, 1024, 2048, 4096]
        scale(ns, 400)
