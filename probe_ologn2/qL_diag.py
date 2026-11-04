# -*- coding: utf-8 -*-
"""qL_diag.py -- epoch-level diagnostics for Lemma Q.

For each tape and each epoch i we record (L,R,M,m*,cost) plus the quantities
that decide the NEXT epoch, so that we can study

  (a) the conditional law of cost_i given F_{i-1}, proxied by history features
      (run length of preceding double cuts, index i, top-contrast value, ...)
  (b) how much the top test moves from one epoch to the next
  (c) the one-step hazard of the *top* test as a function of the run length
  (d) block sums  sum_{i in block} cost_i  vs  sum_{i in block} M_i

Run:  $P qL_diag.py <n> <ntapes> <seed>
"""
import sys
import numpy as np
from q0_walk import prefix, epoch, sample_tapes


def trace(a):
    """Full trace with per-epoch contrast values.  Returns list of dicts."""
    n = len(a)
    P = prefix(a)
    L, R = 0, n - 1
    out = []
    while True:
        if L >= R:
            break
        M = R - L + 1
        if M < 3:
            break
        mstar, cost = epoch(P, L, R)
        rec = dict(L=L, R=R, M=M, mstar=mstar, cost=cost)
        # top contrast: even m closest to M
        mtop = M - (M % 2)
        k = mtop // 2
        rec['mtop'] = mtop
        rec['ktop'] = k
        rec['Dtop'] = int(2 * (P[L + k] - P[L]) - (P[L + 2 * k] - P[L]))  # 2S_k - S_2k = -D_k
        rec['Dtop'] = -rec['Dtop']                                        # D_k = S_2k - 2S_k
        out.append(rec)
        if mstar is None:
            break
        if mstar == M:
            L, R = L + 1, R - 1
        else:
            L += 1
    return out


def main(n, nt, seed):
    rows = []           # per epoch rows
    taperows = []       # per tape
    rng = np.random.default_rng(seed + 1)
    for ti, a in enumerate(sample_tapes(n, nt, seed)):
        tr = trace(a)
        if not tr:
            continue
        T = sum(r['cost'] for r in tr)
        SM = sum(r['M'] for r in tr)
        taperows.append((T, SM, len(tr)))
        run = 0
        for i, r in enumerate(tr):
            prevD = tr[i - 1]['Dtop'] if i > 0 else np.nan
            prevM = tr[i - 1]['M'] if i > 0 else np.nan
            rows.append((i / n, r['M'] / n, r['cost'] / max(r['M'], 1),
                         r['cost'], r['M'], run,
                         r['Dtop'] / np.sqrt(max(r['ktop'], 1)),
                         prevD / np.sqrt(max(tr[i - 1]['ktop'], 1)) if i > 0 else np.nan,
                         T / n ** 2, SM / n ** 2))
            # run of double cuts: a double cut is cost==1 (mstar==M)
            if r['mstar'] is not None and r['mstar'] == r['M']:
                run += 1
            else:
                run = 0
    R = np.array(rows, float)
    (ii, Mn, cM, cost, M, run, dtop, dprev, Tn2, SMn2) = R.T
    T = np.array([t for t, _, _ in taperows], float)
    SM = np.array([s for _, s, _ in taperows], float)
    NE = np.array([e for _, _, e in taperows], float)
    print('n=%d ntapes=%d  E[T]/n^2=%.4f  E[SM]/n^2=%.4f  E[epochs]/n=%.3f  E[cost/M]=%.4f'
          % (n, len(T), T.mean() / n ** 2, SM.mean() / n ** 2, NE.mean() / n,
             cost.sum() / M.sum()))

    print('\n(1) mean cost/M by run length of preceding double cuts')
    for k in range(0, 12):
        sel = (run == k)
        if sel.sum() > 200:
            print('   run=%2d  epochs %7d  mean cost/M %.4f  P(cost>=.1M) %.3f  P(cost>=2) %.3f'
                  % (k, sel.sum(), cM[sel].mean(), np.mean(cost[sel] >= 0.1 * M[sel]),
                     np.mean(cost[sel] >= 2)))
    sel = (run >= 12)
    if sel.sum() > 200:
        print('   run>=12  epochs %7d  mean cost/M %.4f  P(cost>=.1M) %.3f  P(cost>=2) %.3f'
              % (sel.sum(), cM[sel].mean(), np.mean(cost[sel] >= 0.1 * M[sel]),
                 np.mean(cost[sel] >= 2)))

    print('\n(2) mean cost/M by D_top (normalised top contrast at this epoch)')
    for lo, hi in [(-9, -2), (-2, -1), (-1, -0.5), (-0.5, 0), (0, 0.5), (0.5, 1), (1, 2), (2, 9)]:
        sel = (dtop >= lo) & (dtop < hi)
        if sel.sum() > 200:
            print('   Dtop/sqrt(k) in [%5.1f,%5.1f)  epochs %7d  mean cost/M %.4f  P(cost>=.1M) %.3f'
                  % (lo, hi, sel.sum(), cM[sel].mean(), np.mean(cost[sel] >= 0.1 * M[sel])))

    print('\n(3) mean cost/M by epoch index quartile of the window size')
    for lo, hi in [(.5, .6), (.6, .7), (.7, .8), (.8, .9), (.9, 1.01)]:
        sel = (Mn >= lo) & (Mn < hi)
        if sel.sum() > 200:
            print('   M/n in [%.1f,%.1f)  epochs %7d  mean cost/M %.4f  P(cost=1) %.3f'
                  % (lo, hi, sel.sum(), cM[sel].mean(), np.mean(cost[sel] == 1)))

    print('\n(4) block sums: for blocks of b epochs, fraction of blocks with')
    print('    sum cost >= 0.01 * sum M  and  mean')
    for b in (5, 10, 20, 40):
        nb = len(T) // b
        if nb < 1:
            continue
        cs = cost[:nb * b].reshape(nb, b).sum(1)
        ms = M[:nb * b].reshape(nb, b).sum(1)
        ratio = cs / ms
        print('    b=%2d  blocks %5d  mean ratio %.4f  P(ratio>=.01) %.3f  P(ratio>=.05) %.3f  P(ratio>=.1) %.3f'
              % (b, nb, ratio.mean(), np.mean(ratio >= .01), np.mean(ratio >= .05), np.mean(ratio >= .1)))


if __name__ == '__main__':
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 512,
         int(sys.argv[2]) if len(sys.argv) > 2 else 200,
         int(sys.argv[3]) if len(sys.argv) > 3 else 2026)
