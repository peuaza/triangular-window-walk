# -*- coding: utf-8 -*-
"""sl2_measure -- e, g, G and the defect fraction for the `ballast` rule.

PROTOCOL
  200 meeting tapes per cell (n = 512, 1024, 2048, 4096), pcr_lib's T_n
  (two independent uniform h-subsets, h = n//2), seed 2026 + n, tape i = draw i.
  A 4000-tape sample per cell is run as well, because a 1/n defect cannot be seen on
  200 tapes at n = 4096; the two estimates of e are printed side by side so the
  mixture estimator used for the headline is checked, not asserted.
  The defect witnesses are the first three tapes of the search with hash % n == 0.
  The defect fraction is measured on 200 000 hash-only draws per cell (spin up a
  meeting walk only for the 200 + 4000 + 3 tapes).

COSTS, per tape, from the INC accounting (sl2_engine):
  inc = n + 2(n-1) + 2 + n//2  (the window-aggregate build, cr_lib's Win.__init__)
        + 1 per cell taken in or given up + the head drift
        + the re-span at every cut (the reset forces the window back onto the region
          and the O(log n)-bit state cannot store it) + 2 for a double cut.
  msum (SCAN accounting, m reads per query) is stored too, not used for e.

ESTIMATORS
  raw    mean / max / geo-mean of the 200 (and 4000) uniform meeting tapes
  mix    e = (1-q) E[C|normal] + q E[C|defect], with q measured on the 200k draws and
         E[C|defect] on the three witnesses; G = exp((1-q) E[ln C|normal] + q ln C_def)
"""
import json
import math
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
PAPER = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(PAPER, 'probe_class_rerun'))
sys.path.insert(0, HERE)

import pcr_lib as P                                              # noqa: E402
import sl2_engine as E                                           # noqa: E402

NS = tuple(int(v) for v in (os.environ.get('SL2_NS', '512,1024,2048,4096').split(',')))
NT = int(os.environ.get('SL2_NT', 200))
NT_BIG = int(os.environ.get('SL2_NTBIG', 4000))
NDEF = 3
Q_N = int(os.environ.get('SL2_QN', 200000))
OUT = os.path.join(HERE, 'sl2_measure.json')


class LightSt(object):
    """a, Pa, Paa, PR only -- everything the ballast rule's predicates read."""

    def __init__(self, X):
        X = np.asarray(X)
        n = self.n = int(X.shape[0])
        a = (X[:, 0] + X[:, 1]).astype(np.int8)
        self.a = a
        ai = a.astype(np.int64)
        Pa = np.zeros(n + 1, np.int64)
        Pa[1:] = np.cumsum(ai)
        Paa = np.zeros(n + 1, np.int64)
        Paa[1:] = np.cumsum(ai * ai)
        PR = np.zeros(n + 1, np.int64)
        if n > 1:
            PR[1:n] = np.cumsum(ai[:-1] * ai[1:])
            PR[n] = PR[n - 1]
        self.Pa, self.Paa, self.PR = Pa, Paa, PR


def run_one(S, stamp_n):
    R = E.Runner2(S)
    out, steps, msum, inc, cuts, end = R.run('ballast')
    return dict(out=float(out), steps=int(steps), msum=int(msum), inc=int(inc),
                cuts=int(cuts), end=int(end))


def defect_fraction(n, qn, seed):
    """q = P(hash % n == 0), hash-only, in batches (no tape stays resident)."""
    hit = 0
    done = 0
    batch = 500
    while done < qn:
        b = min(batch, qn - done)
        X = P.make_tapes(n, b, seed + done)
        Z = E.tape_hashes(X, n) % n
        hit += int((Z == 0).sum())
        done += b
    return hit, done


def find_defect(n, want, seed):
    """first `want` tapes with hash % n == 0, in draw order."""
    out = []
    done = 0
    batch = 2000
    while len(out) < want and done < 60 * n + 20000:
        b = min(batch, 60 * n + 20000 - done)
        X = P.make_tapes(n, b, seed + done)
        Z = E.tape_hashes(X, n) % n
        for t in np.nonzero(Z == 0)[0]:
            out.append((done + int(t), X[t].copy()))
            if len(out) >= want:
                break
        done += b
    return out, done


def geom(v):
    v = np.asarray(v, float)
    return float(np.exp(np.log(v).mean()))


def fit(xs, ys):
    lx, ly = np.log(np.asarray(xs, float)), np.log(np.asarray(ys, float))
    b, a = np.polyfit(lx, ly, 1)
    return float(b)


def main():
    t0 = time.time()
    cells = []
    for n in NS:
        rec = dict(n=n)
        # --- the 200-tape protocol cell
        X = P.make_tapes(n, NT, 2026 + n)
        rows = [run_one(P.St(X[t]), n) for t in range(NT)]
        meet = [r for r in rows if r['end'] == 0]
        rec['nt'] = NT
        rec['n_meet'] = len(meet)
        rec['reads'] = [r['inc'] for r in meet]
        rec['steps'] = [r['steps'] for r in meet]
        rec['msum'] = [r['msum'] for r in meet]
        rec['delta'] = [r['out'] for r in meet]
        rec['e_raw'] = float(np.mean([r['inc'] for r in meet]))
        rec['g_raw'] = float(np.max([r['inc'] for r in meet]))
        rec['G_raw'] = geom([r['inc'] for r in meet])
        print('n=%-5d 200-tape cell: meet %d/%d  inc/n mean %.3f  max %.3f  G/n %.3f  '
              'steps/n %.3f  delta sd(cells) %.1f'
              % (n, len(meet), NT, rec['e_raw'] / n, rec['g_raw'] / n,
                 rec['G_raw'] / n, np.mean(rec['steps']) / n,
                 float(np.std(rec['delta'])) * n), flush=True)
        # --- the 4000-tape tail check (light struct)
        XL = P.make_tapes(n, NT_BIG, 777 + n)
        Z = E.tape_hashes(XL, n) % n
        inc_big = np.zeros(NT_BIG, np.int64)
        end_big = np.zeros(NT_BIG, np.int8)
        defmask = (Z == 0)
        for t in range(NT_BIG):
            r = run_one(LightSt(XL[t]), n)
            inc_big[t] = r['inc']
            end_big[t] = r['end']
        mm = end_big == 0
        rec['nt_big'] = NT_BIG
        rec['n_meet_big'] = int(mm.sum())
        rec['n_defect_big'] = int((defmask & mm).sum())
        rec['e_big'] = float(inc_big[mm].mean())
        rec['g_big'] = float(inc_big[mm].max())
        rec['G_big'] = geom(inc_big[mm])
        rec['inc_normal_big'] = float(inc_big[mm & ~defmask].mean()) \
            if (mm & ~defmask).any() else float('nan')
        print('  4000-tape check: meet %d  defects among them %d  e/n %.3f  max/n^2 %.4f'
              % (int(mm.sum()), rec['n_defect_big'], rec['e_big'] / n,
                 rec['g_big'] / n ** 2), flush=True)
        # --- defect witnesses
        dw, tried = find_defect(n, NDEF, 31337 + n)
        wr = []
        for t, Xt in dw:
            r = run_one(P.St(Xt), n)
            wr.append(r)
            print('  defect witness draw %-6d end=%d out=%.4f inc=%.4f n^2  cuts=%d'
                  % (t, r['end'], r['out'], r['inc'] / n ** 2, r['cuts']), flush=True)
        rec['witness'] = wr
        rec['defect_tried'] = tried
        # --- defect fraction
        hit, done = defect_fraction(n, Q_N, 515151 + n)
        rec['q_hit'], rec['q_n'] = hit, done
        rec['q'] = hit / float(done)
        rec['q_ci'] = 1.96 * math.sqrt(max(rec['q'] * (1 - rec['q']), 1e-12) / done)
        print('  defect fraction %d/%d = %.6f +- %.6f   (1/n = %.6f)'
              % (hit, done, rec['q'], rec['q_ci'], 1.0 / n), flush=True)
        # --- mixture
        cd = float(np.mean([r['inc'] for r in wr]))
        cn = rec['inc_normal_big']
        q = rec['q']
        rec['c_defect'] = cd
        rec['c_normal'] = cn
        rec['e_mix'] = (1 - q) * cn + q * cd
        rec['g_mix'] = max(rec['g_big'], max(r['inc'] for r in wr))
        lnc = (1 - q) * float(np.log(inc_big[mm & ~defmask]).mean()) + q * math.log(cd)
        rec['G_mix'] = float(math.exp(lnc))
        rec['e_over_n'] = rec['e_mix'] / n
        rec['g_over_n2'] = rec['g_mix'] / n ** 2
        rec['G_over_n'] = rec['G_mix'] / n
        print('  MIX e/n %.3f  g/n^2 %.4f  G/n %.3f  (C_norm/n %.3f, C_def/n^2 %.4f)'
              % (rec['e_over_n'], rec['g_over_n2'], rec['G_over_n'],
                 cn / n, cd / n ** 2), flush=True)
        cells.append(rec)
        print('  [%.1fs]' % (time.time() - t0), flush=True)

    ns = [c['n'] for c in cells]
    ex = dict(
        e_mix=fit(ns, [c['e_mix'] for c in cells]),
        g_mix=fit(ns, [c['g_mix'] for c in cells]),
        G_mix=fit(ns, [c['G_mix'] for c in cells]),
        e_raw=fit(ns, [c['e_raw'] for c in cells]),
        g_raw=fit(ns, [c['g_raw'] for c in cells]),
        c_normal=fit(ns, [c['c_normal'] for c in cells]),
        c_defect=fit(ns, [c['c_defect'] for c in cells]),
        q=fit(ns, [c['q'] for c in cells]),
        one_over_n=fit(ns, [1.0 / n for n in ns]),
    )
    json.dump(dict(cells=cells, ex=ex), open(OUT, 'w'), indent=1)
    print('exponents:', {k: round(v, 3) for k, v in ex.items()})
    print('wrote', OUT)


if __name__ == '__main__':
    main()
