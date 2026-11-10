# -*- coding: utf-8 -*-
"""sl2_content -- does delta depend on the tape, and how much of it?

The trap this probe exists for: a rule whose delta is a coarse summary (a few cells, or
the geometry alone) has e = O(1) and proves nothing.  Four checks, all on the ballast
rule, all measured:

A  SPREAD      delta over the 200-tape cells: sd in cells, quartiles, the fraction of
               tapes with delta in [0.3, 0.7].  A geometry-determined delta has sd 0.
B  PAIR FLIP   pb_content.py's perturbation: one 2-cell -> 1 and one 0-cell -> 1, which
               keeps 2*#2 + #1 = n, i.e. keeps the tape in T_n.  By how many cells does
               delta move?  (the control rules of pb_rule move it by 0-1.)
C  DEFECT FLIP the same pair flip, applied to a DEFECT tape: the hash changes, so the
               branch changes, so delta should move by Theta(n).  Run for every 2-cell of
               one tape, so the count of cells that can move delta is measured, not
               asserted.
D  ERASURE     the direct question behind all three: erase a fraction f of the cells
               (impute 1, the modal value) and re-run.  If a machine that reads only
               (1-f)n cells reproduces delta, then e <= (1-f)n; if it fails on nearly
               every tape, the certificate carries Theta(n) cells and e = Omega(n).

usage: sl2_content.py [n] [ntapes] [npairs]
"""
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
PAPER = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(PAPER, 'probe_class_rerun'))
sys.path.insert(0, HERE)

import pcr_lib as P                                              # noqa: E402
import sl2_engine as E                                           # noqa: E402


def cells(n, nt, seed):
    X = P.make_tapes(n, nt, seed)
    return [P.St(X[t]) for t in range(nt)], X


def walk(S):
    out, steps, msum, inc, cuts, end = E.Runner2(S).run('ballast')
    return out, inc, end


def a_of(X):
    return (X[:, 0] + X[:, 1]).astype(np.int64)


def rebuild(a):
    """a 4-valued tape with this sigma; the label is irrelevant to this rule."""
    Y = np.zeros((len(a), 2), np.int8)
    Y[:, 0] = (a > 0).astype(np.int8)
    Y[:, 1] = (a == 2).astype(np.int8)
    return Y


def pair_flips(a, npairs, rng):
    ins = [i for i, v in enumerate(a) if v == 2 and i != 0]
    outs = [j for j, v in enumerate(a) if v == 0 and j != 0]
    rng.shuffle(ins)
    rng.shuffle(outs)
    return list(zip(ins[:npairs], outs[:npairs]))


def main():
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 512
    nt = int(sys.argv[2]) if len(sys.argv) > 2 else 40
    npairs = int(sys.argv[3]) if len(sys.argv) > 3 else 6
    rng = np.random.default_rng(11)
    SS, X = cells(n, nt, 2026 + n)
    base = []
    for S in SS:
        out, inc, end = walk(S)
        base.append(out if end == 0 else np.nan)
    base = np.array(base)
    print('n=%d  tapes=%d  meeting=%d' % (n, nt, int(np.isfinite(base).sum())))
    ok = np.isfinite(base)
    print('A SPREAD   delta sd %6.1f cells   q25 %.3f  med %.3f  q75 %.3f  '
          '[0.3,0.7] %.2f'
          % (base[ok].std() * n, np.percentile(base[ok], 25), np.median(base[ok]),
             np.percentile(base[ok], 75),
             float(((base[ok] >= 0.3) & (base[ok] <= 0.7)).mean())))

    # ---- B  one-pair flip on the normal tapes
    D = []
    for t in range(nt):
        if not np.isfinite(base[t]):
            continue
        a = a_of(X[t])
        b0 = base[t] * (n - 1)
        for i, j in pair_flips(a, npairs, rng):
            b = a.copy()
            b[i] = 1
            b[j] = 1
            if 2 * int((b == 2).sum()) + int((b == 1).sum()) != n:
                continue
            out, inc, end = walk(P.St(rebuild(b)))
            D.append(n if end != 0 else abs(out * (n - 1) - b0))
    v = np.array(D, float)
    nz = v[v > 0]
    print('B PAIRFLIP %d probes  zero %.1f%%  med %.0f  med(nz) %.0f  p90(nz) %.0f  max %.0f'
          % (len(v), 100.0 * float((v == 0).mean()), np.median(v),
             np.median(nz) if nz.size else 0.0,
             np.percentile(nz, 90) if nz.size else 0.0, v.max()))

    # ---- C  the same flip on a defect tape; every 2-cell gets a partner
    XD, _ = _defect_tape(n)
    aD = a_of(XD)
    outD, incD, endD = walk(P.St(XD))
    zeros = [j for j, v in enumerate(aD) if v == 0 and j != 0]
    twos = [i for i, v in enumerate(aD) if v == 2 and i != 0]
    J = []
    for i in twos:
        j = zeros[rng.integers(len(zeros))]
        b = aD.copy()
        b[i] = 1
        b[j] = 1
        out, inc, end = walk(P.St(rebuild(b)))
        J.append((i, n if end != 0 else abs(out * (n - 1) - outD * (n - 1)),
                  inc))
    J = np.array(J, float)
    print('C DEFECTFLIP base delta %.3f (end %d, %.3f n^2 reads);  %d 2-cells probed: '
          'delta moved by med %.0f cells, p10 %.0f, max %.0f; %d/%d moved > n/4'
          % (outD, endD, incD / n ** 2, len(J), np.median(J[:, 1]),
             np.percentile(J[:, 1], 10), J[:, 1].max(),
             int((J[:, 1] > n / 4).sum()), len(J)))

    # ---- D  erasure: read (1-f) of the cells, impute 1 elsewhere
    for f in (0.001, 0.01, 0.05, 0.25, 0.5):
        hit = 0
        tot = 0
        dd = []
        for t in range(nt):
            if not np.isfinite(base[t]):
                continue
            a = a_of(X[t])
            k = max(1, int(round(f * n)))
            idx = rng.choice(n, k, replace=False)
            b = a.copy()
            b[idx] = 1
            out, inc, end = walk(P.St(rebuild(b)))
            tot += 1
            if end == 0 and abs(out - base[t]) < 0.5 / (n - 1):
                hit += 1
            dd.append(n if end != 0 else abs(out * (n - 1) - base[t] * (n - 1)))
        print('D ERASE f=%.3f (%d cells)  delta reproduced %d/%d  med move %.0f cells'
              % (f, max(1, int(round(f * n))), hit, tot, float(np.median(dd))))


def _defect_tape(n):
    done = 0
    while done < 200000:
        X = P.make_tapes(n, 2000, 90210 + done)
        Z = E.tape_hashes(X, n) % n
        idx = np.nonzero(Z == 0)[0]
        if idx.size:
            return X[idx[0]], done + int(idx[0])
        done += 2000
    raise SystemExit('no defect tape found')


if __name__ == '__main__':
    main()
