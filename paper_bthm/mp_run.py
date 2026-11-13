"""mp -- the machine profile: for each rule, (G, e, g) of several machines, so the three minima
G_rho = min_M exp E[ln C], e_rho = min_M E[C], g_rho = min_M max C can be read off the same table.

Machines implemented here (all O(log n) bits of state, all charged INC):
  pcr            pcr_lib.Runner, the canonical engine
  bandB<B>       hmach3.H3 with buffer B in {16,32,64,128}   (band rules only)
  cesW<W>        kp_lib.ces with window W in {8,16,32,64}    (the rules it supports)
  scan           the same walk re-evaluating its window: steps x m reads (SCAN charge)

Usage:
  python mp_run.py cell RULE N [K]
  python mp_run.py grid --ns 512,1024 --rules a,b [--jobs 1]
Writes probe_minima/mp_cells.json, one record per (rule, machine, n) with the per-tape reads.
"""
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
PAPER = os.path.dirname(HERE)
for d in ('probe_class_rerun', 'probe_kappa_pred', 'probe_triple_nat', 'probe_hybrid3',
          'probe_en_linear', 'probe_peak_fix'):
    sys.path.insert(0, os.path.join(PAPER, d))

import pcr_lib as P            # noqa: E402
import tr_lib as T             # noqa: E402
import el_lib as E             # noqa: E402
try:
    import hmach3             # noqa: E402
except Exception:
    hmach3 = None
try:
    import kp_lib as K         # noqa: E402
except Exception:
    K = None

OUT = os.path.join(HERE, 'mp_cells.json')
BAND_TAG = {'NBLK=1': 'nblk_1', 'NBLK>2': 'nblk_gt2', 'NBLK>3': 'nblk_gt3',
            'AMAX': 'amax_end', 'AMIN': 'amin_end'}
BANDS = (16, 32, 64, 128)
CESW = (8, 16, 32, 64)
CES_RULES = ('pattern', 'runs', 'pal')


def tape(i, n):
    """The ensemble tape used by every probe: two independent h-subsets, h = ceil(n/2)."""
    X = E.tape(n, i)
    s = (np.asarray(X)[:, 0] - np.asarray(X)[:, 1]).astype(np.int64)
    a = (np.asarray(X)[:, 0] + np.asarray(X)[:, 1]).astype(np.int64)
    return X, s, a


def _i(x):
    """int, or None when the engine reports nothing -- never raise on a missing field."""
    return None if x is None else int(x)


def _ok(e):
    """meet in any of the three encodings the engines use (0, 'meet', 'M1b')."""
    return e == 0 or (isinstance(e, str) and e in ('meet', 'M1b'))


def run_one(rule, mname, X, s, a):
    n = len(a)
    if mname == 'pcr':
        r = T.run(rule, a, s)
        return _i(r.get('reads')), 0 if _ok(r.get('end', 1)) else 1
    if mname.startswith('band'):
        B = int(mname[4:])
        r = hmach3.H3([int(v) for v in s], [int(v) for v in a], BAND_TAG[rule],
                      B=B, log=False, deep=False).run()
        return _i(r.get('reads')), 0 if _ok(r.get('end')) else 1
    if mname.startswith('ces'):
        W = int(mname[3:])
        d, reads, info = K.ces(a, rule=rule, W=W)
        return _i(reads), 0 if np.isfinite(d) else 1
    if mname == 'scan':
        r = T.run(rule, a, s)
        steps = r.get('steps')
        if steps is None:            # this engine does not report steps; the scan charge is undefined
            return None, 1
        return int(steps) * n, 0 if _ok(r.get('end', 1)) else 1
    raise KeyError(mname)


LEAN = True          # after the first pass: one buffer, no scan, K=100 -- see mp_run.py header


def machines_for(rule):
    if LEAN:
        ms = ['pcr']
        if rule in BAND_TAG:
            ms.append('band64')
        if rule in CES_RULES:
            ms += ['ces64']
        return ms
    ms = ['pcr', 'scan']
    if rule in BAND_TAG:
        ms += ['band%d' % B for B in BANDS]
    if rule in CES_RULES:
        ms += ['ces%d' % W for W in CESW]
    return ms


def cell(rule, n, K):
    out = json.load(open(OUT)) if os.path.exists(OUT) else []
    have = {(r['rule'], r['machine'], r['n']) for r in out}
    ms = machines_for(rule)
    for mn in ms:
        if (rule, mn, n) in have:
            continue
        reads, ends = [], []
        i = 0
        while len(reads) < K and i < 40 * K:
            X, s, a = tape(i, n)
            try:
                r, e = run_one(rule, mn, X, s, a)
            except Exception as ex:                    # one machine must not kill the grid
                print('%-10s %-8s n=%5d SKIP %s' % (rule, mn, n, ex), flush=True)
                reads = []
                break
            i += 1
            if e == 0 and r is not None:
                reads.append(r)
        rec = dict(rule=rule, machine=mn, n=n, K=len(reads), draws=i,
                   reads=reads, capped=len(reads) < K)
        out.append(rec)
        json.dump(out, open(OUT, 'w'))
        print('%-10s %-8s n=%5d K=%3d mean=%.4g max=%.4g' %
              (rule, mn, n, len(reads), np.mean(reads) if reads else 0,
               max(reads) if reads else 0), flush=True)
    return out


def main():
    if sys.argv[1] == 'cell':
        cell(sys.argv[2], int(sys.argv[3]), int(sys.argv[4]) if len(sys.argv) > 4 else 200)
    elif sys.argv[1] == 'grid':
        ns = [int(v) for v in sys.argv[sys.argv.index('--ns') + 1].split(',')]
        rules = sys.argv[sys.argv.index('--rules') + 1].split(',')
        K = int(sys.argv[sys.argv.index('--K') + 1]) if '--K' in sys.argv else 200
        for rule in rules:
            for n in ns:
                cell(rule, n, K)


if __name__ == '__main__':
    main()
