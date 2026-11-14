# -*- coding: utf-8 -*-
"""mk_diag -- the diagnostics that decide what a one-bit plane could possibly
save: how much of each machine's charge is a re-read, how far apart repeats are,
and what an unbounded certificate store would be worth.

Runs on the same T_n tapes as the sweep.  Writes mk_diag.json.

  pal      reads split b1/b2, comparisons, single-cell reads, seed reads,
           repeat comparisons, distinct (s,i) pairs, lag percentiles, b1 under
           K = 1 and K = inf (the certificate ceiling), scans, clean fraction,
           slot hits, distinct cells touched.
  pattern  inc split (opening / move / head / cut span / cut head), steps, cuts,
           charged cell reads, distinct cells touched, the SCAN column.
"""
import json
import os
import sys
from multiprocessing import Pool

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
PAPER = os.path.dirname(HERE)
for d in ('probe_triple_nat', 'probe_class_rerun', 'probe_e_uniform', 'probe_cert'):
    p = os.path.join(PAPER, d)
    if p not in sys.path:
        sys.path.insert(0, p)

NT_PAL = 300
NT_PAT = 100
NS = (512, 1024, 2048)
P_TABLE = np.array([1, 2, 3, 0], np.int8)      # pcr_lib.TABLE

import paldiag_trace as PT                                      # noqa: E402


class PalCensus(PT.PalTrace):
    """palmach2.PalMach2 (via paldiag_trace) with a census of comparisons,
    cells, repeat lags and scans.  A subclass, so the machine's own `self._scan`
    / `self._cmp` calls land here."""

    def __init__(self, a, K=1):
        PT.PalTrace.__init__(self, a, K=K)
        self.step = 0
        self.lags = []
        self.last = {}
        self.cells = set()
        self.span = []
        self.used = []
        self.clean = []

    def _scan(self, s, i0, i1):
        p, c = PT.PalTrace._scan(self, s, i0, i1)
        self.span.append(i1 - i0 + 1)
        self.used.append(c)
        self.clean.append(1 if p < 0 else 0)
        return p, c

    def _cmp(self, s, i):
        PT.PalTrace._cmp(self, s, i)
        self.cells.add(i)
        self.cells.add(s - i)
        k = (s, i)
        if k in self.last:
            self.lags.append(self.step - self.last[k])
        self.last[k] = self.step

    def _note_cell(self, i):
        PT.PalTrace._note_cell(self, i)
        self.cells.add(i)

    def __call__(self, ell, r):
        self.step += 1
        return PT.PalTrace.__call__(self, ell, r)


def _tape(n, i):
    import tr_lib as T
    for j, a, s, X in T.tapes(n, i + 1):
        if j == i:
            return a, X
    raise IndexError(i)


def pal_worker(job):
    import palmach2
    n, i = job
    a, X = _tape(n, i)
    p = PalCensus(a, K=1)
    d, tot, steps = palmach2.pf.walk(a, p, 'incr')
    o = p
    oi = PT.PalTrace(a, K=10 ** 9)
    palmach2.pf.walk(a, oi, 'incr')
    lag = np.array(p.lags, float)
    return dict(total=tot, b1=o.reads, b2=tot - o.reads, steps=steps,
                cmp=o.cmp, cellreads=o.cellreads, seed=o.seedreads,
                repeats=o.repeats, pairs=len(o.pairs), scans=len(p.span),
                span=float(np.mean(p.span)) if p.span else 0.0,
                per_scan=float(np.mean(p.used)) if p.used else 0.0,
                clean=float(np.mean(p.clean)) if p.clean else 0.0,
                hits=o.hits,
                lag10=float(np.percentile(lag, 10)) if len(lag) else 0.0,
                lag50=float(np.median(lag)) if len(lag) else 0.0,
                lag90=float(np.percentile(lag, 90)) if len(lag) else 0.0,
                cells=len(p.cells), b1_Kinf=oi.reads)


def pat_worker(job):
    import pcr_lib as P
    n, i = job
    a, X = _tape(n, i)
    S = P.St(X)
    d = _pat_trace(S, n)
    out, steps, msum, inc, cuts, end, _ = P.Runner(X).run('pattern')
    assert d['inc'] == inc and d['steps'] == steps, (d['inc'], inc, d['steps'], steps)
    d['scan'] = msum
    d['end'] = end
    return d


def _pat_trace(S, n):
    a, Pa, Ppat = S.a, S.Pa, S.Ppat
    L, R_ = 0, n - 1
    ell, r = L, R_
    y, x = L, R_
    hy, hx = 0, max(1, n // 2) - 1
    steps = 0
    seen = set()
    touched = np.zeros(n, bool)
    d = dict(base=n + 2 * (n - 1) + 2 + max(1, n // 2), move=0, head=0,
             cut_span=0, cut_head=0, cuts=0)
    cap = 40 * n * n + 100

    def touch(u, v):
        if u > v:
            u, v = v, u
        touched[max(0, u):min(n - 1, v) + 1] = True

    while True:
        if ell >= r:
            end = 0
            break
        if R_ - L + 1 < 3:
            end = 1
            break
        steps += 1
        if steps > cap:
            end = 3
            break
        key = ell * (n + 1) + r
        if key in seen:
            end = 2
            break
        seen.add(key)
        m = r - ell + 1
        b1 = 0 if m < 3 else (1 if (Ppat[r - 1] - Ppat[ell]) > 0 else 0)
        k = m >> 1
        if k == 0:
            k = 1
        b2 = 1 if (Pa[r + 1] - Pa[ell]) * k > (Pa[ell + k] - Pa[ell]) * m else 0
        act = P_TABLE[2 * b2 + b1]
        if act >= 4:
            d['cuts'] += 1
            if act == 4:
                L += 1
            elif act == 5:
                R_ -= 1
            else:
                L += 1
                R_ -= 1
            if L > R_:
                end = 1
                break
            ell, r = L, R_
            c = abs(y - L) + abs(x - R_)
            d['cut_span'] += c
            touch(y, L)
            touch(x, R_)
            y, x = L, R_
            k = (x - y + 1) >> 1
            if k == 0:
                k = 1
            ty, tx = y, y + k - 1
            c = abs(hy - ty) + abs(hx - tx)
            d['cut_head'] += c
            touch(hy, ty)
            touch(hx, tx)
            hy, hx = ty, tx
            continue
        dy = 1 if act == 0 else (-1 if act == 1 else 0)
        dx = 1 if act == 2 else (-1 if act == 3 else 0)
        if ell + dy < L or r + dx > R_:
            d['cuts'] += 1
            if ell == L and r == R_:
                L += 1
                R_ -= 1
                y += 1
                x -= 1
                d['cut_span'] += 2
            elif ell + dy < L:
                L += 1
            else:
                R_ -= 1
            if L > R_:
                end = 1
                break
            ell, r = L, R_
            c = abs(y - L) + abs(x - R_)
            d['cut_span'] += c
            touch(y, L)
            touch(x, R_)
            y, x = L, R_
            k = (x - y + 1) >> 1
            if k == 0:
                k = 1
            ty, tx = y, y + k - 1
            c = abs(hy - ty) + abs(hx - tx)
            d['cut_head'] += c
            touch(hy, ty)
            touch(hx, tx)
            hy, hx = ty, tx
            continue
        if ell + dy >= r + dx:
            end = 0
            break
        ell += dy
        r += dx
        y, x = ell, r
        d['move'] += 1
        touch(ell if dy else r, ell if dy else r)
        k = (x - y + 1) >> 1
        if k == 0:
            k = 1
        ty, tx = y, y + k - 1
        c = abs(hy - ty) + abs(hx - tx)
        d['head'] += c
        touch(hy, ty)
        touch(hx, tx)
        hy, hx = ty, tx
    d['steps'] = steps
    d['inc'] = d['base'] + d['move'] + d['head'] + d['cut_span'] + d['cut_head']
    d['charged'] = d['inc'] - d['base']
    d['cells'] = int(touched.sum())
    d['end'] = end
    return d


def agg(rows, extra):
    m = {k: float(np.mean([r[k] for r in rows])) for k in rows[0]}
    m.update(extra)
    return m


def main():
    out = dict(pal=[], pattern=[], nt_pal=NT_PAL, nt_pat=NT_PAT)
    with Pool(min(14, os.cpu_count())) as pool:
        for n in NS:
            rows = pool.map(pal_worker, [(n, i) for i in range(NT_PAL)], chunksize=2)
            out['pal'].append(agg(rows, dict(n=n, tapes=len(rows))))
            print('pal     n=%-5d %s' % (n, {k: round(v, 3) if isinstance(v, float)
                                             else v for k, v in out['pal'][-1].items()}),
                  flush=True)
        for n in NS:
            rows = pool.map(pat_worker, [(n, i) for i in range(NT_PAT)], chunksize=2)
            out['pattern'].append(agg(rows, dict(n=n, tapes=len(rows))))
            print('pattern n=%-5d %s' % (n, {k: round(v, 3) if isinstance(v, float)
                                             else v for k, v in out['pattern'][-1].items()}),
                  flush=True)
    with open(os.path.join(HERE, 'mk_diag.json'), 'w') as f:
        json.dump(out, f)
    print('wrote mk_diag.json')


if __name__ == '__main__':
    main()
