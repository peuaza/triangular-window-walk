# -*- coding: utf-8 -*-
"""ce_lib -- CES exactness probe.

Ground truth: kp_lib.exact (the reference walk) with its per-epoch trace.
An epoch record is (M0, b2, b1, oc, ln, L, R, ee):
    M0 = region width at epoch start, (b2,b1) the two phi bits at the start,
    oc in {'D','R','L'} for a cut epoch, or 'meet'/'narrow'/'cycle' at the
    terminal epoch, ln = number of walk steps in the epoch.

CES's assumption: oc = 'D' if b2 == 0 else ('R' if b1 else 'L').
Two failure modes on the true trajectory:
  cut_mismatch  the epoch cuts, but on the other side than predicted;
  early_term    the epoch produces no cut at all (meet/cycle/narrow) while CES
                would have taken a cut.  Early terminals at M0 <= W are handled
                by CES's exact endgame buffer; only M0 > W ones are failures.
"""
import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
for _p in (os.path.join(ROOT, 'probe_kappa_pred'),
           os.path.join(ROOT, 'probe_ologn2')):
    if _p not in sys.path:
        sys.path.insert(0, _p)

import kp_lib                                                   # noqa: E402


def pred_oc(b2, b1):
    if not b2:
        return 'D'
    return 'R' if b1 else 'L'


def trace(a, rule):
    """Exact walk + per-epoch trace; returns (delta, why, eps)."""
    ex = kp_lib.exact(a, rule)
    return ex['delta'], ex['why'], ex['eps']


def cut_epochs(eps):
    return [e for e in eps if e[3] in ('D', 'R', 'L')]


def measure_tape(a, rule, W=16):
    """Per-tape failure summary of the CES outcome rule, all epochs."""
    d, why, eps = trace(a, rule)
    recs = []
    first_fail = None
    nfail = ncut_mismatch = nterm = 0
    for j, e in enumerate(eps):
        M0, b2, b1, oc, ln, L, R, ee = e
        p = pred_oc(b2, b1)
        terminal = oc not in ('D', 'R', 'L')
        ok = (p == oc)
        # a terminal epoch is a failure of the shortcut iff it happens above W
        fail = (not ok) and (M0 > W)
        nfail += fail
        if not ok and terminal:
            nterm += 1
        if not ok and not terminal:
            ncut_mismatch += 1
        if fail and first_fail is None:
            first_fail = dict(idx=j, M0=int(M0), b2=int(b2), b1=int(b1),
                              pred=p, true=oc, ln=int(ln), L=int(L), R=int(R),
                              kind='term' if terminal else 'cut')
        recs.append(dict(j=j, M0=int(M0), b2=int(b2), b1=int(b1), pred=p,
                         oc=oc, ln=int(ln), ok=bool(ok), fail=bool(fail),
                         terminal=bool(terminal)))
    n_short = sum(1 for r in recs if r['M0'] > W)
    nfail_short = sum(1 for r in recs if r['fail'])
    return dict(delta=None if d is None or d != d else float(d), why=why,
                neps=len(recs), nfail=nfail, first_fail=first_fail,
                nfail_short=nfail_short, n_short=n_short,
                ncut_mismatch=ncut_mismatch, nterm=nterm,
                epochs=recs)


def cell(tapes, rule, W=16):
    out = [measure_tape(a, rule, W) for a in tapes]
    return out
