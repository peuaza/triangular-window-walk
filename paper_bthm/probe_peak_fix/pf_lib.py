# -*- coding: utf-8 -*-
"""pf_lib -- the peak machine of EN_LINEAR plus a literal first epoch.

The base machine is `el_cell.jm_walk` (probe_en_linear, the code that produced
EN_LINEAR.md): state = the region [L,R] + one sliding counter; while M > W2 the
epoch's cut is read off (b2,b1); while M <= W2 the epoch is run literally.

`jm_fixed` is that machine with the first epoch replaced by a literal walk of
the epoch (the only epoch EN_LINEAR's failure analysis blames).  The literal
epoch is charged its worst-case price, n reads, once, so the machine's read
bound is n (opening) + n (first epoch) + n/2 (cuts) + endgame.

Both functions return the machine's own (delta, end) and a charge; neither
consults the truth.  `audit` (optional) compares the machine's epoch sequence
with the trace's while the two have not yet diverged, which is how the first
wrong epoch is located.
"""
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
PAPER = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(PAPER, 'probe_en_linear'))

import el_cell as C
import el_lib as E
import pcr_lib as P

TAB = P.TABLE_NP
CONV = dict(P.CONV)
CONV.setdefault('BAL>m/4', 0)


def b1_at(S, rule, ell, r):
    return C.b1_at(S, rule, ell, r)


def b2_at(S, ell, r):
    return C.b2_at(S, ell, r)


def _literal_epoch(S, rule, n, L, R, seen, audit, step0, stop_at_cut=True):
    """Run one epoch literally, exactly as pcr_lib.walk does.

    Returns (kind, L, R, steps) with kind in
      0 meet, 1 region, 2 cycle, 3 cut.
    `audit` gets (M0, b2, b1, act) per step while it is not None.
    """
    ell, r = L, R
    ln = 0
    while True:
        if ell >= r:
            return 0, L, R, ln
        if R - L + 1 < 3:
            return 1, L, R, ln
        if (ell, r) in seen:
            return 2, L, R, ln
        seen.add((ell, r))
        m = r - ell + 1
        b1 = b1_at(S, rule, ell, r)
        b2 = b2_at(S, ell, r)
        ln += 1
        if audit is not None:
            audit.append((m, int(b2), int(b1)))
        act = TAB[2 * b2 + b1]
        dy = 1 if act == P.LPLUS else (-1 if act == P.LMINUS else 0)
        dx = 1 if act == P.RPLUS else (-1 if act == P.RMINUS else 0)
        if ell + dy < L or r + dx > R:
            if CONV.get(rule, 0) == 1:
                if ell + dy < L:
                    R -= 1
                else:
                    L += 1
            else:
                if ell == L and r == R:
                    L += 1
                    R -= 1
                elif ell + dy < L:
                    L += 1
                else:
                    R -= 1
            if L > R:
                return 1, L, R, ln
            return 3, L, R, ln
        if ell + dy >= r + dx:
            return 0, L, R, ln
        if r + dx < L or ell + dy > R:
            L += 1
            R -= 1
            if L > R:
                return 1, L, R, ln
            ell, r = L, R
            seen.clear()
            if stop_at_cut:
                return 3, L, R, ln
            continue
        ell += dy
        r += dx


def jm_fixed(rule, S, n, W2, mode='sliding', first_charge=None, audit=False):
    """The jump machine with the first epoch run literally.

    Returns dict(delta, end, reads, epochs, cuts, eg, first_len,
                 first_steps, nscan, ep_audit).
    `reads` charges n (opening) + first_charge (the literal first epoch, n by
    default) + 1 per cell a cut gives up + endgame steps (1 per step, +m per
    step for a scanning statistic).
    `first_steps` is the measured length of the literal first epoch.
    """
    if first_charge is None:
        first_charge = n
    sv = S.s
    L, R = 0, n - 1
    reads = n + first_charge
    epochs = 0
    cuts = 0
    eg = 0
    nscan = 0
    first_len = 0
    kind = None
    first_audit = [] if audit else None
    Ssum = 0
    for j in range(n):
        Ssum += int(sv[j])

    # ---- the literal first epoch (only when the region is above the buffer)
    if n > W2:
        seen = set()
        kind, L, R, first_len = _literal_epoch(
            S, rule, n, L, R, seen, first_audit, 0)
        epochs += 1
        if kind in (0, 2):
            d = (L / (n - 1.0)) if kind == 0 else None
            return dict(delta=d, end=(0 if kind == 0 else 2), reads=reads,
                        epochs=epochs, cuts=cuts, eg=eg, first_len=first_len,
                        first_steps=first_len, nscan=nscan,
                        ep_audit=first_audit or [])
        if kind == 1:
            return dict(delta=None, end=1, reads=reads, epochs=epochs,
                        cuts=cuts, eg=eg, first_len=first_len,
                        first_steps=first_len, nscan=nscan,
                        ep_audit=first_audit or [])
        # kind == 3: the epoch ended with a cut, L,R are the new region
        if L > R:
            return dict(delta=None, end=1, reads=reads, epochs=epochs,
                        cuts=cuts, eg=eg, first_len=first_len,
                        first_steps=first_len, nscan=nscan,
                        ep_audit=first_audit or [])
        # charge the cut the epoch gave up
        cuts += 1
        reads += 2 if (n - (R - L + 1)) == 2 else 1

    ep_audit = [] if audit else None
    lazy = (mode == 'lazy')
    ismax = rule in ('MAXP', 'AMAX')
    p = -1
    pv = 0
    if lazy:
        qhi = (R + 1) if rule == 'MAXP' else R
        acc = 0
        best = 0
        p = L
        for j in range(L, qhi + 1):
            if j > L:
                acc += int(sv[j - 1])
            if (acc > best) if ismax else (acc < best):
                best = acc
                p = j
        pv = best

    while True:
        M = R - L + 1
        if M < 3:
            return dict(delta=None, end=1, reads=reads, epochs=epochs,
                        cuts=cuts, eg=eg, first_len=first_len,
                        first_steps=first_len, nscan=nscan,
                        ep_audit=ep_audit or [])
        if M <= W2:
            ell, r = L, R
            seen = set()
            while True:
                if ell >= r:
                    return dict(delta=ell / (n - 1.0), end=0, reads=reads,
                                epochs=epochs, cuts=cuts, eg=eg,
                                first_len=first_len, first_steps=first_len,
                                nscan=nscan, ep_audit=ep_audit or [])
                if R - L + 1 < 3:
                    return dict(delta=None, end=1, reads=reads, epochs=epochs,
                                cuts=cuts, eg=eg, first_len=first_len,
                                first_steps=first_len, nscan=nscan,
                                ep_audit=ep_audit or [])
                if (ell, r) in seen:
                    return dict(delta=None, end=2, reads=reads, epochs=epochs,
                                cuts=cuts, eg=eg, first_len=first_len,
                                first_steps=first_len, nscan=nscan,
                                ep_audit=ep_audit or [])
                seen.add((ell, r))
                m = r - ell + 1
                b1 = b1_at(S, rule, ell, r)
                b2 = b2_at(S, ell, r)
                reads += 1 + (0 if mode == 'sliding' else m)
                eg += 1
                act = TAB[2 * b2 + b1]
                dy = 1 if act == P.LPLUS else (-1 if act == P.LMINUS else 0)
                dx = 1 if act == P.RPLUS else (-1 if act == P.RMINUS else 0)
                if ell + dy < L or r + dx > R:
                    if CONV.get(rule, 0) == 1:
                        if ell + dy < L:
                            R -= 1
                        else:
                            L += 1
                    else:
                        if ell == L and r == R:
                            L += 1
                            R -= 1
                        elif ell + dy < L:
                            L += 1
                        else:
                            R -= 1
                    if L > R:
                        return dict(delta=None, end=1, reads=reads,
                                    epochs=epochs, cuts=cuts, eg=eg,
                                    first_len=first_len,
                                    first_steps=first_len, nscan=nscan,
                                    ep_audit=ep_audit or [])
                    ell, r = L, R
                    seen.clear()
                    continue
                if ell + dy >= r + dx:
                    return dict(delta=(ell + dy) / (n - 1.0), end=0,
                                reads=reads, epochs=epochs, cuts=cuts, eg=eg,
                                first_len=first_len, first_steps=first_len,
                                nscan=nscan, ep_audit=ep_audit or [])
                if r + dx < L or ell + dy > R:
                    L += 1
                    R -= 1
                    if L > R:
                        return dict(delta=None, end=1, reads=reads,
                                    epochs=epochs, cuts=cuts, eg=eg,
                                    first_len=first_len,
                                    first_steps=first_len, nscan=nscan,
                                    ep_audit=ep_audit or [])
                    ell, r = L, R
                    seen.clear()
                    continue
                ell += dy
                r += dx
        else:
            if lazy:
                qhi = (R + 1) if rule == 'MAXP' else R
                if p < L or p > qhi:
                    reads += qhi - L + 1
                    nscan += 1
                    acc = 0
                    best = 0
                    p = L
                    for j in range(L, qhi + 1):
                        if j > L:
                            acc += int(sv[j - 1])
                        if (acc > best) if ismax else (acc < best):
                            best = acc
                            p = j
                    pv = best
                if rule == 'MAXP':
                    b1 = 1 if pv * pv > M else 0
                elif rule == 'AMAX':
                    b1 = 1 if Ssum > pv else 0
                else:
                    b1 = 1 if Ssum < pv else 0
            else:
                b1 = b1_at(S, rule, L, R)
                if mode == 'scan':
                    reads += M
            b2 = b2_at(S, L, R)
            epochs += 1
            pr = C.pred_of(rule, b2, b1)
            if ep_audit is not None:
                ep_audit.append((M, int(b2), int(b1), pr))
            cuts += 1
            if pr == 3:
                reads += 2
                dw = int(sv[L])
                Ssum -= dw + int(sv[R])
                if lazy:
                    pv -= dw
                L += 1
                R -= 1
            elif pr == 2:
                reads += 1
                Ssum -= int(sv[R])
                R -= 1
            else:
                reads += 1
                dw = int(sv[L])
                Ssum -= dw
                if lazy:
                    pv -= dw
                L += 1


def jm_base_audit(rule, S, n, W2, mode='sliding'):
    """el_cell.jm_walk, copied, with the machine's epoch audit attached."""
    a = S.a
    sv = S.s
    L, R = 0, n - 1
    reads = n
    epochs = 0
    cuts = 0
    eg = 0
    nscan = 0
    audit = []
    Ssum = 0
    for j in range(n):
        Ssum += int(sv[j])
    while True:
        M = R - L + 1
        if M < 3:
            return dict(delta=None, end=1, reads=reads, epochs=epochs,
                        cuts=cuts, eg=eg, nscan=nscan, audit=audit)
        if M <= W2:
            ell, r = L, R
            seen = set()
            while True:
                if ell >= r:
                    return dict(delta=ell / (n - 1.0), end=0, reads=reads,
                                epochs=epochs, cuts=cuts, eg=eg,
                                nscan=nscan, audit=audit)
                if R - L + 1 < 3:
                    return dict(delta=None, end=1, reads=reads, epochs=epochs,
                                cuts=cuts, eg=eg, nscan=nscan, audit=audit)
                if (ell, r) in seen:
                    return dict(delta=None, end=2, reads=reads, epochs=epochs,
                                cuts=cuts, eg=eg, nscan=nscan, audit=audit)
                seen.add((ell, r))
                m = r - ell + 1
                b1 = b1_at(S, rule, ell, r)
                b2 = b2_at(S, ell, r)
                reads += 1 + (0 if mode == 'sliding' else m)
                eg += 1
                act = TAB[2 * b2 + b1]
                dy = 1 if act == P.LPLUS else (-1 if act == P.LMINUS else 0)
                dx = 1 if act == P.RPLUS else (-1 if act == P.RMINUS else 0)
                if ell + dy < L or r + dx > R:
                    if CONV.get(rule, 0) == 1:
                        if ell + dy < L:
                            R -= 1
                        else:
                            L += 1
                    else:
                        if ell == L and r == R:
                            L += 1
                            R -= 1
                        elif ell + dy < L:
                            L += 1
                        else:
                            R -= 1
                    if L > R:
                        return dict(delta=None, end=1, reads=reads,
                                    epochs=epochs, cuts=cuts, eg=eg,
                                    nscan=nscan, audit=audit)
                    ell, r = L, R
                    seen.clear()
                    continue
                if ell + dy >= r + dx:
                    return dict(delta=(ell + dy) / (n - 1.0), end=0,
                                reads=reads, epochs=epochs, cuts=cuts,
                                eg=eg, nscan=nscan, audit=audit)
                if r + dx < L or ell + dy > R:
                    L += 1
                    R -= 1
                    if L > R:
                        return dict(delta=None, end=1, reads=reads,
                                    epochs=epochs, cuts=cuts, eg=eg,
                                    nscan=nscan, audit=audit)
                    ell, r = L, R
                    seen.clear()
                    continue
                ell += dy
                r += dx
        else:
            b1 = b1_at(S, rule, L, R)
            if mode == 'scan':
                reads += M
            b2 = b2_at(S, L, R)
            epochs += 1
            pr = C.pred_of(rule, b2, b1)
            audit.append((M, int(b2), int(b1), pr))
            cuts += 1
            if pr == 3:
                reads += 2
                Ssum -= int(sv[L]) + int(sv[R])
                L += 1
                R -= 1
            elif pr == 2:
                reads += 1
                Ssum -= int(sv[R])
                R -= 1
            else:
                reads += 1
                Ssum -= int(sv[L])
                L += 1


def exact_of(m, t, n):
    """The EN_LINEAR exactness test: machine (delta,end) against the trace."""
    d, en = m['delta'], m['end']
    return (int(en) == int(t['end'])) and (
        (d is None and t['delta'] < 0) or
        (d is not None and t['delta'] >= 0 and abs(d - t['delta']) < 1e-12))


def first_divergence(m, t):
    """Index and M0 of the first epoch where the machine's cut differs.

    Valid only while the machine has not diverged: its region at epoch k is the
    trace's, so its predicted `out` can be compared with epout[k].
    """
    aud = m.get('audit')
    if aud is None:
        return None
    for k, (M, b2, b1, pr) in enumerate(aud):
        if k >= t['nep']:
            return (k, M, pr, -1)
        if int(t['epout'][k]) != pr:
            return (k, M, pr, int(t['epout'][k]))
    return None
