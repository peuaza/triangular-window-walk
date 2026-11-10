# -*- coding: utf-8 -*-
"""pf_var -- machines in which the literal/jump choice is a per-epoch decision.

`jm_sel(rule,S,n,W2,sel,pred)` evaluates (b2,b1) at each epoch of the region and
asks `sel(M, b2, b1, k)` whether to run that epoch literally.  A literal epoch
costs its steps (1 read per step, the record's endgame charge); a jumped epoch
costs 1 read per cell the cut gives up.  The opening pass is n reads.

The point is to price the candidates that could repair EN_LINEAR's peak
failures without giving up Theta(n):
  * first epoch only (`sel = k == 0`), charged n once;
  * the leading run of b2=1,b1=0 scans (`sel = b2==1 and b1==0 and still leading`);
  * a flipped top-half predictor (no literal epochs; pred override).
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import pf_lib as F
import el_cell as C
import el_lib as E
import pcr_lib as P

TAB = P.TABLE_NP
CONV = F.CONV


def jm_sel(rule, S, n, W2, sel, pred=None, first_charge=None, audit=False):
    """The jump machine with a per-epoch literal decision.

    sel(M, b2, b1, k) -> True to run the epoch literally.
    pred(M, b2, b1)   -> the predicted outcome (default el_cell.pred_of).
    first_charge      -> if not None, epoch 0 is charged this many reads instead
                         of its measured length (the "cost at most n, once"
                         convention).
    """
    if pred is None:
        def pred(M, b2, b1):
            return C.pred_of(rule, b2, b1)
    sv = S.s
    L, R = 0, n - 1
    reads = n
    epochs = 0
    cuts = 0
    eg = 0
    lit_steps = 0
    lit_epochs = 0
    audit_ep = [] if audit else None
    while True:
        M = R - L + 1
        if M < 3:
            return dict(delta=None, end=1, reads=reads, epochs=epochs,
                        cuts=cuts, eg=eg, lit_steps=lit_steps,
                        lit_epochs=lit_epochs, audit=audit_ep)
        if M <= W2:
            # ---- the endgame: the epoch is run literally, as in jm_walk
            ell, r = L, R
            seen = set()
            while True:
                if ell >= r:
                    return dict(delta=ell / (n - 1.0), end=0, reads=reads,
                                epochs=epochs, cuts=cuts, eg=eg,
                                lit_steps=lit_steps, lit_epochs=lit_epochs,
                                audit=audit_ep)
                if R - L + 1 < 3:
                    return dict(delta=None, end=1, reads=reads, epochs=epochs,
                                cuts=cuts, eg=eg, lit_steps=lit_steps,
                                lit_epochs=lit_epochs, audit=audit_ep)
                if (ell, r) in seen:
                    return dict(delta=None, end=2, reads=reads, epochs=epochs,
                                cuts=cuts, eg=eg, lit_steps=lit_steps,
                                lit_epochs=lit_epochs, audit=audit_ep)
                seen.add((ell, r))
                m = r - ell + 1
                bb1 = F.b1_at(S, rule, ell, r)
                bb2 = F.b2_at(S, ell, r)
                reads += 1
                eg += 1
                act = TAB[2 * bb2 + bb1]
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
                                    lit_steps=lit_steps,
                                    lit_epochs=lit_epochs, audit=audit_ep)
                    ell, r = L, R
                    seen.clear()
                    continue
                if ell + dy >= r + dx:
                    return dict(delta=(ell + dy) / (n - 1.0), end=0,
                                reads=reads, epochs=epochs, cuts=cuts, eg=eg,
                                lit_steps=lit_steps, lit_epochs=lit_epochs,
                                audit=audit_ep)
                if r + dx < L or ell + dy > R:
                    L += 1
                    R -= 1
                    if L > R:
                        return dict(delta=None, end=1, reads=reads,
                                    epochs=epochs, cuts=cuts, eg=eg,
                                    lit_steps=lit_steps,
                                    lit_epochs=lit_epochs, audit=audit_ep)
                    ell, r = L, R
                    seen.clear()
                    continue
                ell += dy
                r += dx

        b1 = F.b1_at(S, rule, L, R)
        b2 = F.b2_at(S, L, R)
        do_lit = bool(sel(M, int(b2), int(b1), epochs))
        if audit_ep is not None:
            audit_ep.append((M, int(b2), int(b1), 1 if do_lit else 0))
        if do_lit:
            ln = 0
            seen = set()
            ell, r = L, R
            kind = None
            meetpos = -1
            while True:
                if ell >= r:
                    kind = 0
                    meetpos = ell
                    break
                if R - L + 1 < 3:
                    kind = 1
                    break
                if (ell, r) in seen:
                    kind = 2
                    break
                seen.add((ell, r))
                m = r - ell + 1
                bb1 = F.b1_at(S, rule, ell, r)
                bb2 = F.b2_at(S, ell, r)
                reads += 1
                ln += 1
                act = TAB[2 * bb2 + bb1]
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
                    kind = 3
                    break
                if ell + dy >= r + dx:
                    kind = 0
                    meetpos = ell + dy
                    break
                if r + dx < L or ell + dy > R:
                    L += 1
                    R -= 1
                    kind = 3
                    break
                ell += dy
                r += dx
            epochs += 1
            lit_epochs += 1
            if first_charge is not None and epochs == 1:
                reads += first_charge - ln   # replace the measured steps by the cap
                lit_steps += first_charge
            else:
                lit_steps += ln
            if kind == 0:
                return dict(delta=meetpos / (n - 1.0),
                            end=0, reads=reads, epochs=epochs, cuts=cuts,
                            eg=eg, lit_steps=lit_steps, lit_epochs=lit_epochs,
                            audit=audit_ep)
            if kind in (1, 2):
                return dict(delta=None, end=(1 if kind == 1 else 2),
                            reads=reads, epochs=epochs, cuts=cuts, eg=eg,
                            lit_steps=lit_steps, lit_epochs=lit_epochs,
                            audit=audit_ep)
            if L > R:
                return dict(delta=None, end=1, reads=reads, epochs=epochs,
                            cuts=cuts, eg=eg, lit_steps=lit_steps,
                            lit_epochs=lit_epochs, audit=audit_ep)
            cuts += 1
            reads += 2 if (M - (R - L + 1)) == 2 else 1
            continue
        # ---- jump epoch
        epochs += 1
        pr = pred(M, int(b2), int(b1))
        cuts += 1
        if pr == 3:
            reads += 2
            L += 1
            R -= 1
        elif pr == 2:
            reads += 1
            R -= 1
        else:
            reads += 1
            L += 1
        if L > R:
            return dict(delta=None, end=1, reads=reads, epochs=epochs,
                        cuts=cuts, eg=eg, lit_steps=lit_steps,
                        lit_epochs=lit_epochs, audit=audit_ep)


def _lit_epoch(S, rule, n, L, R):
    """One literal epoch; returns (kind, L, R, steps, meetpos).

    kind: 0 meet, 1 region, 2 cycle, 3 cut.
    """
    ell, r = L, R
    ln = 0
    seen = set()
    meetpos = -1
    while True:
        if ell >= r:
            return 0, L, R, ln, ell
        if R - L + 1 < 3:
            return 1, L, R, ln, -1
        if (ell, r) in seen:
            return 2, L, R, ln, -1
        seen.add((ell, r))
        m = r - ell + 1
        b1 = F.b1_at(S, rule, ell, r)
        b2 = F.b2_at(S, ell, r)
        ln += 1
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
            return 3, L, R, ln, -1
        if ell + dy >= r + dx:
            return 0, L, R, ln, ell + dy
        if r + dx < L or ell + dy > R:
            L += 1
            R -= 1
            return 3, L, R, ln, -1
        ell += dy
        r += dx


def jm_adapt(rule, S, n, W2, first_charge=None, audit=False):
    """Literal epochs until the law would have been right, then jump.

    The machine runs epoch 0 literally.  After a literal epoch whose actual cut
    equals pred_of(b2,b1) (the law would have got it right), it hands over to
    the jump machine.  While the literal cut disagrees with the law it keeps
    running epochs literally -- the repairs are only paid where the law is
    observed to be wrong at this tape's own region.
    """
    sv = S.s
    L, R = 0, n - 1
    reads = n
    epochs = 0
    cuts = 0
    eg = 0
    lit_steps = 0
    lit_epochs = 0
    audit_ep = [] if audit else None
    while True:
        M = R - L + 1
        if M < 3:
            return dict(delta=None, end=1, reads=reads, epochs=epochs,
                        cuts=cuts, eg=eg, lit_steps=lit_steps,
                        lit_epochs=lit_epochs, audit=audit_ep)
        if M <= W2:
            break
        b1 = int(F.b1_at(S, rule, L, R))
        b2 = int(F.b2_at(S, L, R))
        pr = C.pred_of(rule, b2, b1)
        if audit_ep is not None:
            audit_ep.append((M, b2, b1, pr, 1))
        kind, L2, R2, ln, meetpos = _lit_epoch(S, rule, n, L, R)
        epochs += 1
        lit_epochs += 1
        if first_charge is not None and epochs == 1:
            reads += first_charge
            lit_steps += first_charge
        else:
            reads += ln
            lit_steps += ln
        if kind == 0:
            return dict(delta=meetpos / (n - 1.0), end=0, reads=reads,
                        epochs=epochs, cuts=cuts, eg=eg,
                        lit_steps=lit_steps, lit_epochs=lit_epochs,
                        audit=audit_ep)
        if kind == 1:
            return dict(delta=None, end=1, reads=reads, epochs=epochs,
                        cuts=cuts, eg=eg, lit_steps=lit_steps,
                        lit_epochs=lit_epochs, audit=audit_ep)
        if kind == 2:
            return dict(delta=None, end=2, reads=reads, epochs=epochs,
                        cuts=cuts, eg=eg, lit_steps=lit_steps,
                        lit_epochs=lit_epochs, audit=audit_ep)
        # the cut's side
        if L2 > R2:
            return dict(delta=None, end=1, reads=reads, epochs=epochs,
                        cuts=cuts, eg=eg, lit_steps=lit_steps,
                        lit_epochs=lit_epochs, audit=audit_ep)
        dm = M - (R2 - L2 + 1)
        act = 3 if dm == 2 else (1 if L2 > L else 2)
        L, R = L2, R2
        cuts += 1
        reads += dm
        if act != pr:
            continue          # the law was wrong here: keep going literally
        break                 # the law was right: hand over to the jump path

    # ---- the jump path, from the region the literal phase left
    while True:
        M = R - L + 1
        if M < 3:
            return dict(delta=None, end=1, reads=reads, epochs=epochs,
                        cuts=cuts, eg=eg, lit_steps=lit_steps,
                        lit_epochs=lit_epochs, audit=audit_ep)
        if M <= W2:
            ell, r = L, R
            seen = set()
            while True:
                if ell >= r:
                    return dict(delta=ell / (n - 1.0), end=0, reads=reads,
                                epochs=epochs, cuts=cuts, eg=eg,
                                lit_steps=lit_steps, lit_epochs=lit_epochs,
                                audit=audit_ep)
                if R - L + 1 < 3:
                    return dict(delta=None, end=1, reads=reads, epochs=epochs,
                                cuts=cuts, eg=eg, lit_steps=lit_steps,
                                lit_epochs=lit_epochs, audit=audit_ep)
                if (ell, r) in seen:
                    return dict(delta=None, end=2, reads=reads, epochs=epochs,
                                cuts=cuts, eg=eg, lit_steps=lit_steps,
                                lit_epochs=lit_epochs, audit=audit_ep)
                seen.add((ell, r))
                m = r - ell + 1
                bb1 = F.b1_at(S, rule, ell, r)
                bb2 = F.b2_at(S, ell, r)
                reads += 1
                eg += 1
                act2 = TAB[2 * bb2 + bb1]
                dy = 1 if act2 == P.LPLUS else (-1 if act2 == P.LMINUS else 0)
                dx = 1 if act2 == P.RPLUS else (-1 if act2 == P.RMINUS else 0)
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
                                    lit_steps=lit_steps,
                                    lit_epochs=lit_epochs, audit=audit_ep)
                    ell, r = L, R
                    seen.clear()
                    continue
                if ell + dy >= r + dx:
                    return dict(delta=(ell + dy) / (n - 1.0), end=0,
                                reads=reads, epochs=epochs, cuts=cuts, eg=eg,
                                lit_steps=lit_steps, lit_epochs=lit_epochs,
                                audit=audit_ep)
                if r + dx < L or ell + dy > R:
                    L += 1
                    R -= 1
                    if L > R:
                        return dict(delta=None, end=1, reads=reads,
                                    epochs=epochs, cuts=cuts, eg=eg,
                                    lit_steps=lit_steps,
                                    lit_epochs=lit_epochs, audit=audit_ep)
                    ell, r = L, R
                    seen.clear()
                    continue
                ell += dy
                r += dx
        b1 = int(F.b1_at(S, rule, L, R))
        b2 = int(F.b2_at(S, L, R))
        pr = C.pred_of(rule, b2, b1)
        if audit_ep is not None:
            audit_ep.append((M, b2, b1, pr, 0))
        epochs += 1
        cuts += 1
        if pr == 3:
            reads += 2
            L += 1
            R -= 1
        elif pr == 2:
            reads += 1
            R -= 1
        else:
            reads += 1
            L += 1
        if L > R:
            return dict(delta=None, end=1, reads=reads, epochs=epochs,
                        cuts=cuts, eg=eg, lit_steps=lit_steps,
                        lit_epochs=lit_epochs, audit=audit_ep)
