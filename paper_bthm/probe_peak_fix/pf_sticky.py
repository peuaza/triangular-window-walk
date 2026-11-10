# -*- coding: utf-8 -*-
"""pf_sticky -- peak machine: literal epoch 0, then a one-entry cut memory.

If epoch 0's literal cut disagrees with pred_of(b2,b1), the machine records the
signature (b2,b1) -> the side the cut actually fell on, and while that signature
repeats it predicts the observed side; a signature change clears the record.
The rest of the machine is el_cell.jm_walk: jump above W2, literal below.

Cost: n (opening) + epoch 0 (its measured length, or `first_charge` if given)
+ 1 read per cell a later cut gives up + 1 per endgame step.
State: the region, the one-entry memory (2 + 1 bits), the counter.
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(os.path.dirname(HERE), 'probe_en_linear'))

import pf_lib as F
import pf_var as V
import el_cell as C
import pcr_lib as P

TAB = P.TABLE_NP


def jm_sticky(rule, S, n, W2, first_charge=None):
    L, R = 0, n - 1
    reads = n
    cuts = 0
    epochs = 0
    lit_steps = 0
    b1 = int(F.b1_at(S, rule, L, R))
    b2 = int(F.b2_at(S, L, R))
    pr = C.pred_of(rule, b2, b1)
    kind, L2, R2, ln, meetpos = V._lit_epoch(S, rule, n, L, R)
    epochs += 1
    lit_steps += (first_charge if first_charge is not None else ln)
    reads += (first_charge if first_charge is not None else ln)
    if kind == 0:
        return dict(delta=meetpos / (n - 1.0), end=0, reads=reads,
                    cuts=cuts, epochs=epochs, lit_steps=lit_steps,
                    lit_epochs=1, mode=None)
    if kind == 1:
        return dict(delta=None, end=1, reads=reads, cuts=cuts, epochs=epochs,
                    lit_steps=lit_steps, lit_epochs=1, mode=None)
    if kind == 2:
        return dict(delta=None, end=2, reads=reads, cuts=cuts, epochs=epochs,
                    lit_steps=lit_steps, lit_epochs=1, mode=None)
    if L2 > R2:
        return dict(delta=None, end=1, reads=reads, cuts=cuts, epochs=epochs,
                    lit_steps=lit_steps, lit_epochs=1, mode=None)
    dm = (R - L + 1) - (R2 - L2 + 1)
    act = 3 if dm == 2 else (1 if L2 > L else 2)
    L, R = L2, R2
    cuts += 1
    reads += dm
    mode = None
    if act != pr and act != 3:
        mode = [(b2, b1), act]
    while True:
        M = R - L + 1
        if M < 3:
            return dict(delta=None, end=1, reads=reads, cuts=cuts,
                        epochs=epochs, lit_steps=lit_steps, lit_epochs=1,
                        mode=mode)
        if M <= W2:
            ell, r = L, R
            seen = set()
            while True:
                if ell >= r:
                    return dict(delta=ell / (n - 1.0), end=0, reads=reads,
                                cuts=cuts, epochs=epochs, lit_steps=lit_steps,
                                lit_epochs=1, mode=mode)
                if R - L + 1 < 3:
                    return dict(delta=None, end=1, reads=reads, cuts=cuts,
                                epochs=epochs, lit_steps=lit_steps,
                                lit_epochs=1, mode=mode)
                if (ell, r) in seen:
                    return dict(delta=None, end=2, reads=reads, cuts=cuts,
                                epochs=epochs, lit_steps=lit_steps,
                                lit_epochs=1, mode=mode)
                seen.add((ell, r))
                m = r - ell + 1
                bb1 = int(F.b1_at(S, rule, ell, r))
                bb2 = int(F.b2_at(S, ell, r))
                reads += 1
                a = TAB[2 * bb2 + bb1]
                dy = 1 if a == P.LPLUS else (-1 if a == P.LMINUS else 0)
                dx = 1 if a == P.RPLUS else (-1 if a == P.RMINUS else 0)
                if ell + dy < L or r + dx > R:
                    if ell == L and r == R:
                        L += 1
                        R -= 1
                    elif ell + dy < L:
                        L += 1
                    else:
                        R -= 1
                    if L > R:
                        return dict(delta=None, end=1, reads=reads,
                                    cuts=cuts, epochs=epochs,
                                    lit_steps=lit_steps, lit_epochs=1,
                                    mode=mode)
                    ell, r = L, R
                    seen.clear()
                    continue
                if ell + dy >= r + dx:
                    return dict(delta=(ell + dy) / (n - 1.0), end=0,
                                reads=reads, cuts=cuts, epochs=epochs,
                                lit_steps=lit_steps, lit_epochs=1, mode=mode)
                if r + dx < L or ell + dy > R:
                    L += 1
                    R -= 1
                    if L > R:
                        return dict(delta=None, end=1, reads=reads,
                                    cuts=cuts, epochs=epochs,
                                    lit_steps=lit_steps, lit_epochs=1,
                                    mode=mode)
                    ell, r = L, R
                    seen.clear()
                    continue
                ell += dy
                r += dx
        b1 = int(F.b1_at(S, rule, L, R))
        b2 = int(F.b2_at(S, L, R))
        pr = C.pred_of(rule, b2, b1)
        if mode is not None:
            if (b2, b1) == mode[0]:
                pr = mode[1]
            else:
                mode = None
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
            return dict(delta=None, end=1, reads=reads, cuts=cuts,
                        epochs=epochs, lit_steps=lit_steps, lit_epochs=1,
                        mode=mode)
