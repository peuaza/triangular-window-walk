# -*- coding: utf-8 -*-
"""sl2_engine -- pcr_lib's walk, plus one new rule, gated line-for-line on pcr_lib.

WHY A COPY
----------
pcr_lib.walk is a closed @njit over twenty hard-wired rule ids.  The rule this probe
needs does not exist in it.  `add_rule` was not available, so the walk is re-stated
here with the SAME action set, the SAME cut-side semantics (conv 0, "a": the cut
falls on the side the move would have left, double cut on a full window), the SAME
INC accounting and the SAME repeat stamp -- and then CHECKED, not asserted, against
pcr_lib.Runner on the canonical rule: sl2_gate.py demands step counts, m-sum, INC
reads, cut counts and the ending to agree exactly, tape by tape.  A hand-rolled walk
that does not reproduce pcr_lib is void; the gate is what makes this one legal.

TWO ACCOUNTINGS (pcr_lib's, unchanged)
--------------------------------------
INC   the O(log n)-bit machine: one read per cell taken in or given up, one read per
      cell of a re-span after a cut (the reset forces the window back onto the whole
      region and the O(log n)-bit state cannot store it).  The initial build of the
      window aggregates is n + 2(n-1) + 2 + n//2, cr_lib's Win.__init__ charge.
SCAN  m reads per query -- "a pass over a window of length m costs m reads".

THE NEW RULE  `ballast`
-----------------------
varphi(w) = (c(w), z(w)) with

  c(w)  the canonical bit of the window (pcr_lib's R_CANON branch, unchanged);
  z(w)  1 iff  hash(w) = sum_i a_i 3^i mod (2^31-1)  is divisible by n.

f(varphi, n, ell, r, L, R), i.e. a memoryless strategy that sees the geometry:

  (D) L = 0, ell = 0, r = n-1, z = 1   ->  left cut.        the one-shot defect
  (N) L = 0                            ->  ell+1 if c = 1 else r-1.
  (X) L > 0                            ->  the long-epoch walk of the defect.

Read the branch off the geometry: the normal branch never cuts, so L = 0 for the whole
of it; the defect branch cuts at step 1, so L >= 1 for the whole of it.  One bit of
memory, carried by L, which is already part of the state -- the rule stays a
memoryless f(varphi, n, ell, r, L, R) as FORMAL_MODEL.md:151-160 requires, and no
bit has to be added to the tuple.

NORMAL BRANCH COST.  Monotone inward: exactly n-2 steps, one read each plus the
head's drift, no cut, so no re-span.  Theta(n) reads and delta is the meeting point
of a content-driven sequence -- not a constant.
DEFECT BRANCH COST.  Epochs: shrink the window from the r side until it is half the
region, then cut on the side c(w) names; the cut is a jump, so the machine re-reads
the whole new window.  Epoch of width W costs ~W/2 moves + ~W/2 re-span reads, and the
region loses one cell per epoch, so the branch totals sum_W Theta(W) = Theta(n^2).
The cut side is read off the window at every epoch, so the re-span is not a
bookkeeping artefact: the decision needs the window's aggregates, and after a jump the
O(log n)-bit state has to re-acquire them.
"""
import numpy as np
from numba import njit

LPLUS, LMINUS, RPLUS, RMINUS, CUTL, CUTR, CUTB, HALT = range(8)
TABLE_NP = np.array([LMINUS, RPLUS, RMINUS, LPLUS], np.int8)   # q = 2*b2 + b1

R_CANON, R_BALLAST = 0, 1
P31 = 2147483647                    # 2^31 - 1, the hashmod prime


# ------------------------------------------------------------------ predicates
@njit(cache=True, inline='always')
def canon_b1(a, Pa, Paa, PR, ell, r):
    """pcr_lib.b1_of, branch R_CANON, verbatim."""
    m = r - ell + 1
    if m == 1:
        return 1 if a[ell] > 1 else 0
    if m == 2:
        return 1 if a[r] > a[ell] else 0
    S = Pa[r + 1] - Pa[ell]
    Q = Paa[r + 1] - Paa[ell]
    Rr = PR[r] - PR[ell]
    A = S - a[r]
    B = S - a[ell]
    np_ = m - 1
    if np_ * (Q - a[r] * a[r]) - A * A <= 0 or np_ * (Q - a[ell] * a[ell]) - B * B <= 0:
        return 0
    return 1 if np_ * Rr - A * B > 0 else 0


@njit(cache=True)
def tape_hash(a, n):
    """sum_i a_i 3^i mod 2^31-1, one pass, O(log n) bits of accumulator."""
    h = 0
    pw = 1
    for i in range(n):
        h = (h + pw * a[i]) % P31
        pw = (pw * 3) % P31
    return h


@njit(cache=True)
def tape_hashes(X, n):
    nt = X.shape[0]
    out = np.zeros(nt, np.int64)
    for t in range(nt):
        h = 0
        pw = 1
        for i in range(n):
            h = (h + pw * (X[t, i, 0] + X[t, i, 1])) % P31
            pw = (pw * 3) % P31
        out[t] = h
    return out


# ------------------------------------------------------------------------ walk
@njit(cache=True)
def walk2(rule, conv, a, Pa, Paa, PR, stamp, genbox, n, cap):
    """pcr_lib.walk for two rules.  Returns (out, steps, msum, inc, cuts, end)."""
    L, R_ = 0, n - 1
    ell, r = L, R_
    y, x = L, R_                                  # the machine's window
    hy, hx = 0, max(1, n // 2) - 1                # its head
    gen = genbox[0] + 1
    steps = 0
    msum = 0
    cuts = 0
    inc = n + 2 * (n - 1) + 2 + max(1, n // 2)
    out = -1.0
    end = 1
    z = 0                                          # the one-shot defect bit
    while True:
        if ell >= r:
            out = ell / (n - 1.0)
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
        if stamp[key] == gen:
            end = 2
            break
        stamp[key] = gen
        m = r - ell + 1
        msum += m
        # ---- the action
        if rule == R_CANON:
            b1 = canon_b1(a, Pa, Paa, PR, ell, r)
            b2 = 0
            k = m >> 1
            if k == 0:
                k = 1
            b2 = 1 if (Pa[r + 1] - Pa[ell]) * k > (Pa[ell + k] - Pa[ell]) * m else 0
            act = TABLE_NP[2 * b2 + b1]
            cv = conv
        else:
            if L == 0:
                if ell == 0 and r == n - 1:
                    # the opening window is the whole tape: the defect is decided once
                    z = 1 if tape_hash(a, n) % n == 0 else 0
                if z == 1:
                    act = CUTL
                else:
                    act = LPLUS if canon_b1(a, Pa, Paa, PR, ell, r) == 1 else RMINUS
            else:
                W = R_ - L + 1
                if W <= 4:
                    act = LPLUS if canon_b1(a, Pa, Paa, PR, ell, r) == 1 else RMINUS
                elif 2 * m > W:
                    act = RMINUS
                else:
                    act = CUTL if canon_b1(a, Pa, Paa, PR, ell, r) == 0 else CUTR
            cv = 0
        if act >= CUTL:
            cuts += 1
            if act == CUTL:
                L += 1
            elif act == CUTR:
                R_ -= 1
            else:
                L += 1
                R_ -= 1
            if L > R_:
                end = 1
                break
            ell, r = L, R_
            gen += 1
            inc += abs(y - L) + abs(x - R_)
            y, x = L, R_
            kk = (x - y + 1) >> 1
            if kk == 0:
                kk = 1
            ty, tx = y, y + kk - 1
            inc += abs(hy - ty) + abs(hx - tx)
            hy, hx = ty, tx
            continue
        dy = 1 if act == LPLUS else (-1 if act == LMINUS else 0)
        dx = 1 if act == RPLUS else (-1 if act == RMINUS else 0)
        if ell + dy < L or r + dx > R_:
            cuts += 1
            if cv == 1:
                if ell + dy < L:
                    R_ -= 1
                else:
                    L += 1
            else:
                if ell == L and r == R_:
                    L += 1
                    R_ -= 1
                    y += 1
                    x -= 1
                    inc += 2
                elif ell + dy < L:
                    L += 1
                else:
                    R_ -= 1
            if L > R_:
                end = 1
                break
            ell, r = L, R_
            gen += 1
            inc += abs(y - L) + abs(x - R_)
            y, x = L, R_
            kk = (x - y + 1) >> 1
            if kk == 0:
                kk = 1
            ty, tx = y, y + kk - 1
            inc += abs(hy - ty) + abs(hx - tx)
            hy, hx = ty, tx
            continue
        if ell + dy >= r + dx:
            out = (ell + dy) / (n - 1.0)
            end = 0
            break
        ell += dy
        r += dx
        y, x = ell, r
        inc += 1
        kk = (x - y + 1) >> 1
        if kk == 0:
            kk = 1
        ty, tx = y, y + kk - 1
        inc += abs(hy - ty) + abs(hx - tx)
        hy, hx = ty, tx
    genbox[0] = gen
    return out, steps, msum, inc, cuts, end


class Runner2(object):
    """Same contract as pcr_lib.Runner, plus the ballast rule."""

    def __init__(self, S):
        self.S = S
        self.stamp = np.zeros((S.n + 1) * (S.n + 1), np.int32)
        self.genbox = np.zeros(1, np.int64)

    def run(self, rule, cap=None, conv=None):
        S = self.S
        n = S.n
        if cap is None:
            cap = 40 * n * n + 100
        cv = 1 if rule == 'canonical' else 0
        if conv is not None:
            cv = conv
        rid = R_CANON if rule == 'canonical' else R_BALLAST
        return walk2(rid, cv, S.a, S.Pa, S.Paa, S.PR, self.stamp, self.genbox,
                     n, cap)
