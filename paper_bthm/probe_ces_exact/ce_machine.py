# -*- coding: utf-8 -*-
"""ce_machine -- CES with patchable exact-epoch simulation, read-counted.

One machine, three knobs:

  W      endgame buffer: while M > W, shortcut with the (b2,b1) outcome rule;
         at M <= W buffer the window and run the literal walk on it.
  W2     if M <= W2 (and > W) run the epoch exactly instead of shortcutting.
         W2=None -> pure CES.  W2=n -> the literal walk, O(n^2) reads.
  K      prefix verification: at every shortcut epoch, simulate the first K
         steps on a copy; if the epoch terminates within K use that outcome,
         otherwise take the CES cut.  K=None -> no check.

Cycle detection inside an exact epoch is Brent's, so the machine keeps O(log n)
bits (no visited set).  All tape access goes through read counters.

usage (self-test): ce_machine.py <n> <nt>
"""
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
for _p in (os.path.join(ROOT, 'probe_e_uniform'),
           os.path.join(ROOT, 'probe_kappa_pred'),
           os.path.join(ROOT, 'probe_ologn2')):
    if _p not in sys.path:
        sys.path.insert(0, _p)


class Win(object):
    """Window [ell, r] with S, C, H, PC maintained incrementally.

    A copy() is O(1) and costs no reads; the cells are shared.
    """

    __slots__ = ('a', 'n', 'rule', 'ell', 'r', 'M', 'k', 'S', 'C', 'H', 'PC', 'reads')

    def __init__(self, a, L, R, rule, reads=None):
        self.a = a
        self.n = len(a)
        self.rule = rule
        self.ell, self.r = L, R
        M = R - L + 1
        self.M = M
        self.k = max(1, M >> 1)
        S = C = H = PC = 0
        prev = None
        for i in range(L, R + 1):
            v = int(a[i])
            S += v
            if i < L + self.k:
                H += v
            if prev is not None and prev != v:
                C += 1
            prev = v
        if rule == 'pattern':
            for j in range(L, max(L, R - 1)):
                if a[j] == 2 and a[j + 1] == 1 and a[j + 2] == 1:
                    PC += 1
        self.S, self.C, self.H, self.PC = S, C, H, PC
        self.reads = M          # the opening pass: one read per cell

    def copy(self):
        w = Win.__new__(Win)
        for s in self.__slots__:
            setattr(w, s, getattr(self, s))
        return w

    def rd(self, i):
        self.reads += 1
        return int(self.a[i])

    # ---------------------------------------------------------- b2 / b1
    def b2(self):
        return self.S * self.k > self.H * self.M

    def b1(self):
        if self.rule == 'b0':
            return False
        if self.rule == 'runs':
            return self.C * 2 > self.M
        if self.rule == 'pattern':
            return self.PC > 0
        raise ValueError(self.rule)

    # ---------------------------------------------------------- moves
    def trim_left(self):
        """drop ell."""
        a, ell, r, M = self.a, self.ell, self.r, self.M
        x = self.rd(ell)
        y = self.rd(ell + 1) if M >= 2 else 0
        self.S -= x
        if M >= 2:
            self.C -= (x != y)
        if M % 2 == 0:
            self.H -= x
        else:
            self.H += self.rd(ell + self.k) - x
        if self.rule == 'pattern' and M >= 3:
            if x == 2 and y == 1 and self.rd(ell + 2) == 1:
                self.PC -= 1
        self.ell = ell + 1
        self.M = M - 1
        self.k = max(1, self.M >> 1)

    def trim_right(self):
        """drop r."""
        a, ell, r, M = self.a, self.ell, self.r, self.M
        x = self.rd(r)
        y = self.rd(r - 1) if M >= 2 else 0
        self.S -= x
        if M >= 2:
            self.C -= (x != y)
        if M % 2 == 0:
            self.H -= self.rd(ell + self.k - 1)
        if self.rule == 'pattern' and M >= 3:
            if y == 1 and x == 1 and self.rd(r - 2) == 2:
                self.PC -= 1
        self.r = r - 1
        self.M = M - 1
        self.k = max(1, self.M >> 1)

    def ext_left(self):
        """add ell-1."""
        ell, r, M = self.ell, self.r, self.M
        x = self.rd(ell - 1)
        self.S += x
        self.C += (x != self.rd(ell))
        if M % 2 == 0:
            self.H += x - self.rd(ell + self.k - 1)
        else:
            self.H += x
        if self.rule == 'pattern':
            if x == 2 and self.rd(ell) == 1 and self.rd(ell + 1) == 1:
                self.PC += 1
        self.ell = ell - 1
        self.M = M + 1
        self.k = max(1, self.M >> 1)

    def ext_right(self):
        """add r+1."""
        ell, r, M = self.ell, self.r, self.M
        x = self.rd(r + 1)
        self.S += x
        self.C += (x != self.rd(r))
        if M % 2 == 1:
            self.H += self.rd(ell + self.k)
        if self.rule == 'pattern':
            if self.rd(r - 1) == 2 and self.rd(r) == 1 and x == 1:
                self.PC += 1
        self.r = r + 1
        self.M = M + 1
        self.k = max(1, self.M >> 1)

    # ------------------------------------------------------- region cuts
    def cut_left(self):
        a, L, M, k = self.a, self.ell, self.M, self.k
        x = self.rd(L)
        y = self.rd(L + 1) if M >= 2 else 0
        self.S -= x
        if M >= 2:
            self.C -= (x != y)
        if self.rule == 'pattern' and M >= 3:
            if x == 2 and y == 1 and self.rd(L + 2) == 1:
                self.PC -= 1
        if M % 2 == 0:
            self.H -= x
        else:
            self.H += self.rd(L + k) - x
        self.ell = L + 1
        self.r = self.ell + M - 2
        self.M = M - 1
        self.k = max(1, self.M >> 1)

    def cut_right(self):
        a, L, M, k = self.a, self.ell, self.M, self.k
        x = self.rd(L + M - 1)
        y = self.rd(L + M - 2) if M >= 2 else 0
        self.S -= x
        if M >= 2:
            self.C -= (x != y)
        if self.rule == 'pattern' and M >= 3:
            if y == 1 and x == 1 and self.rd(L + M - 3) == 2:
                self.PC -= 1
        if M % 2 == 0:
            self.H -= self.rd(L + k - 1)
        self.ell = L
        self.r = L + M - 2
        self.M = M - 1
        self.k = max(1, self.M >> 1)

    def double_cut(self):
        self.cut_left()
        # after cut_left the region is [L+1, R]; cut its right end
        self.cut_right()


def _event(w, L, R):
    """One move of the reference walk; returns (kind, payload) or ('move',)."""
    m = w.M
    b2 = w.b2(); b1 = w.b1()
    if b2:
        dy, dx = (1, 0) if b1 else (0, -1)
    else:
        dy, dx = (0, 1) if b1 else (-1, 0)
    ell, r = w.ell, w.r
    if ell + dy < L or r + dx > R:
        if ell == L and r == R:
            return ('D', None)
        return ('L' if ell + dy < L else 'R', None)
    if ell + dy >= r + dx:
        return ('meet', ell + dy)
    if r + dx < L or ell + dy > R:
        return ('D', None)
    return ('move', (dy, dx))


def sim_epoch(w, L, R, budget=None):
    """Run one epoch from w (window == region), Brent cycle detection.

    Returns (oc, meet_at, steps); oc in {'D','L','R','meet','cycle','budget'}.
    The caller passes a copy when it wants to keep the region state.
    """
    steps = 0
    power = 1
    t_ell, t_r = w.ell, w.r
    while True:
        if w.ell >= w.r:
            return ('meet', w.ell, steps)
        ev, payload = _event(w, L, R)
        if ev == 'meet':
            return ('meet', payload, steps + 1)
        if ev in ('D', 'L', 'R'):
            return (ev, None, steps + 1)
        steps += 1
        dy, dx = payload
        if dy == 1:
            w.trim_left()
        elif dy == -1:
            w.ext_left()
        elif dx == 1:
            w.ext_right()
        else:
            w.trim_right()
        if budget is not None and steps >= budget:
            return ('budget', None, steps)
        if w.ell == t_ell and w.r == t_r:
            return ('cycle', None, steps)
        if steps == power:
            t_ell, t_r = w.ell, w.r
            power *= 2


def apply_cut(w, oc):
    if oc == 'D':
        w.double_cut()
    elif oc == 'R':
        w.cut_right()
    elif oc == 'L':
        w.cut_left()
    else:
        raise ValueError(oc)


def buffer_end(w, n, rule):
    """Exact walk on the M <= W window: W reads, O(W) bits, free compute."""
    L0, M = w.ell, w.M
    a = w.a
    buf = [int(a[L0 + i]) for i in range(M)]
    reads = M + w.reads
    if M == 1:
        return L0 / (n - 1), reads
    if M < 3:
        return np.nan, reads
    P = [0] * (M + 1)
    for i, v in enumerate(buf):
        P[i + 1] = P[i] + v
    ch = [0] * M
    for i in range(1, M):
        ch[i] = ch[i - 1] + (buf[i - 1] != buf[i])
    lw, rw = 0, M - 1
    ell, r = 0, M - 1
    seen = set()
    while True:
        if ell >= r:
            return (L0 + ell) / (n - 1), reads
        if rw - lw + 1 < 3:
            return np.nan, reads
        if (ell, r) in seen:
            return np.nan, reads
        seen.add((ell, r))
        m = r - ell + 1
        k = max(1, m >> 1)
        b2 = (P[ell + m] - P[ell]) * k > (P[ell + k] - P[ell]) * m
        if rule == 'runs':
            b1 = (ch[r] - ch[ell]) * 2 > m
        elif rule == 'pattern':
            b1 = any(buf[t] == 2 and buf[t + 1] == 1 and buf[t + 2] == 1
                     for t in range(ell, r - 1))
        else:
            b1 = False
        if b2:
            dy, dx = (1, 0) if b1 else (0, -1)
        else:
            dy, dx = (0, 1) if b1 else (-1, 0)
        if ell + dy < lw or r + dx > rw:
            if ell == lw and r == rw:
                lw, rw = lw + 1, rw - 1
            elif ell + dy < lw:
                lw += 1
            else:
                rw -= 1
            ell, r = lw, rw
            seen.clear()
            continue
        if ell + dy >= r + dx:
            return (L0 + ell + dy) / (n - 1), reads
        if r + dx < lw or ell + dy > rw:
            lw, rw = lw + 1, rw - 1
            ell, r = lw, rw
            seen.clear()
            continue
        ell, r = ell + dy, r + dx


def machine(a, rule, W=14, W2=None, K=None, oracle_term=None):
    """Run the machine; return (delta_or_nan, reads, info)."""
    n = len(a)
    w = Win(a, 0, n - 1, rule)
    ep = 0
    info = dict(epochs=0, exact_epochs=0, prefix_cuts=0, forced=0,
                terminal=None)
    while w.M > max(W, 2):
        M = w.M
        base = w.reads
        if W2 is not None and M <= W2:
            trial = w.copy()
            oc, meet, st = sim_epoch(trial, w.ell, w.r)
            w.reads += trial.reads - base
            info['exact_epochs'] += 1
            if oc in ('meet', 'cycle'):
                info['terminal'] = oc
                d = meet / (n - 1) if oc == 'meet' else np.nan
                return d, w.reads, info
            apply_cut(w, oc)
        else:
            b2 = w.b2()
            b1 = w.b1()
            oc = 'D' if not b2 else ('R' if b1 else 'L')
            if oracle_term is not None and ep in oracle_term:
                trial = w.copy()
                toc, tmeet, st = sim_epoch(trial, w.ell, w.r)
                w.reads += trial.reads - base
                info['forced'] += 1
                if toc in ('meet', 'cycle'):
                    info['terminal'] = toc
                    d = tmeet / (n - 1) if toc == 'meet' else np.nan
                    return d, w.reads, info
                oc = toc
            elif K is not None:
                trial = w.copy()
                toc, tmeet, tst = sim_epoch(trial, w.ell, w.r, budget=K)
                w.reads += trial.reads - base
                if toc != 'budget':
                    info['prefix_cuts'] += 1
                    if toc in ('meet', 'cycle'):
                        info['terminal'] = toc
                        d = tmeet / (n - 1) if toc == 'meet' else np.nan
                        return d, w.reads, info
                    oc = toc
            apply_cut(w, oc)
        ep += 1
        info['epochs'] = ep
    d, reads = buffer_end(w, n, rule)
    info['terminal'] = 'buffer'
    return d, reads, info
