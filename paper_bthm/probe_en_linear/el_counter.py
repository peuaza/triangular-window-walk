# -*- coding: utf-8 -*-
"""el_counter -- the sliding devices, checked against the rules window by window.

For dense2, peak, trip, oddrun, canonical, LAS and KCHG the report claims the
rule's bit is a function of O(1) counters over the region, maintainable with one
read per cell given up.  This file implements the counters and drags the region
from the whole tape down to a point along a random cut sequence, comparing the
counter's bit with pcr_lib's b1 at every step.

The three non-obvious ones are written out:

  KCHG   b1 = [2*(x - ell) > m] with x the lower median sign-change position.
         Equivalently [cnt_left <= (c+1)//2 - 1] where c is the number of sign
         changes in [ell+1, r] and cnt_left those with p - ell <= m//2.  So the
         median test is a comparison of two counts.
  LAS    las = 1 + (#sign alternations between consecutive nonzero s in the
         window); the window's first/last nonzero and the alternation count are
         the counters, and a scan to the next nonzero is charged to the zeros it
         crosses (the region only shrinks, so each zero is crossed once).
  oddrun b1 = [2*c > m] with c = #odd maximal constant runs of a fully inside
         the window plus the parities of the two partial boundary runs.

Usage: python el_counter.py     -> prints mismatches per rule
"""
import sys

import numpy as np

sys.path.insert(0, '.')
import el_cell as C
import el_lib as E
import pcr_lib as P


class Dev(object):
    """Counters over the region [L,R]; `bit()` must equal pcr_lib's b1."""

    def __init__(self, rule, S, n):
        self.rule = rule
        self.S = S
        self.n = n
        self.L = 0
        self.R = n - 1
        a = S.a
        self.a = a
        self.s = S.s
        self.reads = 0
        m = n
        if rule == 'dense2':
            self.c2 = int((a == 2).sum())
        elif rule == 'peak':
            self.pk = 0
            for i in range(1, n - 1):
                if a[i - 1] < a[i] and a[i] > a[i + 1]:
                    self.pk += 1
        elif rule == 'trip':
            self.tp = 0
            for i in range(0, n - 2):
                if a[i] == a[i + 1] == a[i + 2]:
                    self.tp += 1
        elif rule == 'oddrun':
            self._build_runs()
        elif rule == 'canonical':
            self.Sa = int(a.sum())
            self.Q = int((a.astype(np.int64) ** 2).sum())
            self.Rr = int((a[:-1].astype(np.int64) * a[1:]).sum()) if n > 1 else 0
        elif rule == 'LAS':
            self._build_las()
        elif rule == 'KCHG':
            self._build_kchg()
        else:
            raise ValueError(rule)

    # ---- oddrun
    def _build_runs(self):
        a = self.a
        self.runs = []          # (start, end) maximal constant runs
        i = 0
        n = self.n
        while i < n:
            j = i
            while j + 1 < n and a[j + 1] == a[i]:
                j += 1
            self.runs.append((i, j))
            i = j + 1
        self.ri = 0
        self.rj = len(self.runs) - 1

    # ---- LAS
    def _build_las(self):
        s = self.s
        self.nz = [i for i in range(self.n) if s[i] != 0]
        self.alt = 0
        for k in range(1, len(self.nz)):
            if s[self.nz[k]] != s[self.nz[k - 1]]:
                self.alt += 1
        self.i0 = 0
        self.i1 = len(self.nz) - 1

    # ---- KCHG
    def _build_kchg(self):
        s = self.s
        self.ch = [p for p in range(1, self.n)
                   if s[p - 1] != 0 and s[p] != 0 and s[p - 1] * s[p] < 0]

    # ------------------------------------------------------------- bits
    def bit(self):
        rule = self.rule
        L, R = self.L, self.R
        m = R - L + 1
        a = self.a
        if rule == 'dense2':
            return 1 if 4 * self.c2 > m else 0
        if rule == 'peak':
            return 1 if (m >= 3 and 6 * self.pk > m) else 0
        if rule == 'trip':
            return 1 if (m >= 3 and 6 * self.tp > m) else 0
        if rule == 'oddrun':
            c = 0
            ri, rj = self.ri, self.rj
            if ri == rj:
                return 1 if (2 * (1 if (m & 1) else 0) > m) else 0
            if ((self.runs[ri][1] - L + 1) & 1):
                c += 1
            if ((R - self.runs[rj][0] + 1) & 1):
                c += 1
            for k in range(ri + 1, rj):
                s0, e0 = self.runs[k]
                if ((e0 - s0 + 1) & 1):
                    c += 1
            return 1 if 2 * c > m else 0
        if rule == 'canonical':
            if m == 1:
                return 1 if a[L] > 1 else 0
            if m == 2:
                return 1 if a[R] > a[L] else 0
            S = self.Sa
            Q = self.Q
            Rr = self.Rr
            A = S - int(a[R])
            B = S - int(a[L])
            np_ = m - 1
            if np_ * (Q - int(a[R]) ** 2) - A * A <= 0:
                return 0
            if np_ * (Q - int(a[L]) ** 2) - B * B <= 0:
                return 0
            return 1 if np_ * Rr - A * B > 0 else 0
        if rule == 'LAS':
            if self.i1 < self.i0:
                return 0
            las = 1 + self.alt
            return 1 if 4 * las > m else 0
        if rule == 'KCHG':
            c = 0
            cnt_left = 0
            half = ell_half = L + m // 2
            for p in self.ch:
                if L + 1 <= p <= R:
                    c += 1
                    if p <= half:
                        cnt_left += 1
            if c == 0:
                return 0
            k = (c + 1) // 2
            return 1 if cnt_left <= k - 1 else 0
        raise ValueError(rule)

    # ---------------------------------------------------------- updates
    def remove_left(self):
        rule = self.rule
        L = self.L
        a = self.a
        m = self.R - L + 1
        if rule == 'dense2':
            if a[L] == 2:
                self.c2 -= 1
        elif rule == 'peak':
            for i in (L + 1,):
                if L + 1 <= i <= self.R - 1:
                    if a[i - 1] < a[i] and a[i] > a[i + 1]:
                        self.pk -= 1
        elif rule == 'trip':
            for i in (L,):
                if L <= i <= self.R - 2:
                    if a[i] == a[i + 1] == a[i + 2]:
                        self.tp -= 1
        elif rule == 'oddrun':
            if self.ri == self.rj:
                pass
            elif L >= self.runs[self.ri][1]:
                self.ri += 1
                if self.ri > self.rj:
                    self.rj = self.ri
        elif rule == 'canonical':
            if m >= 3:
                self.Sa -= int(a[L])
                self.Q -= int(a[L]) ** 2
                self.Rr -= int(a[L]) * int(a[L + 1])
        elif rule == 'LAS':
            if self.i1 >= self.i0 and self.nz[self.i0] == L:
                if self.i0 + 1 <= self.i1:
                    if self.s[self.nz[self.i0 + 1]] != self.s[self.nz[self.i0]]:
                        self.alt -= 1
                    self.i0 += 1
                else:
                    self.i0 += 1
        self.L += 1

    def remove_right(self):
        rule = self.rule
        R = self.R
        a = self.a
        m = R - self.L + 1
        if rule == 'dense2':
            if a[R] == 2:
                self.c2 -= 1
        elif rule == 'peak':
            i = R - 1
            if self.L + 1 <= i <= R - 1:
                if a[i - 1] < a[i] and a[i] > a[i + 1]:
                    self.pk -= 1
        elif rule == 'trip':
            i = R - 2
            if self.L <= i <= R - 2:
                if a[i] == a[i + 1] == a[i + 2]:
                    self.tp -= 1
        elif rule == 'oddrun':
            if self.ri == self.rj:
                pass
            elif R <= self.runs[self.rj][0]:
                self.rj -= 1
                if self.rj < self.ri:
                    self.ri = self.rj
        elif rule == 'canonical':
            if m >= 3:
                self.Sa -= int(a[R])
                self.Q -= int(a[R]) ** 2
                self.Rr -= int(a[R - 1]) * int(a[R])
        elif rule == 'LAS':
            if self.i1 >= self.i0 and self.nz[self.i1] == R:
                if self.i0 <= self.i1 - 1:
                    if self.s[self.nz[self.i1 - 1]] != self.s[self.nz[self.i1]]:
                        self.alt -= 1
                    self.i1 -= 1
                else:
                    self.i1 -= 1
        self.R -= 1


def check(rule, n, nt=3, seed=0):
    rng = np.random.default_rng(seed)
    bad = 0
    tot = 0
    for i in range(nt):
        X = E.tape(n, i)
        S = P.St(X)
        D = Dev(rule, S, n)
        while D.L < D.R:
            tot += 1
            b = D.bit()
            t = C.b1_at(S, rule, D.L, D.R)
            if b != t:
                bad += 1
                if bad < 5:
                    print('   %s [%d,%d] dev %d true %d' % (rule, D.L, D.R, b, t))
            # a random single cut, keeping L <= R
            if D.R - D.L <= 1 or rng.random() < 0.5:
                D.remove_left()
            else:
                D.remove_right()
    return bad, tot


if __name__ == '__main__':
    for rule in ('dense2', 'peak', 'trip', 'oddrun', 'canonical', 'LAS', 'KCHG'):
        b, t = check(rule, 64)
        print('%-9s windows %6d  mismatches %d' % (rule, t, b))
