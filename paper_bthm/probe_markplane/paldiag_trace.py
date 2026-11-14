# -*- coding: utf-8 -*-
"""paldiag_trace -- palmach2.PalMach2, copied, with per-comparison bookkeeping.

The machine is unchanged (same _known/_push/_seed/_scan/_put/__call__, same
palfast.walk).  Added counters only:
    cmp        comparisons made, i.e. b1 reads
    pairs      dict (s, i) -> times compared     (repeat-comparison census)
    repeats    comparisons of a pair already compared
    cellreads  single-cell value reads (the "one side known" branch and _seed)
    cellrep    single-cell reads of a cell that was read before
    seedreads  reads spent in _seed
Everything else must equal palmach2 exactly; mk_gate checks that.
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
PAPER = os.path.dirname(HERE)
for d in ('probe_cert',):
    p = os.path.join(PAPER, d)
    if p not in sys.path:
        sys.path.insert(0, p)

import palmach2                                                 # noqa: E402


class PalTrace(object):

    def __init__(self, a, K=1, B=4, tau=4, budget_mult=4, lmax=None):
        self.a = a
        self.K = K
        self.B = B
        self.tau = tau
        self.slots = []
        self.blks = []
        self.budget = budget_mult * len(a)
        self.reads = 0
        self.calls = 0
        self.scans = 0
        self.hits = 0
        self.skips = 0
        self.seeds = 0
        self.zero = 0
        self.maxslots = 0
        self.lmax = lmax
        self.lfast = 0
        # instrumentation
        self.cmp = 0
        self.pairs = {}
        self.repeats = 0
        self.cellreads = 0
        self.cellseen = set()
        self.cellrep = 0
        self.seedreads = 0

    # ---------------- known values ----------------
    def _known(self, i):
        for lo, hi, c in self.blks:
            if lo <= i <= hi:
                return c
        return None

    def _push(self, lo, hi, c):
        for j in range(len(self.blks)):
            if self.blks[j][0] <= lo and hi <= self.blks[j][1]:
                return
        self.blks.insert(0, (lo, hi, c))
        if len(self.blks) > self.B:
            self.blks.pop()

    def _seed(self, p):
        if self.budget <= 0 or self._known(p) is not None:
            return
        a = self.a
        n = len(a)
        c = a[p]
        self.reads += 1
        self.seedreads += 1
        self._note_cell(p)
        self.budget -= 1
        lo = hi = p
        while hi + 1 < n and self.budget > 0:
            self.reads += 1
            self.seedreads += 1
            self._note_cell(hi + 1)
            self.budget -= 1
            if a[hi + 1] != c:
                break
            hi += 1
        while lo - 1 >= 0 and self.budget > 0:
            self.reads += 1
            self.seedreads += 1
            self._note_cell(lo - 1)
            self.budget -= 1
            if a[lo - 1] != c:
                break
            lo -= 1
        self.seeds += 1
        self._push(lo, hi, c)

    def _note_cell(self, i):
        self.cellreads += 1
        if i in self.cellseen:
            self.cellrep += 1
        else:
            self.cellseen.add(i)

    def _cmp(self, s, i):
        self.cmp += 1
        k = (s, i)
        if k in self.pairs:
            self.repeats += 1
            self.pairs[k] += 1
        else:
            self.pairs[k] = 1

    # ---------------- scan ----------------
    def _scan(self, s, i0, i1):
        a = self.a
        blks = self.blks
        self.scans += 1
        c = 0
        p = -1
        for i in range(i0, i1 + 1):
            j = s - i
            ci = None
            for lo, hi, v in blks:
                if lo <= i <= hi:
                    ci = v
                    break
            cj = None
            for lo, hi, v in blks:
                if lo <= j <= hi:
                    cj = v
                    break
            if ci is not None and cj is not None:
                if ci != cj:
                    p = i
                    break
                self.skips += 1
                self.zero += 1
                continue
            c += 1
            if ci is not None:
                self._note_cell(j)
                if a[j] != ci:
                    p = i
                    break
            elif cj is not None:
                self._note_cell(i)
                if a[i] != cj:
                    p = i
                    break
            else:
                self._cmp(s, i)
                if a[i] != a[j]:
                    p = i
                    break
        self.reads += c
        return p, c

    def _put(self, s, lo, pos):
        sl = self.slots
        for j in range(len(sl)):
            if sl[j][0] == s:
                sl.pop(j)
                break
        sl.insert(0, (s, lo, pos))
        if len(sl) > self.K:
            sl.pop()
        if len(sl) > self.maxslots:
            self.maxslots = len(sl)

    # ---------------- query ----------------
    def __call__(self, ell, r):
        self.calls += 1
        s = ell + r
        mid = (s - 1) >> 1
        if ell > mid:
            return True
        if self.lmax is not None and r - ell + 1 > self.lmax:
            self.reads += 1
            self.lfast += 1
            return False
        blks = self.blks
        for lo, hi, v in blks:
            if lo <= ell and r <= hi:
                self.hits += 1
                return True
        sl = self.slots
        idx = -1
        for j in range(len(sl)):
            if sl[j][0] == s:
                idx = j
                break
        cost = 0
        if idx >= 0:
            _, lo, pos = sl[idx]
            if ell >= lo:
                if idx:
                    sl.insert(0, sl.pop(idx))
                if pos < 0:
                    self.hits += 1
                    return True
                if ell <= pos:
                    self.hits += 1
                    return False
                p, cost = self._scan(s, ell, mid)
            else:
                p, cost = self._scan(s, ell, lo - 1)
                if p < 0:
                    if pos < 0:
                        self.hits += 1
                        return True
                    p = pos
        else:
            p, cost = self._scan(s, ell, mid)
        if p < 0:
            self._put(s, ell, -1)
        else:
            self._put(s, ell, p)
        if cost >= self.tau:
            self._seed(ell)
            self._seed(r)
        return p < 0


def walk(a, o, b2_mode='incr'):
    # palmach2.walk_mach2 calls exactly this: palmach2.pf.walk(a, oracle, mode)
    return palmach2.pf.walk(a, o, b2_mode)
