# -*- coding: utf-8 -*-
"""mk_pal -- palmach2.PalMach2 with a one-bit mark plane, four policies.

The machine is palmach2.PalMach2 (probe_cert/palmach2.py): diagonal certificate
in one slot, up to B=4 constant-run blocks, 4n seed budget, charged 1 read per
comparison ("diagonal certificate M: 1 per cell touched").  Copied, not edited;
the only change is that every charged access goes through a Plane.

policies
  none        identical to the baseline machine; writes nothing.  mk_gate checks
              that reads == palmach2.walk_mach2 on all 300 tapes at every n.
  cert        "mark cells already certified so a later pass can skip them":
              a comparison that comes out equal marks both cells (1 ink each).
              A later scan may consult those marks instead of comparing -- but
              only while the marks are known to belong to the diagonal it is
              scanning, which the machine tracks with one flag.  Comparisons
              skipped by a mark consultation are counted.
  cert_naive  the same marks, trusted by diagonal-agnostic position only (the
              obvious implementation).  Expected to lose exactness; run to
              measure how often, then discarded.
  run         "mark the cells whose constant run has been established": _seed
              walks a constant run cell by cell, so it marks every cell whose
              value continues the run (1 ink each).  A scan may then read a mark
              (1 read) instead of the cell's value (1 read) to decide whether
              the cell continues a certified run.
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
PAPER = os.path.dirname(HERE)
for d in ('probe_cert', 'probe_triple_nat'):
    p = os.path.join(PAPER, d)
    if p not in sys.path:
        sys.path.insert(0, p)

import numpy as np                                              # noqa: E402
import palmach2                                                 # noqa: E402
import mk_lib                                                   # noqa: E402


class PalMarked(object):

    def __init__(self, a, policy='none', K=1, B=4, tau=4, budget_mult=4,
                 lmax=None):
        self.policy = policy
        self.pl = mk_lib.Tape(a) if policy == 'none' else mk_lib.Plane(a)
        self.a = a
        self.K = K
        self.B = B
        self.tau = tau
        self.slots = []
        self.blks = []
        self.budget = budget_mult * len(a)
        self.calls = 0
        self.scans = 0
        self.hits = 0
        self.skips = 0
        self.seeds = 0
        self.zero = 0
        self.maxslots = 0
        self.lmax = lmax
        self.lfast = 0
        # plane bookkeeping (O(1) numbers, i.e. O(log n) bits)
        self.cl_s = -1            # diagonal of the last clean scan
        self.cl_lo = 0
        self.cl_mid = -1
        self.cl_clean = False
        self.consults = 0         # scans that consulted the marks
        self.mark_skips = 0       # comparisons skipped because a mark answered
        self.mark_writes = 0
        self.stale_seen = 0       # scans that found marks but could not use them

    @property
    def reads(self):
        return self.pl.reads

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
        self.pl.val(p)
        self.budget -= 1
        lo = hi = p
        if self.policy == 'run':
            self.pl.wr(p, 1)                 # the run's base cell
            self.mark_writes += 1
        while hi + 1 < n and self.budget > 0:
            v = self.pl.val(hi + 1)
            self.budget -= 1
            if self.policy == 'run':
                self.pl.wr(hi + 1, 1 if v == c else 0)
                self.mark_writes += 1
            if v != c:
                break
            hi += 1
        while lo - 1 >= 0 and self.budget > 0:
            v = self.pl.val(lo - 1)
            self.budget -= 1
            if self.policy == 'run':
                self.pl.wr(lo - 1, 1 if v == c else 0)
                self.mark_writes += 1
            if v != c:
                break
            lo -= 1
        self.seeds += 1
        self._push(lo, hi, c)

    # ---------------- scan ----------------
    def _scan(self, s, i0, i1):
        a = self.a
        blks = self.blks
        self.scans += 1
        pol = self.policy
        consult = False
        if pol == 'cert' and self.cl_clean and self.cl_s == s \
                and i0 >= self.cl_lo and i1 <= self.cl_mid:
            # the previous scan of this same diagonal was clean over
            # [cl_lo, cl_mid]; a window nested in it at the same s is certified.
            # The certificate is only trusted after one mark read.
            consult = True
            self.consults += 1
        elif pol == 'cert_naive':
            consult = True
            self.consults += 1
        c = 0
        p = -1
        r0 = self.pl.reads
        if consult and pol == 'cert':
            if self.pl.rd(i0)[1]:
                self.mark_skips += 1
                return -1, self.pl.reads - r0
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
            if consult:
                mi = self.pl.rd(i)[1]
                if mi:
                    mj = self.pl.rd(j)[1]
                    if mj:
                        self.mark_skips += 1      # both certified: 0 comparisons
                        continue
            if ci is not None:
                if self.pl.val(j) != ci:
                    p = i
                    break
            elif cj is not None:
                if self.pl.val(i) != cj:
                    p = i
                    break
            else:
                if not self.pl.cmp(i, j):
                    p = i
                    break
                if pol == 'cert' or pol == 'cert_naive':
                    self.pl.wr(i, 1)
                    self.pl.wr(j, 1)
                    self.mark_writes += 2
        if pol == 'cert':
            # remember only the certificate the slot could hold anyway
            self.cl_s = s
            self.cl_lo = i0
            self.cl_mid = i1
            self.cl_clean = (p < 0 and i1 == (s - 1) >> 1)
        return p, self.pl.reads - r0

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
            self.pl.val(ell)
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


def run_pal(a, policy='none', **kw):
    """palmach2.walk_mach2's own control flow: it builds a PalMach2 and calls
    palmach2.pf.walk(a, oracle, b2_mode).  This passes the plane machine to the
    same function palmach2 itself calls; mk_gate proves the result equals
    palmach2.walk_mach2 tape by tape."""
    o = PalMarked(a, policy=policy, **kw)
    d, tot, steps = palmach2.pf.walk(np.asarray(a, np.int64), o, 'incr')
    ink = 0 if policy == 'none' else o.pl.ink
    return dict(reads=tot, ink=ink, delta=d, end=(0 if d is not None else 1),
                steps=steps, consults=o.consults, skips=o.mark_skips,
                writes=o.mark_writes, seeds=o.seeds, hits=o.hits)


def run_base(a):
    d, cost, steps, o = palmach2.walk_mach2(a)
    return dict(reads=cost, ink=0, delta=d, end=(0 if d is not None else 1),
                steps=steps, consults=0, skips=0, writes=0, seeds=o.seeds,
                hits=o.hits)


POLICIES = ('none', 'cert', 'cert_naive', 'run')
