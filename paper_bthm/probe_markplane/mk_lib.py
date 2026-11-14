# -*- coding: utf-8 -*-
"""mk_lib -- the one-bit mark plane: model, tape, charging.

MODEL (declared once, applied to both rules)
--------------------------------------------
tape    the same T_n as probe_triple_nat/tr_lib: two independent uniform h-subsets
        of [n], h = n//2, rng = default_rng(5000 + n + 100003 i), tapes i = 0..299.
        A cell carries its data a[i] in {0,1,2}, read-only, and ONE writable bit,
        the mark m[i].  Nothing else about the data plane changes.
read    rd(i) returns (a[i], m[i]) and costs 1 read.  val(i) costs 1 read and
        returns a[i].  cmp(i,j) costs 1 read and returns a[i]==a[j] -- this is the
        pal machine's declared primitive ("one comparison = one read", the charge
        rule printed in every tr_lib table for pal); for the pcr machine a read is
        one cell, and the two are never mixed inside one machine.
write   wr(i,b) costs 1 ink, 0 reads.  The mark is overwritable and may be set to
        0 or 1 at any time.  It is the machine's own bookkeeping; it is NOT a
        record of whether the cell was read (nothing in either machine consults a
        mark to decide whether to read).
trim    trim(lo,hi) removes those cells and their marks.  It costs no read: the
        machines' own accounting already charges the cells they give up (pcr's
        INC "1 per cell taken in or given up", and the re-span term).  Nothing is
        read after a trim: `alive` is asserted on every access.
state   the machine's work tape stays O(log n) bits.  A variant may hold counters
        or one certificate, exactly as the baseline machine does; nothing
        per-cell may be carried in state.
report  reads and ink are counted and reported separately.
"""
import numpy as np


class Trimmed(Exception):
    pass


class Plane(object):
    """One tape: data (read-only), one mark bit per cell, two counters."""

    __slots__ = ('n', 'a', 'm', 'alive', 'reads', 'ink', 'nread_cell')

    def __init__(self, a):
        self.a = np.asarray(a, np.int8)
        self.n = int(self.a.shape[0])
        self.m = np.zeros(self.n, np.int8)
        self.alive = np.ones(self.n, bool)
        self.reads = 0
        self.ink = 0
        self.nread_cell = 0          # distinct cells ever read

    # ---- data plane -----------------------------------------------------
    def rd(self, i):
        if not self.alive[i]:
            raise Trimmed(i)
        self.reads += 1
        return int(self.a[i]), int(self.m[i])

    def val(self, i):
        if not self.alive[i]:
            raise Trimmed(i)
        self.reads += 1
        return int(self.a[i])

    def cmp(self, i, j):
        """pal's charged primitive: one comparison of two cells = one read."""
        if not self.alive[i] or not self.alive[j]:
            raise Trimmed((i, j))
        self.reads += 1
        return int(self.a[i]) == int(self.a[j])

    # ---- mark plane -----------------------------------------------------
    def wr(self, i, b):
        if not self.alive[i]:
            raise Trimmed(i)
        self.ink += 1
        self.m[i] = 1 if b else 0

    def mark(self, i):
        if not self.alive[i]:
            raise Trimmed(i)
        return int(self.m[i])

    def trim(self, lo, hi):
        if hi < lo:
            return
        self.alive[lo:hi + 1] = False
        self.m[lo:hi + 1] = 0


class Tape(object):
    """The baseline's tape: data only, no plane.  Used to re-run the baseline
    machines under the same call shape (they never touch a mark)."""

    __slots__ = ('n', 'a', 'reads')

    def __init__(self, a):
        self.a = np.asarray(a, np.int8)
        self.n = int(self.a.shape[0])
        self.reads = 0

    def val(self, i):
        self.reads += 1
        return int(self.a[i])

    def cmp(self, i, j):
        self.reads += 1
        return int(self.a[i]) == int(self.a[j])
