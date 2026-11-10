# -*- coding: utf-8 -*-
"""sl_machine.py -- a LITERAL in-class machine for the threshold-run rule

    b1(L,m) = [ the longest run of 2s inside the window [L,L+m-1] exceeds t(n) ]
    b2(L,m) = [ mean(window) > mean(its first floor(m/2) cells) ]

Model (the billing criterion of INCLASS_COSTS.md section 1): the machine holds the
read-only tape plus O(log n) bits.  Every cell read is charged.  No array indexed by
position, no precomputed table.

Two exact algorithms, selected per tape by the one number rho = longest 2-run of the
whole tape (computed in the opening pass, O(log n) bits):

  rho <= t   b1 == 0 on EVERY window (a run longer than t cannot fit in any window).
             The walk is therefore the b1=0 walk, and the b0 shortcut answers one
             epoch with one top test plus a short upward scan of window lengths.
             Per-epoch cost is O(1 + first zero of the row), ~1.15 tests per cell.

  rho > t    no tape-wide bound is available: b1 is decided by scanning the window
             left to right and stopping at the first run of t+1 twos (or at the end).

b2 is maintained, never searched: S = sum(window), Q = sum of its first k cells,
plus the same two registers for the region [L,R] so that a cut resets the window in
O(1).  The bookkeeping is the one validated in probe_ologn2/rv_machine.py.
"""
import sys
import numpy as np

sys.path.insert(0, '../probe_ologn2')
from rs_common import ref_walk, sample_tapes  # noqa: E402


class _M(object):
    def __init__(self, a, t, store=False, K0=32):
        self.A = [int(x) for x in a]
        self.n = len(a)
        self.t = t
        self.reads = 0
        self.cache = {}
        self.maxcache = 0
        self.scans = 0
        self.scan_reads = 0
        self.store = store
        self.K0 = K0
        self.runs = []          # maximal 2-runs longer than t, capped at K0
        self.overflow = False

    # ---- reads
    def raw(self, i):
        self.reads += 1
        return self.A[i]

    def get(self, i):
        if i in self.cache:
            return self.cache[i]
        v = self.raw(i)
        self.cache[i] = v
        if len(self.cache) > self.maxcache:
            self.maxcache = len(self.cache)
        return v

    def prune(self, keep):
        for k in [q for q in self.cache if q not in keep]:
            del self.cache[k]
        if len(self.cache) > 16:
            raise RuntimeError('register file grew to %d' % len(self.cache))

    # ---- b1 by scan (only called when m > t and rho > t)
    def b1_scan(self, L, m):
        t = self.t
        self.scans += 1
        run = 0
        for i in range(L, L + m):
            if self.raw(i) == 2:
                self.scan_reads += 1
                run += 1
                if run > t:
                    return True
            else:
                self.scan_reads += 1
                run = 0
        return False

    def b1(self, L, m):
        """Exact b1.  With store=True the precomputed list of maximal 2-runs longer
        than t answers the overlap query in O(K0) register operations (0 reads).
        A run [s,e] puts b1=1 iff |[s,e] cap [L,L+m-1]| > t."""
        if m <= self.t:
            return False
        if not self.store:
            return self.b1_scan(L, m)
        r = L + m - 1
        for (s, e) in self.runs:
            lo = s if s > L else L
            hi = e if e < r else r
            if hi - lo + 1 > self.t:
                return True
        if self.overflow:
            return self.b1_scan(L, m)
        return False

    # ---- opening pass: rho, total sum, first-half sum, long-run witness list
    def open(self):
        n = self.n
        rho = 0
        cur = 0
        S = 0
        start = 0
        for i in range(n):
            v = self.raw(i)
            S += v
            if v == 2:
                if cur == 0:
                    start = i
                cur += 1
                if cur > rho:
                    rho = cur
            else:
                if cur > self.t and self.store:
                    if len(self.runs) < self.K0:
                        self.runs.append((start, i - 1))
                    else:
                        self.overflow = True
                cur = 0
        if cur > self.t and self.store:
            if len(self.runs) < self.K0:
                self.runs.append((start, n - 1))
            else:
                self.overflow = True
        k = max(1, n // 2)
        Q = 0
        for i in range(k):
            Q += self.raw(i)
        return rho, S, Q, k


def solve(a, t, want=None, store=False, K0=32):
    """Returns dict(delta, reads, steps, mode, rho, maxcache, nlong, overflow).
    store=False: b1 by window scan on every uncertified test.
    store=True : b1 from the precomputed list of maximal 2-runs longer than t."""
    m = _M(a, t, store=store, K0=K0)
    n = m.n
    rho, Sreg, Qreg, kreg = m.open()
    L, R = 0, n - 1
    M = n
    reads0 = m.reads
    steps = 0
    mode = 'shortcut' if rho <= t else 'literal'

    if mode == 'shortcut':
        # ---- b0 shortcut, registers (Sreg,Qreg,kreg) for the region
        while True:
            if L >= R:
                return dict(delta=L / (n - 1), reads=m.reads, steps=steps,
                            mode=mode, rho=rho, maxcache=m.maxcache, init=reads0, scans=m.scans, scan_reads=m.scan_reads, nlong=len(m.runs), overflow=m.overflow)
            if M < 3:
                return dict(delta=None, reads=m.reads, steps=steps, mode=mode,
                            rho=rho, maxcache=m.maxcache, init=reads0, scans=m.scans, scan_reads=m.scan_reads, nlong=len(m.runs), overflow=m.overflow)
            steps += 1
            if Sreg * kreg > Qreg * M:          # top bit 1 -> descend inside
                found = False
                # running sums for the upward scan of window lengths
                A = m.raw(L)                    # sum(L)
                kprev = 1
                B = A                           # first half of m=1 is cell L
                zero = None
                for mm in range(2, M):
                    A += m.raw(L + mm - 1)
                    kk = max(1, mm // 2)
                    if kk > kprev:
                        B += m.raw(L + kk - 1)
                        kprev = kk
                    if not (A * kk > B * mm):
                        zero = mm
                        break
                if zero is not None:            # left cut
                    vL = m.get(L)
                    kold = kreg
                    L += 1
                    M = R - L + 1
                    kreg = max(1, M // 2)
                    Sreg -= vL
                    if kreg == kold:
                        Qreg = Qreg - vL + m.get(L - 1 + kold)
                    else:
                        Qreg -= vL
                    m.prune({L, R, L + kreg - 1, L + kreg})
                    continue
                return dict(delta=L / (n - 1), reads=m.reads, steps=steps,
                            mode=mode, rho=rho, maxcache=m.maxcache, init=reads0, scans=m.scans, scan_reads=m.scan_reads, nlong=len(m.runs), overflow=m.overflow)
            # top bit 0 -> double cut
            vL = m.get(L)
            vR = m.get(R)
            kold = kreg
            Sreg -= vL + vR
            L, R = L + 1, R - 1
            M = R - L + 1
            kreg = max(1, M // 2)
            if kreg == kold:
                Qreg = Qreg - vL + m.get(L)
            else:
                Qreg -= vL
            m.prune({L, R, L + kreg - 1, L + kreg})

    # ---- literal walk, b2 maintained, b1 scanned only when m > t
    ell, r = L, R
    Sm, Qm = Sreg, Qreg
    km = kreg
    seen = set()
    cap = 200 * n * n
    # region cut helper, exact rv_machine bookkeeping
    def region_step(kind):
        nonlocal L, R, M, kreg, Sreg, Qreg, ell, r, Sm, Qm, km
        vL = m.get(L)
        vR = m.get(R)
        if kind == 'double':
            Sreg = Sreg - vL - vR
            kold = kreg
            L, R = L + 1, R - 1
            M = R - L + 1
            kreg = max(1, M // 2)
            if kreg == kold:
                Qreg = Qreg - vL + m.get(L)
            else:
                Qreg = Qreg - vL
        elif kind == 'left':
            Sreg -= vL
            kold = kreg
            L += 1
            M = R - L + 1
            kreg = max(1, M // 2)
            if kreg == kold:
                Qreg = Qreg - vL + m.get(L - 1 + kold)
            else:
                Qreg = Qreg - vL
        else:
            Sreg -= vR
            kold = kreg
            R -= 1
            M = R - L + 1
            kreg = max(1, M // 2)
            if kreg != kold:
                Qreg = Qreg - m.get(L + kold - 1)
        ell, r = L, R
        Sm, Qm, km = Sreg, Qreg, kreg
        seen.clear()
        m.prune({L, R, ell, r, ell + km - 1, ell + km, L + kreg - 1, L + kreg})

    while True:
        if ell >= r:
            return dict(delta=ell / (n - 1), reads=m.reads, steps=steps,
                        mode=mode, rho=rho, maxcache=m.maxcache, init=reads0, scans=m.scans, scan_reads=m.scan_reads, nlong=len(m.runs), overflow=m.overflow)
        if M < 3:
            return dict(delta=None, reads=m.reads, steps=steps, mode=mode,
                        rho=rho, maxcache=m.maxcache, init=reads0, scans=m.scans, scan_reads=m.scan_reads, nlong=len(m.runs), overflow=m.overflow)
        steps += 1
        if steps > cap:
            return dict(delta=None, reads=m.reads, steps=steps, mode=mode,
                        rho=rho, maxcache=m.maxcache, init=reads0, cap=True,
                        scans=m.scans, scan_reads=m.scan_reads, nlong=len(m.runs), overflow=m.overflow)
        key = (ell, r)
        if key in seen:
            return dict(delta=None, reads=m.reads, steps=steps, mode=mode,
                        rho=rho, maxcache=m.maxcache, init=reads0, cycle=True,
                        scans=m.scans, scan_reads=m.scan_reads, nlong=len(m.runs), overflow=m.overflow)
        seen.add(key)
        mm = r - ell + 1
        b2v = Sm * km > Qm * mm
        b1v = False if mm <= t else m.b1(ell, mm)
        if b2v:
            dy, dx = (1, 0) if b1v else (0, -1)
        else:
            dy, dx = (0, 1) if b1v else (-1, 0)
        if ell + dy < L or r + dx > R:
            if ell == L and r == R:
                region_step('double')
            elif ell + dy < L:
                region_step('left')
            else:
                region_step('right')
            continue
        if ell + dy >= r + dx:
            return dict(delta=(ell + dy) / (n - 1), reads=m.reads, steps=steps,
                        mode=mode, rho=rho, maxcache=m.maxcache, init=reads0, scans=m.scans, scan_reads=m.scan_reads, nlong=len(m.runs), overflow=m.overflow)
        if r + dx < L or ell + dy > R:
            region_step('double')
            continue
        # ---- move with the validated S,Q bookkeeping
        mold = mm
        if dy == 1:
            x = m.get(ell + 1)
            v = m.get(ell)
            Sm -= v
            mm -= 1
            knew = max(1, mm // 2)
            if mold % 2 == 0:
                Qm -= v
            else:
                Qm = Qm - v + m.get(ell + km)
            ell += 1
            km = knew
        elif dx == -1:
            y = m.get(r - 1)
            v = m.get(r)
            Sm -= v
            mm -= 1
            knew = max(1, mm // 2)
            if mold % 2 == 0 and mold > 2:
                Qm -= m.get(ell + km - 1)
            r -= 1
            km = knew
        elif dx == 1:
            z = m.get(r + 1)
            v = m.get(r)
            Sm += z
            mm += 1
            knew = max(1, mm // 2)
            if mold % 2 == 1:
                Qm += m.get(ell + km)
            r += 1
            km = knew
        else:
            w = m.get(ell - 1)
            v = m.get(ell)
            Sm += w
            mm += 1
            knew = max(1, mm // 2)
            if mold % 2 == 0:
                Qm = Qm + w - m.get(ell + km - 1)
            else:
                Qm += w
            ell -= 1
            km = knew
        m.prune({L, R, ell, r, ell + km - 1, ell + km, L + kreg - 1, L + kreg})


def t_of(n, c=-1):
    import math
    return int(math.floor(math.log2(n))) + c


if __name__ == '__main__':
    import time
    cmd = sys.argv[1] if len(sys.argv) > 1 else 'bench'
    if cmd == 'bench':
        ns = [int(x) for x in sys.argv[2].split(',')]
        N = int(sys.argv[3]) if len(sys.argv) > 3 else 200
        seed = int(sys.argv[4]) if len(sys.argv) > 4 else 2026
        c = int(sys.argv[5]) if len(sys.argv) > 5 else -1
        print('%5s %4s %5s %12s %8s %12s %8s %10s %8s %8s' %
              ('n', 't', 'rare', 'mean_rd', 'e/n', 'max_rd', 'g/n',
               'exc_mean', 'exc/n2', 'maxcache'))
        for n in ns:
            t = t_of(n, c)
            rows = []
            for a in sample_tapes(n, N, seed):
                rows.append(solve(a, t))
            rd = np.array([r['reads'] for r in rows], float)
            rho = np.array([r['rho'] for r in rows])
            exc = rd[rho > t]
            print('%5d %4d %5d %12.1f %8.2f %12.1f %8.1f %10.1f %8.3f %8d' %
                  (n, t, len(exc), rd.mean(), rd.mean() / n, rd.max(),
                   rd.max() / n, exc.mean() if len(exc) else float('nan'),
                   (exc.mean() / n ** 2) if len(exc) else float('nan'),
                   max(r['maxcache'] for r in rows)), flush=True)
