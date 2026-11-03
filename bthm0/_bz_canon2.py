# -*- coding: utf-8 -*-
"""规范规则（起点双切、边单切、成环终止）的增量模拟。

原版 _bz_inc.py 是现行规则（成环双切继续）。本文件只改切割与成环两处。


Two machines, both deterministic, of the canonical tape-and-window form:

  M_code  = the published implementation (_onetable.py walk): b1 = np.corrcoef(..) > 0.
            At an exactly-zero autocorrelation the sign is decided by float roundoff.
  M_exact = the mathematical rule: b1 = [autocorrelation > 0] in exact arithmetic
            (tie -> False).  This is the literal reading of "corr > 0".

Claim under test: M_exact's meeting position is computable with O(1) cell reads per
step and O(log n) bits of state, i.e. O(n) reads per tape -- not Theta(n^2).

Why: a cut fires only when the moving header already sits on the interval boundary
(ny < L forces y = L; nx > R forces x = R).  The new window is then exactly the old
window with one end cell deleted, so every accumulator updates in O(1) reads.
No rebuild of length m is ever needed.

Tape cells are held scaled by 2 (values 0,1,2) so every comparison is exact integer
arithmetic:  mean(w) > mean(head) <=> S*k > Sh*m, and cov > 0 <=> np*Rr - A*B > 0.

Usage:
    python _bz_inc.py check <n> <tapes> [seed]   # vs brute-force exact + vs published code
    python _bz_inc.py bench <n> <tapes> [seed]   # read-count distribution
"""
import sys
import time

import numpy as np


# ---------------------------------------------------------------- published code
def walk_ref(a, n):
    """Verbatim copy of the published walk (_onetable.py).  Float arithmetic."""
    L, R = 0, n - 1
    ell, r = 0, n - 1
    steps = 0
    while True:
        if ell >= r:
            return ell / (n - 1)
        if R - L + 1 < 3:
            return None
        steps += 1
        if steps > 20 * n:
            return None
        m = r - ell + 1
        k = max(1, m // 2)
        z = a[ell:r + 1]
        b2 = z.mean() > z[:k].mean()
        if m == 1:
            b1 = a[ell] > 0.5
        elif m == 2:
            b1 = a[r] > a[ell]
        else:
            u, v = z[:-1], z[1:]
            b1 = (np.corrcoef(u, v)[0, 1] > 0) if (u.std() > 1e-12 and v.std() > 1e-12) else False
        o = ("B" if b1 else "D") if b2 else ("C" if b1 else "A")
        ny, nx = {"A": (ell - 1, r), "B": (ell + 1, r), "C": (ell, r + 1), "D": (ell, r - 1)}[o]
        if ny < L or nx > R:
            if ny < L:
                R = r - 1
            else:
                L = ell + 1
            ell, r = L, R
            continue
        if ny >= nx:
            return ny / (n - 1)
        if nx < L or ny > R:
            L, R = L + 1, R - 1
            ell, r = L, R
            continue
        ell, r = ny, nx


# ---------------------------------------------------------------- exact brute force
def walk_exact_bf(c, n):
    """Independent oracle for M_exact: recompute every window statistic from
    scratch in exact integer arithmetic, O(m) reads per step."""
    L, R = 0, n - 1
    y, x = 0, n - 1
    steps = 0
    while True:
        if y >= x:
            return y / (n - 1)
        if R - L + 1 < 3:
            return None
        steps += 1
        if steps > 20 * n:
            return None
        m = x - y + 1
        k = max(1, m // 2)
        z = c[y:x + 1].astype(np.int64)
        b2 = int(z.sum()) * k > int(z[:k].sum()) * m
        if m == 1:
            b1 = int(z[0]) > 1
        elif m == 2:
            b1 = int(z[1]) > int(z[0])
        else:
            np_ = m - 1
            S = int(z.sum())
            Q = int((z * z).sum())
            Rr = int((z[:-1] * z[1:]).sum())
            A = S - int(z[-1])
            B = S - int(z[0])
            C_ = Q - int(z[-1]) ** 2
            D_ = Q - int(z[0]) ** 2
            if np_ * C_ - A * A <= 0 or np_ * D_ - B * B <= 0:
                b1 = False
            else:
                b1 = (np_ * Rr - A * B) > 0
        o = ("B" if b1 else "D") if b2 else ("C" if b1 else "A")
        dy, dx = {"A": (-1, 0), "B": (1, 0), "C": (0, 1), "D": (0, -1)}[o]
        ny, nx = y + dy, x + dx
        if ny < L or nx > R:
            if ny < L:
                R = x - 1
            else:
                L = y + 1
            y, x = L, R
            continue
        if ny >= nx:
            return ny / (n - 1)
        if nx < L or ny > R:
            L, R = L + 1, R - 1
            y, x = L, R
            continue
        y, x = ny, nx


# ---------------------------------------------------------------- instrumented tape
class Tap(object):
    __slots__ = ('c', 'reads')

    def __init__(self, c):
        self.c = c
        self.reads = 0

    def __getitem__(self, i):
        self.reads += 1
        return self.c[i]


# ---------------------------------------------------------------- O(1)/step state
class Win(object):
    """State: y, x (window); S, Q, Rr over it; cy=c[y], cx=c[x]; head [hy,hx]
    with Sh.  All integer, O(log n) bits each."""

    def __init__(self, T, y, x):
        self.T = T
        self.y, self.x = y, x
        S = Q = Rr = 0
        for i in range(y, x + 1):
            v = T[i]
            S += v
            Q += v * v
        for i in range(y, x):
            Rr += T[i] * T[i + 1]
        self.S, self.Q, self.Rr = S, Q, Rr
        self.cy, self.cx = T[y], T[x]
        m = x - y + 1
        k = max(1, m // 2)
        self.hy, self.hx = y, y + k - 1
        self.Sh = 0
        for i in range(self.hy, self.hx + 1):
            self.Sh += T[i]

    def add_left(self):
        v = self.T[self.y - 1]
        self.Rr += v * self.cy
        self.S += v
        self.Q += v * v
        self.y -= 1
        self.cy = v

    def del_left(self):
        v = self.cy
        nv = self.T[self.y + 1]
        self.Rr -= v * nv
        self.S -= v
        self.Q -= v * v
        self.y += 1
        self.cy = nv

    def add_right(self):
        v = self.T[self.x + 1]
        self.Rr += self.cx * v
        self.S += v
        self.Q += v * v
        self.x += 1
        self.cx = v

    def del_right(self):
        v = self.cx
        pv = self.T[self.x - 1]
        self.Rr -= pv * v
        self.S -= v
        self.Q -= v * v
        self.x -= 1
        self.cx = pv

    def sync_head(self):
        m = self.x - self.y + 1
        k = max(1, m // 2)
        ty, tx = self.y, self.y + k - 1
        while self.hy > ty:
            self.hy -= 1
            self.Sh += self.T[self.hy]
        while self.hy < ty:
            self.Sh -= self.T[self.hy]
            self.hy += 1
        while self.hx < tx:
            self.hx += 1
            self.Sh += self.T[self.hx]
        while self.hx > tx:
            self.Sh -= self.T[self.hx]
            self.hx -= 1

    def move(self, dy, dx):
        if dy < 0:
            self.add_left()
        elif dy > 0:
            self.del_left()
        if dx < 0:
            self.del_right()
        elif dx > 0:
            self.add_right()
        self.sync_head()


# ---------------------------------------------------------------- incremental walk
def walk_inc(T, n, verify=None):
    """Exact simulation of M_exact in O(1) reads/step.  verify (integer tape) turns
    on a from-scratch check of the state and of both decision bits at every step."""
    L, R = 0, n - 1
    w = Win(T, 0, n - 1)
    steps = 0
    cuts = 0
    ties = 0
    seen = set()
    while True:
        if w.y >= w.x:
            out = w.y / (n - 1)
            break
        if R - L + 1 < 3:
            out = None
            break
        steps += 1
        if (w.y, w.x) in seen:
            out = None
            break
        seen.add((w.y, w.x))
        if steps > 20 * n:
            out = None
            break
        m = w.x - w.y + 1
        k = max(1, m // 2)
        b2 = (w.S * k > w.Sh * m)
        if m == 1:
            b1 = w.cy > 1
        elif m == 2:
            b1 = w.cx > w.cy
        else:
            np_ = m - 1
            A = w.S - w.cx
            B = w.S - w.cy
            C_ = w.Q - w.cx * w.cx
            D_ = w.Q - w.cy * w.cy
            if np_ * C_ - A * A <= 0 or np_ * D_ - B * B <= 0:
                b1 = False
            else:
                d = np_ * w.Rr - A * B
                if d == 0:
                    ties += 1
                b1 = d > 0
        if verify is not None:
            z = verify[w.y:w.x + 1]
            kk = max(1, m // 2)
            assert w.S == int(z.sum()), ('S', w.y, w.x)
            assert w.Sh == int(z[:kk].sum()), ('Sh', w.y, w.x, m)
            assert b2 == (int(z.sum()) * kk > int(z[:kk].sum()) * m), ('b2', w.y, w.x)
            if m >= 3:
                zz = z.astype(np.int64)
                assert (w.S, w.Q, w.Rr) == (int(zz.sum()), int((zz * zz).sum()),
                                            int((zz[:-1] * zz[1:]).sum())), ('stat', w.y, w.x, m)
        o = ("B" if b1 else "D") if b2 else ("C" if b1 else "A")
        dy, dx = {"A": (-1, 0), "B": (1, 0), "C": (0, 1), "D": (0, -1)}[o]
        ny, nx = w.y + dy, w.x + dx
        if ny < L or nx > R:
            cuts += 1
            if (w.y, w.x) == (L, R):
                # 起点即出界 => 双切
                L, R = L + 1, R - 1
                w.del_left(); w.del_right()
                w.sync_head()
                seen.clear()
                continue
            if ny < L:
                assert w.y == L, ('cutL', w.y, L)
                R = w.x - 1
                w.del_right()
            else:
                assert w.x == R, ('cutR', w.x, R)
                L = w.y + 1
                w.del_left()
            w.sync_head()
            seen.clear()
            if verify is not None:
                assert (w.y, w.x) == (L, R), ('cutwin', w.y, w.x, L, R)
            continue
        if ny >= nx:
            out = ny / (n - 1)
            break
        if nx < L or ny > R:
            L, R = L + 1, R - 1
            w = Win(T, L, R)
            continue
        w.move(dy, dx)
    return out, T.reads, steps, cuts, ties


# ---------------------------------------------------------------- tape generator
def make_tape(rng, n):
    h = n // 2
    B = np.zeros((n, 2), dtype=np.int8)
    for c in range(2):
        B[rng.choice(n, h, replace=False), c] = 1
    return B.mean(1)


def r2(pred, true):
    pred = np.asarray(pred, float)
    true = np.asarray(true, float)
    return float(1 - ((true - pred) ** 2).sum() / ((true - true.mean()) ** 2).sum())


def main():
    task = sys.argv[1] if len(sys.argv) > 1 else 'check'
    n = int(sys.argv[2]) if len(sys.argv) > 2 else 256
    NW = int(sys.argv[3]) if len(sys.argv) > 3 else 200
    seed = int(sys.argv[4]) if len(sys.argv) > 4 else 2026
    rng = np.random.default_rng(seed)
    budget = n * n / np.log(n)
    print('n=%d tapes=%d seed=%d   budget n^2/ln n = %.0f  (n^2 = %d)' % (n, NW, seed, budget, n * n))
    rows = []
    t0 = time.time()
    n_exact = n_code = 0
    both = []
    for _ in range(NW):
        a = make_tape(rng, n)
        c = np.round(a * 2).astype(np.int32)
        T = Tap(c)
        d_inc, reads, steps, cuts, ties = walk_inc(T, n, verify=c if task == 'check' else None)
        if task == 'check':
            d_bf = walk_exact_bf(c, n)
            assert (d_bf is None) == (d_inc is None), ('bf-none', d_bf, d_inc)
            if d_bf is not None:
                assert abs(d_bf - d_inc) < 1e-12, ('bf', d_bf, d_inc)
        d_ref = walk_ref(a, n)
        if d_inc is not None:
            n_exact += 1
        if d_ref is not None:
            n_code += 1
        if d_inc is not None and d_ref is not None:
            both.append((d_ref, d_inc))
        rows.append((d_ref is not None, d_inc is not None, reads, steps, cuts, ties))
    rows = np.array(rows, dtype=float)
    meet = rows[:, 1] > 0
    R = rows[meet]
    print('meeting: M_exact %d/%d (%.3f)   M_code %d/%d (%.3f)'
          % (n_exact, NW, n_exact / NW, n_code, NW, n_code / NW))
    print('%-10s %10s %10s %10s' % ('statistic', 'mean', 'median', 'max'))
    for i, nm in enumerate(('reads', 'steps', 'cuts', 'ties')):
        v = R[:, i + 2]
        print('%-10s %10.1f %10.1f %10.0f' % (nm, v.mean(), np.median(v), v.max()))
    print('reads/n = %.2f   reads/n^2 = %.4f   reads/(n^2/ln n) = %.3f   max/budget = %.3f'
          % (R[:, 2].mean() / n, R[:, 2].mean() / n / n, R[:, 2].mean() / budget, R[:, 2].max() / budget))
    print('meeting tapes with reads < budget: %.3f' % (R[:, 2] < budget).mean())
    if both:
        t = np.array([b[0] for b in both])
        p = np.array([b[1] for b in both])
        print('M_exact vs M_code, %d commonly-meeting tapes: exact agree %.4f, R2 = %.5f, max|d| = %.4f'
              % (len(both), (np.abs(t - p) < 1e-12).mean(), r2(p, t), np.abs(t - p).max()))
    print('elapsed %.1fs' % (time.time() - t0))


if __name__ == '__main__':
    main()
