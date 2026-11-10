# -*- coding: utf-8 -*-
"""sl_screen.py -- screen O(1)-maintainable b1 bits for a cheap, non-degenerate walk.

The walk is the literal rs_common control flow.  b1 candidates are evaluated from
prefix structures (out-of-class on purpose: this is a *screen* for trajectory
length and output spread, not a cost measurement).
"""
import sys
import numpy as np
from collections import Counter

sys.path.insert(0, '../probe_ologn2')
from rs_common import sample_tapes, manacher, longest_run_ref  # noqa: E402


class P:
    def __init__(self, a):
        self.a = [int(x) for x in a]
        n = len(a)
        self.n = n
        self.P = [0] * (n + 1)
        self.P2 = [0] * (n + 1)   # count of 2s
        self.PP = [0] * (n + 1)   # adjacent products
        for i in range(n):
            self.P[i + 1] = self.P[i] + self.a[i]
            self.P2[i + 1] = self.P2[i] + (1 if self.a[i] == 2 else 0)
            if i:
                self.PP[i + 1] = self.PP[i] + self.a[i - 1] * self.a[i]
            else:
                self.PP[i + 1] = 0
        self.ch = [0] * (n + 1)
        for i in range(1, n):
            self.ch[i + 1] = self.ch[i] + (1 if self.a[i] != self.a[i - 1] else 0)
        # runs of 2
        rs, re_ = [], []
        i = 0
        while i < n:
            if self.a[i] == 2:
                j = i
                while j + 1 < n and self.a[j + 1] == 2:
                    j += 1
                rs.append(i); re_.append(j); i = j + 1
            else:
                i += 1
        self.rs, self.re = rs, re_

    def S(self, L, m):
        return self.P[L + m] - self.P[L]

    def cnt2(self, L, m):
        return self.P2[L + m] - self.P2[L]

    def changes(self, L, m):
        if m <= 1:
            return 0
        return self.ch[L + m - 1] - self.ch[L]

    def Rr(self, L, m):
        if m <= 1:
            return 0
        return self.PP[L + m - 1] - self.PP[L + 1]

    def lrun(self, L, m):
        # longest run of 2s inside [L, L+m-1]
        if m <= 0:
            return 0
        best = 0
        for s, e in zip(self.rs, self.re):
            if e < L:
                continue
            if s > L + m - 1:
                break
            v = min(e, L + m - 1) - max(s, L) + 1
            if v > best:
                best = v
        return best


def b2(P, L, m):
    k = max(1, m // 2)
    return P.S(L, m) * k > P.S(L, k) * m


def make_b1(P, rule, t=0):
    if rule == 'b0':
        return lambda L, m: False
    if rule == 'b1':
        return lambda L, m: True
    if rule == 'canon':
        def f(L, m):
            if m == 1:
                return False
            S = P.S(L, m)
            A = S - P.a[L + m - 1]
            B = S - P.a[L]
            return (m - 1) * P.Rr(L, m) - A * B > 0
        return f
    if rule == 'c2last':
        return lambda L, m: P.a[L + m - 1] == 2
    if rule == 'c2first':
        return lambda L, m: P.a[L] == 2
    if rule == 'mean':
        return lambda L, m: P.S(L, m) > m
    if rule == 'cnt2':
        return lambda L, m: 4 * P.cnt2(L, m) > m
    if rule == 'runs':
        return lambda L, m: 2 * P.changes(L, m) > m
    if rule == 'parity':
        return lambda L, m: (P.changes(L, m) & 1) == 1
    if rule == 'runthr':
        return lambda L, m: P.lrun(L, m) > t
    if rule == 'runlog':
        return lambda L, m: P.lrun(L, m) > max(0, int(np.floor(np.log2(m))) - 1)
    if rule == 'longrun':
        return lambda L, m: 4 * P.lrun(L, m) > m
    if rule == 'ge2':
        return lambda L, m: P.S(L, m) >= m + 1
    if rule == 'le0':
        return lambda L, m: P.S(L, m) <= m - 1
    if rule == 'all1':
        return lambda L, m: P.S(L, m) == m
    if rule == 'ends_eq':
        return lambda L, m: P.a[L] == P.a[L + m - 1]
    if rule == 'ends2':
        return lambda L, m: P.a[L] == 2 and P.a[L + m - 1] == 2
    if rule == 'par':
        return lambda L, m: (P.S(L, m) & 1) == 1
    if rule == 'lenpar':
        return lambda L, m: (m & 1) == 0
    if rule == 'cnt2lo':
        return lambda L, m: 4 * P.cnt2(L, m) < m
    if rule == 'second':
        return lambda L, m: P.S(L, m) - P.S(L, max(1, m // 2)) > m - max(1, m // 2)
    raise ValueError(rule)


def walk(a, rule, t=0, cap_mult=40):
    n = len(a)
    Q = P(a)
    b1f = make_b1(Q, rule, t)
    L, R = 0, n - 1
    ell, r = L, R
    seen = set()
    steps = 0
    cap = cap_mult * n * n
    while True:
        if ell >= r:
            return ell / (n - 1), steps
        if R - L + 1 < 3:
            return None, steps
        steps += 1
        if steps > cap:
            return None, steps
        key = (ell, r)
        if key in seen:
            return None, steps
        seen.add(key)
        m = r - ell + 1
        b2v = b2(Q, ell, m)
        b1v = b1f(ell, m)
        if b2v:
            dy, dx = (1, 0) if b1v else (0, -1)
        else:
            dy, dx = (0, 1) if b1v else (-1, 0)
        if ell + dy < L or r + dx > R:
            if ell == L and r == R:
                L, R = L + 1, R - 1
            elif ell + dy < L:
                L += 1
            else:
                R -= 1
            ell, r = L, R
            seen.clear()
            continue
        if ell + dy >= r + dx:
            return (ell + dy) / (n - 1), steps
        if r + dx < L or ell + dy > R:
            L, R = L + 1, R - 1
            ell, r = L, R
            seen.clear()
            continue
        ell, r = ell + dy, r + dx


def stat(n, rule, N, seed, t=0):
    st = []
    outs = []
    none = 0
    for a in sample_tapes(n, N, seed):
        d, s = walk(a, rule, t)
        st.append(s)
        if d is None:
            none += 1
        else:
            outs.append(round(d, 9))
    c = Counter(outs)
    tot = len(outs)
    H = -sum((v / tot) * np.log2(v / tot) for v in c.values()) if tot else 0.0
    return dict(steps=float(np.mean(st)), maxsteps=int(max(st)), Pnone=none / N,
                distinct=len(c), H=H, modal=max(c.values()) / tot if tot else 0.0)


if __name__ == '__main__':
    rules = sys.argv[1].split(',')
    ns = [int(x) for x in sys.argv[2].split(',')]
    N = int(sys.argv[3]) if len(sys.argv) > 3 else 100
    seed = int(sys.argv[4]) if len(sys.argv) > 4 else 2026
    print('%10s %5s %10s %10s %8s %9s %6s %6s' %
          ('rule', 'n', 'steps', 'maxsteps', 'P(none)', 'distinct', 'H', 'modal'))
    for rule in rules:
        for n in ns:
            t = int(np.floor(np.log2(n))) if rule == 'runthr' else 0
            s = stat(n, rule, N, seed, t)
            print('%10s %5d %10.1f %10d %8.3f %9d %6.2f %6.3f' %
                  (rule, n, s['steps'], s['maxsteps'], s['Pnone'],
                   s['distinct'], s['H'], s['modal']), flush=True)
