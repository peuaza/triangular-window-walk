# -*- coding: utf-8 -*-
"""q0: fast exact b1=0 (largest-concavity) walk, verified against the references.

Model recap (Lemma A of REPORT.md):
    with A_j = mean(a[L..L+j-1]) and m = r-L+1,
    b2 = 1  <=>  A_m > A_{floor(m/2)}          (second block mean > first block mean)
    b2 = 0  <=>  A_m <= A_{floor(m/2)}.

b1:=0 rule: an *epoch* starts at (L,R); scan m = M, M-1, ..., 2 (M = R-L+1).
The largest m with b2=0 is m*.  If none: output L/(n-1).  If m*=M: double cut.
Else: left cut.  An epoch therefore costs M - m* + 1 b2 tests (or M-1 if it meets).

All windows in one epoch start at the same L, so one prefix-sum array of a
answers every test in O(1); an epoch is O(M), a whole run O(n^2).
"""
import numpy as np


def prefix(a):
    P = np.zeros(len(a) + 1, dtype=np.int64)
    np.cumsum(a, out=P[1:])
    return P


def epoch(P, L, R):
    """Return (m*, cost) for the epoch at window [L,R]; m*=None means 'meet'."""
    M = R - L + 1
    if M < 2:
        return None, 0
    m = np.arange(2, M + 1, dtype=np.int64)
    k = m // 2
    Sm = P[L + m] - P[L]
    Sk = P[L + k] - P[L]
    # b2 = 0  <=>  Sm/m <= Sk/k  <=>  Sk*m >= Sm*k
    conc = Sk * m >= Sm * k
    idx = np.flatnonzero(conc)
    if idx.size == 0:
        return None, M - 1
    mstar = int(m[idx[-1]])
    return mstar, M - mstar + 1


def walk_b0_fast(a, want_diag=False):
    """Exact b1=0 walk.  Returns (delta or None, n_tests, diagnostics)."""
    n = len(a)
    P = prefix(a)
    L, R = 0, n - 1
    tests = 0
    diag = []
    seen = set()
    while True:
        M = R - L + 1
        if L >= R:
            return L / (n - 1), tests, diag
        if M < 3:
            return None, tests, diag
        if (L, R) in seen:
            return None, tests, diag
        seen.add((L, R))
        mstar, cost = epoch(P, L, R)
        tests += cost
        if want_diag:
            diag.append((L, R, M, mstar, cost))
        if mstar is None:
            return L / (n - 1), tests, diag
        if mstar == M:
            L, R = L + 1, R - 1
        else:
            L += 1


def walk_b0_ref(a, n):
    """Reference implementation, same control flow as p2_b0.walk_b0."""
    L, R = 0, n - 1
    ell, r = 0, n - 1
    steps = 0
    seen = set()
    while True:
        if ell >= r:
            return ell / (n - 1), steps
        if R - L + 1 < 3:
            return None, steps
        steps += 1
        if steps > 20 * n * n:
            return None, steps
        if (ell, r) in seen:
            return None, steps
        seen.add((ell, r))
        z = a[ell:r + 1]
        m = len(z)
        k = max(1, m // 2)
        b2 = int(z.sum()) * k > int(z[:k].sum()) * m
        dy, dx = (0, -1) if b2 else (-1, 0)
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


# ------------------------------------------------------------------ ensemble
def sample_tapes(n, ntapes, seed):
    """Uniform over ensemble tapes: two binary columns with ceil(n/2) ones each.

    Equivalent to: uniform over a in {0,1,2}^n with sum(a) = 2*ceil(n/2)
    (n even: #0=#2; n odd: #0=#2+1).
    """
    h = n // 2 if n % 2 == 0 else (n + 1) // 2
    rng = np.random.default_rng(seed)
    for _ in range(ntapes):
        B = np.zeros((n, 2), dtype=np.int16)
        for c in range(2):
            B[rng.choice(n, h, replace=False), c] = 1
        yield B.sum(1).astype(np.int64)


if __name__ == '__main__':
    import sys
    sys.path.insert(0, '../probe_ologn')
    from p0_tapes import enum_tapes

    # exhaustive agreement with the reference walk (b1 forced to 0)
    for n in (6, 7, 8, 9, 10, 11, 12):
        tapes = list(enum_tapes(n))
        bad = 0
        for a in tapes:
            d1, t1, _ = walk_b0_fast(a)
            d2, t2 = walk_b0_ref(a, n)
            ok = ((d1 is None) == (d2 is None)) and (d1 is None or abs(d1 - d2) < 1e-12)
            if not ok:
                bad += 1
        print('n=%2d tapes=%6d  agreement %d/%d' % (n, len(tapes), len(tapes) - bad, len(tapes)))
