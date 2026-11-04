# -*- coding: utf-8 -*-
"""q13_meet.py -- does the b1=0 run ever MEET (scan finds no concavity at all)?

The question this file settles numerically, and the quantities it defines:

  window [L,R], M=R-L+1.  C = {m in [2,M] : A_floor(m/2) >= A_m}   (concavities)
  meet  <=>  C = empty   <=>  cost = M-1   <=>  the run outputs L.
  f(M)  := P(C = empty) for a window whose content is a fresh uniform
           arrangement  (iid b in {-1,0,1}, or the two-column ensemble).

Proved here (see REPORT_MEET.md):
  T1  P(run meets at some epoch i < (n-1)/4)  <=  sum_{i<(n-1)/4} f(n-2i)
      -- uses L_i = i, M_i >= n-2i, and the fact that the test set of a window
         is an initial segment in m, so C=empty survives truncation of R.
  T2  f(2h) = E[ 1{E_h(W)} q(W) ]  with E_h the inner event and q the exact
      barrier probability of the second half  (identity, checked below).
  T3  E_M = empty when b_1 = 1, and E_M with b_1 = 0 forces X_{2^j} >= 2^j - 1
      (near-explosion): P <= 3^{-M/2}(2M+1).

Modes:  exact | rec | run | dip | mc <M> <N>
"""
import sys
import numpy as np

# ---------------------------------------------------------------- exact counts
def f_exact_brute(M, chunk=3000000):
    """Exact count of b in {-1,0,1}^M with C=empty.  criterion m*X_k < k*X_m."""
    m = np.arange(2, M + 1); k = m // 2
    tot = 3 ** M; idx = 0; hits = 0
    while idx < tot:
        n = min(chunk, tot - idx)
        ids = np.arange(idx, idx + n)
        B = np.empty((n, M), dtype=np.int8); x = ids.copy()
        for j in range(M):
            B[:, j] = (x % 3) - 1; x //= 3
        X = np.zeros((n, M + 1), dtype=np.int16)
        np.cumsum(B, axis=1, out=X[:, 1:])
        ok = np.ones(n, bool)
        for a, b in zip(m, k):
            ok &= (a * X[:, b] < b * X[:, a])
        hits += int(ok.sum()); idx += n
    return hits / tot, hits, tot

def E_ok(X, M):
    ok = np.ones(X.shape[0], bool)
    for m in range(2, M + 1):
        j = m // 2
        ok &= (j * X[:, m] > m * X[:, j])
    return ok

def recursion(h, batch=150000):
    """Exact f(2h) by enumerating the 3^h first halves and computing q exactly.

    Returns (f2h, P(E_h), mean q | E_h, max q, min over E_h of min_m X_m).
    q(W) = P( V_t >= th(t) for t=1..h ),  th(t) = floor((m*X_j - j*X_h)/j)+1,
    m=h+t, j=m//2  [from j*(X_h+V_t) > m*X_j].
    """
    tot = 3 ** h; s_q = 0.0; nE = 0; mx = -1e9; mn = 10 ** 9; sumq = 0.0
    for start in range(0, tot, batch):
        n = min(batch, tot - start)
        ids = np.arange(start, start + n)
        B = np.empty((n, h), dtype=np.int8); x = ids.copy()
        for j in range(h):
            B[:, j] = (x % 3) - 1; x //= 3
        X = np.zeros((n, h + 1), dtype=np.int64)
        np.cumsum(B, axis=1, out=X[:, 1:])
        ok = E_ok(X, h)
        if not ok.any():
            continue
        Xs = X[ok]
        th = np.zeros((Xs.shape[0], h + 1), dtype=np.int64)
        for t in range(1, h + 1):
            m = h + t; j = m // 2
            th[:, t] = np.floor_divide(m * Xs[:, j] - j * Xs[:, h], j) + 1
        dp = np.zeros((Xs.shape[0], 2 * h + 1)); dp[:, h] = 1.0
        vv = np.arange(-h, h + 1)
        for t in range(1, h + 1):
            new = np.zeros_like(dp)
            new[:, :-1] = dp[:, 1:] / 3.0
            new[:, 1:] += dp[:, :-1] / 3.0
            new[:, :] += dp[:, :] / 3.0
            new *= (vv[None, :] >= th[:, t][:, None])
            dp = new
        q = dp.sum(1)
        s_q += q.sum(); sumq += q.sum(); nE += int(q.size)
        mx = max(mx, float(q.max())); mn = min(mn, int(Xs[:, 1:].min()))
    return s_q / tot, nE / tot, sumq / max(nE, 1), mx, mn

# ------------------------------------------------------------------ run itself
def walk_diag(a):
    """b1=0 walk with the diagnostics Theorem T1 needs."""
    n = len(a)
    P = np.zeros(n + 1, dtype=np.int64); np.cumsum(a, out=P[1:])
    L, R = 0, n - 1; i = 0; ncuts = 0; bad = 0; meets = []
    while True:
        if L >= R: return dict(end='L>=R', meets=meets, bad=bad, i=i)
        M = R - L + 1
        if M < 3: return dict(end='M<3', meets=meets, bad=bad, i=i)
        if L != i: bad += 1
        m = np.arange(2, M + 1, dtype=np.int64); k = m // 2
        Sm = P[L + m] - P[L]; Sk = P[L + k] - P[L]
        conc = Sk * m >= Sm * k
        idx = np.flatnonzero(conc)
        if idx.size == 0:
            meets.append((i, M)); return dict(end='meet', meets=meets, bad=bad, i=i)
        mstar = int(m[idx[-1]])
        if mstar == M: L, R = L + 1, R - 1; ncuts += 1
        else: L += 1
        i += 1

def tapes(n, ntape, seed, kind='2col'):
    rng = np.random.default_rng(seed)
    if kind == 'iid':
        for _ in range(ntape):
            yield rng.integers(0, 3, n)
    else:
        h = n // 2 if n % 2 == 0 else (n + 1) // 2
        for _ in range(ntape):
            B = np.zeros((n, 2), dtype=np.int16)
            for c in range(2):
                B[rng.choice(n, h, replace=False), c] = 1
            yield B.sum(1).astype(np.int64)

# ---------------------------------------------------------------------- driver
def main():
    mode = sys.argv[1] if len(sys.argv) > 1 else 'exact'
    if mode == 'exact':
        print('Exact f(M) by full enumeration of 3^M words (iid ensemble)')
        for M in range(2, 17):
            p, h, t = f_exact_brute(M)
            print('  M=%2d  f=%.6e   (%d/%d)' % (M, p, h, t)); sys.stdout.flush()
    elif mode == 'rec':
        print('Exact f(2h) from the recursion  f(2h)=E[1{E_h} q]  (iid ensemble)')
        print('   h    M        f(M)      P(E_h)    E[q|E_h]  max q  min_m X on E_h')
        for h in range(4, 15):
            p, pE, mq, mxq, mnX = recursion(h)
            print('  %2d  %3d  %.6e  %.4e  %8.4f  %.3f  %d'
                  % (h, 2 * h, p, pE, mq, mxq, mnX)); sys.stdout.flush()
    elif mode == 'run':
        print('Run-level check of T1: L_i=i, M_i>=n-2i, and early meets')
        for n in (64, 128, 256, 512, 1024):
            nt = max(100, int(2e6 / n))
            bad = 0; early = 0; anymeet = 0; Ms = []
            for a in tapes(n, nt, 4242, '2col'):
                d = walk_diag(a); bad += d['bad']
                if d['meets']:
                    anymeet += 1
                    if d['meets'][0][0] < (n - 1) / 4: early += 1
            print('  n=%4d tapes=%5d  P(meet)=%.4f  P(meet before (n-1)/4)=%.1e'
                  '  L!=i: %d' % (n, nt, anymeet / nt, early / nt, bad)); sys.stdout.flush()
    elif mode == 'dip':
        print('Depth of the walk on E_M (iid ensemble, exact)')
        m = None
        for M in (10, 12, 14, 16):
            m = np.arange(2, M + 1); k = m // 2
            tot = 3 ** M; idx = 0; hist = {}
            while idx < tot:
                n = min(2000000, tot - idx)
                ids = np.arange(idx, idx + n)
                B = np.empty((n, M), dtype=np.int8); x = ids.copy()
                for j in range(M):
                    B[:, j] = (x % 3) - 1; x //= 3
                X = np.zeros((n, M + 1), dtype=np.int16)
                np.cumsum(B, axis=1, out=X[:, 1:])
                ok = np.ones(n, bool)
                for a, b in zip(m, k): ok &= (a * X[:, b] < b * X[:, a])
                mn = X[ok][:, 1:].min(1)
                for v, c in zip(*np.unique(mn, return_counts=True)):
                    hist[int(v)] = hist.get(int(v), 0) + int(c)
                idx += n
            tot_hits = sum(hist.values())
            print('  M=%2d hits=%7d  min histogram %s' % (M, tot_hits, sorted(hist.items())))
            print('        P(min<=-2)=%.4f P(min<=-3)=%.4f P(min<=-4)=%.4f'
                  % (sum(c for v, c in hist.items() if v <= -2) / tot_hits,
                     sum(c for v, c in hist.items() if v <= -3) / tot_hits,
                     sum(c for v, c in hist.items() if v <= -4) / tot_hits))
            sys.stdout.flush()
    elif mode == 'barrier':
        print('Barrier statistics on E_h (iid, exact); checks the three "easy" bounds of REPORT_MEET §8')
        print('   h    nE   E[dip] P(dip>=3)  E[Lam] maxLam P(Lam>=5)  E[q]  maxq | beta2/q  maxbeta2')
        for h in (8, 10, 12, 14):
            tot = 3 ** h; recs = []
            for start in range(0, tot, 120000):
                n = min(120000, tot - start)
                ids = np.arange(start, start + n)
                B = np.empty((n, h), dtype=np.int8); x = ids.copy()
                for j in range(h):
                    B[:, j] = (x % 3) - 1; x //= 3
                X = np.zeros((n, h + 1), dtype=np.int64)
                np.cumsum(B, axis=1, out=X[:, 1:])
                ok = E_ok(X, h)
                if not ok.any():
                    continue
                Xs = X[ok]; N = Xs.shape[0]
                th = np.zeros((N, h + 1), dtype=np.int64)
                for t in range(1, h + 1):
                    m = h + t; j = m // 2
                    th[:, t] = np.floor_divide(m * Xs[:, j] - j * Xs[:, h], j) + 1
                OFF = h; W = 2 * h + 1; vv = np.arange(-h, h + 1)
                dp = np.zeros((N, W)); dp[:, OFF] = 1.0
                for t in range(1, h + 1):
                    new = np.zeros_like(dp)
                    new[:, :-1] = dp[:, 1:] / 3.0
                    new[:, 1:] += dp[:, :-1] / 3.0
                    new[:, :] += dp[:, :] / 3.0
                    new *= (vv[None, :] >= th[:, t][:, None])
                    dp = new
                q = dp.sum(1)
                dp = np.zeros((N, W)); dp[:, OFF] = 1.0
                for t in range(1, h + 1):
                    new = np.zeros_like(dp)
                    new[:, :-1] = dp[:, 1:] / 3.0
                    new[:, 1:] += dp[:, :-1] / 3.0
                    new[:, :] += dp[:, :] / 3.0
                    if t == h // 2:
                        new *= (vv[None, :] >= th[:, h // 2][:, None])
                    new *= (vv[None, :] >= th[:, t][:, None])
                    dp = new
                b2 = dp.sum(1)
                lam = -th[:, 1:].min(1)
                dip = -Xs[:, 1:].min(1)
                recs.append(np.stack([dip, lam, q, b2], 1))
            R = np.concatenate(recs, 0)
            dip, lam, q, b2 = R[:, 0], R[:, 1], R[:, 2], R[:, 3]
            print('  %2d %6d  %6.2f  %7.3f  %6.2f %6d  %8.3f  %.3f  %.3f | %5.2f  %8.3f'
                  % (h, len(R), dip.mean(), (dip >= 3).mean(), lam.mean(), lam.max(),
                     (lam >= 5).mean(), q.mean(), q.max(), b2.sum() / q.sum(), b2.max()))
            sys.stdout.flush()
    elif mode == 'mc':
        M = int(sys.argv[2]); N = int(sys.argv[3])
        rng = np.random.default_rng(98765)
        m = np.arange(2, M + 1); k = m // 2
        hits = 0; tot = 0
        while tot < N:
            n = min(1000000, N - tot)
            B = rng.integers(-1, 2, size=(n, M), dtype=np.int8)
            X = np.zeros((n, M + 1), dtype=np.int16)
            np.cumsum(B, axis=1, out=X[:, 1:])
            ok = np.ones(n, bool)
            for a, b in zip(m, k): ok &= (a * X[:, b] < b * X[:, a])
            hits += int(ok.sum()); tot += n
        print('M=%d N=%d f=%.4e hits=%d relse=%.2f%%'
              % (M, N, hits / tot, hits, 100 * hits ** -0.5))

if __name__ == '__main__':
    main()
