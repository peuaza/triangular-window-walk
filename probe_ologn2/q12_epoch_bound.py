# -*- coding: utf-8 -*-
"""q12: the epoch lower bound  P(cost_0 >= c M) >= p  for the b1=0 concavity scan.

Deliverables
  A  exact finite-M covariance of the tape prefix sums  (Lemma 1, exact formula)
  B  FCLT sanity: Var of the dyadic contrast == theta * M/(M-1) exactly (Lemma 2)
  C  discrete first-epoch law: P(cost=1), P(cost>=cM), conditional means
  D  Brownian counterpart  V(theta) = W(theta) - 2 W(theta/2):
       Var V = theta,  Cov V(theta),V(theta') = (2 th< - th>)_+   (Lemma 2 limit)
       p(c) = P(sup_{[1-c,1]} V < 0)   and   p_cond(c) = 2 p(c)
  E  constants: simulated (c,p) vs the rigorous Borell+Dudley bound
  F  ensemble-independence check (sigma^2 = 1/2 vs 2/3)

Conventions (identical to q0_walk.py):
  A_j = mean(a[L..L+j-1]) = 1 + X_j/j,  k = floor(m/2),
  b2 = 1  <=>  A_m > A_k  <=>  X_m - 2 X_k > 0     (even m; odd m: > X_k/k)
  cost = M - m* + 1 where m* = largest m in [2,M] with b2 = 0 at m
         (cost = M-1 if no such m).
"""
import sys
import numpy as np

from q0_walk import prefix, epoch, sample_tapes


# ---------------------------------------------------------------- A / B
def cov_check(n=1024, ntapes=400, seed=2026):
    """Exact Cov(X_i,X_j) = s^2 i (M-j)/(M-1) and Var of the dyadic contrast."""
    print('\n== A. covariance of prefix sums (exact formula) ==')
    th = np.array([0.125, 0.25, 0.375, 0.5, 0.75, 1.0])
    D = []
    V = []
    for a in sample_tapes(n, ntapes, seed):
        X = np.concatenate(([0], np.cumsum(a - 1))).astype(float)
        s2 = np.mean((a - 1.0) ** 2)
        # normalised contrast at thetai and thetai/2 (exact integer indices)
        g = []
        for t in th:
            j = int(round(t * n))
            k = j // 2
            g.append((X[j] - 2 * X[k]) / np.sqrt(n * s2))
        V.append(g)
        D.append([(X[int(round(t1 * n))] * X[int(round(t2 * n))]) / (n * s2)
                  for t1 in th for t2 in th])
    V = np.array(V)
    cv = np.cov(V.T)
    print('  empirical Corr of D(theta) = (X_t - 2X_{t/2})/sqrt(n s^2)')
    print('        theta :', ' '.join('%7.3f' % t for t in th))
    pred = np.array([[2 * min(a, b) - max(a, b) for b in th] for a in th])
    pred = np.where(pred > 0, pred, 0.0) / np.sqrt(np.outer(th, th))
    for i in range(len(th)):
        print('   rho_emp %7.3f:' % th[i], ' '.join('%7.3f' % x for x in cv[i] / np.sqrt(cv[i, i] * cv.diagonal())))
    print('   rho_pred (2 a^-b^)_+/sqrt(ab):')
    for i in range(len(th)):
        print('              ', ' '.join('%7.3f' % x for x in pred[i]))
    sd = np.sqrt(np.diag(cv))
    print('   sd_emp :', ' '.join('%7.3f' % x for x in sd))
    print('   sqrt(th):', ' '.join('%7.3f' % x for x in np.sqrt(th)))
    print('   ratio sd^2/theta (should be M/(M-1)=%.5f)' % (n / (n - 1.0)))
    print('              ', ' '.join('%7.5f' % (x ** 2 / t) for x, t in zip(sd, th)))


# ---------------------------------------------------------------- C
def first_epoch(ns=(64, 128, 256, 512, 1024, 2048, 4096, 8192),
                ntapes=800, seed=2026,
                cs=(0.05, 0.1, 0.15, 0.2, 0.3, 0.4)):
    """Discrete first-epoch law from the full region [0,n-1]."""
    print('\n== C. first-epoch cost, full region, %d tapes/n, seed %d ==' % (ntapes, seed))
    print('%6s %9s %9s %8s %8s %8s %8s   %s'
          % ('n', 'E[cost]/M', 'se', 'P(c=1)', 'P(c>=2)', 'E[c/M|c>=2]',
             'P(c>=.1M)', ' '.join('P>=%.2f' % c for c in cs)))
    rows = {}
    for n in ns:
        nt = ntapes if n <= 4096 else max(200, ntapes // 4)
        cost, M = [], []
        for a in sample_tapes(n, nt, seed):
            P = prefix(a)
            ms, c = epoch(P, 0, n - 1)
            cost.append(c)
            M.append(n)
        cost = np.array(cost, float)
        M = np.array(M, float)
        r = cost / M
        big = cost >= 2
        ps = [np.mean(cost >= c * M) for c in cs]
        rows[n] = (r.mean(), ps, np.mean(big), np.mean(r[big]) if big.any() else np.nan)
        print('%6d %9.4f %9.4f %8.3f %8.3f %8.4f %8.3f   %s'
              % (n, r.mean(), r.std(ddof=1) / np.sqrt(len(r)), np.mean(cost == 1),
                 np.mean(big), r[big].mean() if big.any() else np.nan,
                 np.mean(cost >= 0.1 * M),
                 ' '.join('%.2f' % p for p in ps)))
        sys.stdout.flush()
    return rows


def ensemble_B(n, ntapes, seed):
    """Uniform over words with sum(a) = 2 ceil(n/2)  (different sigma^2)."""
    from math import lgamma
    h = n // 2 if n % 2 == 0 else (n + 1) // 2
    tot = 2 * h
    rng = np.random.default_rng(seed + 1)
    js = np.arange(0, tot // 2 + 1)
    lw = np.array([lgamma(n + 1) - 2 * lgamma(j + 1) - lgamma(n - 2 * j + 1) for j in js])
    w = np.exp(lw - lw.max())
    w /= w.sum()
    for _ in range(ntapes):
        j = int(rng.choice(js, p=w))
        a = np.concatenate([np.zeros(j, int), np.ones(n - 2 * j, int), np.full(j, 2, int)])
        rng.shuffle(a)
        yield a


# ---------------------------------------------------------------- D
def brownian(N=4096, reps=40000, seed=11, cs=(0.02, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6)):
    """V(theta) = W(theta) - 2 W(theta/2) at theta_j = j/N.

    Var V = theta, Cov = (2 th< - th>)_+ ;  cost/M -> 1 - theta^*.
    """
    rng = np.random.default_rng(seed)
    cs = np.array(cs)
    ncc = len(cs)
    hits = np.zeros(ncc)          # path with sup_{[1-c,1]} V < 0
    v1neg = 0
    m1 = np.zeros(ncc)
    tail = []                      # E[1-theta^*]
    cond_tail = []                 # E[1-theta^* | V(1)<0]
    chunk = 500
    done = 0
    supdev = {c: [] for c in (0.05, 0.1, 0.2, 0.3)}
    while done < reps:
        B = min(chunk, reps - done)
        dW = rng.normal(0.0, 1.0 / np.sqrt(2 * N), size=(B, 2 * N))
        Wf = np.concatenate([np.zeros((B, 1)), np.cumsum(dW, axis=1)], axis=1)
        V = Wf[:, 2 * np.arange(N + 1)] - 2.0 * Wf[:, np.arange(N + 1)]
        th = np.arange(N + 1) / float(N)
        neg = V <= 0.0
        # largest theta with V >= 0
        idx = np.where(~neg, np.arange(N + 1)[None, :], -1).max(axis=1)
        val = 1.0 - idx / float(N)
        val[idx < 0] = 1.0
        tail.append(val)
        isneg = V[:, N] < 0
        v1neg += int(isneg.sum())
        cond_tail.append(val[isneg])
        # sup over the top slice
        for i, c in enumerate(cs):
            lo = int(np.ceil((1.0 - c) * N))
            mx = V[:, lo:].max(axis=1)
            hits[i] += int(np.sum(mx < 0.0))
        for c in supdev:
            lo = int(np.ceil((1.0 - c) * N))
            supdev[c].append(np.abs(V[:, lo:] - V[:, N:N + 1]).max(axis=1))
        done += B
    tail = np.concatenate(tail)
    cond_tail = np.concatenate(cond_tail)
    hits /= reps
    print('\n== D. Brownian V(theta)=W(theta)-2W(theta/2), N=%d, %d paths ==' % (N, reps))
    print('  E[1-theta^*] = %.4f   (report q5 value 0.1365)' % tail.mean())
    print('  median %.4f   P(1-theta^*=0) = %.3f   P(V(1)<0) = %.4f'
          % (np.median(tail), np.mean(tail == 0.0), v1neg / reps))
    print('  E[1-theta^* | V(1)<0] = %.4f' % cond_tail.mean())
    print('\n   c     p(c)=P(sup<0)   p_cond=2p(c)   [c*p]     E[sup|V-V(1)|]  (Dudley 2.17 sqrt c)')
    for i, c in enumerate(cs):
        d = np.concatenate(supdev[c]) if c in supdev else None
        print('  %.2f    %8.4f       %8.4f     %7.4f    %s'
              % (c, hits[i], 2 * hits[i], c * hits[i],
                 ('%8.3f        %8.3f' % (d.mean(), 1.2533 * np.sqrt(3.0 * c))) if d is not None else ''))
    return hits, tail


def brownian_grid_check():
    """grid-refinement: is the discrete sup close to the continuous one?"""
    print('\n== D2. grid sensitivity of p(c) ==')
    for N in (1024, 2048, 4096, 8192):
        rng = np.random.default_rng(5)
        reps = 6000
        out = []
        for c in (0.1, 0.5):
            hits = 0
            for _ in range(reps):
                dW = rng.normal(0.0, 1.0 / np.sqrt(2 * N), size=2 * N)
                Wf = np.concatenate(([0.0], np.cumsum(dW)))
                V = Wf[2 * np.arange(N + 1)] - 2 * Wf[np.arange(N + 1)]
                lo = int(np.ceil((1 - c) * N))
                hits += V[lo:].max() < 0
            out.append(hits / reps)
        print('  N=%5d  p(0.1)=%.4f  p(0.5)=%.4f' % (N, out[0], out[1]))


# ---------------------------------------------------------------- E
def rigor_p(cs=(0.02, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4)):
    """Rigorous lower bound for p(c) = P(sup_{[1-c,1]} V < 0).

    Condition on V(1) = -a (a>0).  Then V(theta) = -a(2theta-1)_+ + Z(theta) with Z
    centred Gaussian, Var Z(theta|a) = theta - (2theta-1)^2 <= 9/16, and
    E|Z(th)-Z(th')|^2 <= 3|th-th'|, so E sup_{[1-c,1]} Z <= 2.17 sqrt(c) (Dudley:
    D=sqrt(3c), N(eps)=3c/eps^2, int_0^D sqrt(2 log(D/eps)) d eps = 1.2533 D).
    sup V <= sup Z - a(1-2c)  =>  P(sup V >= 0 | a) <= exp(-((a(1-2c)-m)_+)^2/(2*9/16)).
    p(c) >= 1 - E_a[ that ].
    """
    from math import sqrt, exp
    import numpy as np
    s2 = 9.0 / 16.0
    print('\n== E. rigorous constants (Borell/Tsirelson + Dudley) ==')
    print('    c     m_c=EsupZ   p_rig     c*p_rig')
    best = (0, 0, 0)
    a = np.linspace(0, 8, 4001)
    phi = np.exp(-a * a / 2.0) * np.sqrt(2.0 / np.pi)      # half-normal density
    for c in cs:
        m = 1.2533 * sqrt(3.0 * c)
        k = 1.0 - 2.0 * c
        bad = np.exp(-np.maximum(a * k - m, 0.0) ** 2 / (2 * s2))
        Ebad = np.trapezoid(bad * phi, a)
        p = max(0.0, 1.0 - Ebad)
        if c * p > best[2]:
            best = (c, p, c * p)
        print('  %.2f    %.3f      %.4f    %.4f' % (c, m, p, c * p))
    print('  best on this grid: c=%.2f  p=%.3f  c*p=%.4f' % best)
    return best


# ---------------------------------------------------------------- F
def ensemble_independence(ns=(512, 2048), ntapes=600):
    print('\n== F. ensemble independence of E[cost_0]/M (sigma^2 = 1/2 vs 2/3) ==')
    for n in ns:
        out = []
        for gen, lab in ((lambda: sample_tapes(n, ntapes, 2026), 'two balanced cols'),
                         (lambda: ensemble_B(n, ntapes, 2026), 'uniform sum=2h')):
            c = []
            s2 = []
            for a in gen():
                P = prefix(a)
                ms, cc = epoch(P, 0, n - 1)
                c.append(cc / n)
                s2.append(np.mean((a - 1.0) ** 2))
            c = np.array(c)
            out.append((lab, c.mean(), c.std(ddof=1) / np.sqrt(len(c)), np.mean(s2)))
        for lab, m, se, s2 in out:
            print('  n=%5d %-18s sigma^2=%.4f  E[cost]/M = %.4f +- %.4f'
                  % (n, lab, s2, m, se))


# ---------------------------------------------------------------- G
def validate(n=4096, ntapes=3000, seed=4242,
             cs=(0.05, 0.1, 0.2, 0.3, 0.5),
             pp=None):
    """Discrete P(cost >= cM) vs the Brownian p(c); conditional version too."""
    print('\n== G. discrete vs Brownian, n=%d, %d tapes ==' % (n, ntapes))
    cost, M = [], []
    sumM, nep = [], []
    for a in sample_tapes(n, ntapes, seed):
        P = prefix(a)
        ms, c = epoch(P, 0, n - 1)
        cost.append(c)
        M.append(n)
    cost = np.array(cost, float)
    M = np.array(M, float)
    top = cost >= 2
    print('   c    P(cost>=cM) +- se   2P  (Wald)   p(c) Brownian   ratio   E[c/M | c>=cM]')
    for i, c in enumerate(cs):
        ind = cost >= c * M
        p = ind.mean()
        se = np.sqrt(p * (1 - p) / len(ind))
        pb = pp[i] if pp is not None else np.nan
        print('  %.2f    %.4f +- %.4f      %.4f        %.4f      %6.3f    %.4f'
              % (c, p, se, 2 * p, pb, p / pb if pb else np.nan,
                 cost[ind].mean() / n if ind.any() else np.nan))
    print('  P(cost=1)=%.4f  P(cost>=2)=%.4f  E[c/M|c>=2]=%.4f  E[c/M]=%.4f'
          % (np.mean(cost == 1), top.mean(), cost[top].mean() / n, cost.mean() / n))
    for c in (0.1, 0.2):
        ind = cost >= c * M
        print('  P(cost>=%.1fM | cost>=2) = %.4f   [Brownian 2p(c) = %.4f]'
              % (c, ind.sum() / top.sum(), 2 * (pp[cs.index(c)] if pp is not None else np.nan)))


def window_geometry(ns=(512, 1024), ntapes=150):
    """Deterministic facts: L_i = i, M_i = n - i - delta_i, sum M_i >= n^2/4."""
    from q0_walk import walk_b0_fast
    print('\n== H. epoch geometry (deterministic identities) ==')
    for n in ns:
        sm, ne, ls, bad = [], [], [], 0
        for a in sample_tapes(n, ntapes, 909):
            d, T, diag = walk_b0_fast(a, want_diag=True)
            if not diag:
                continue
            Ls = np.array([x[0] for x in diag])
            Ms = np.array([x[2] for x in diag])
            ms = np.array([x[3] if x[3] is not None else x[2] + 1 for x in diag])
            if not np.array_equal(Ls, np.arange(len(diag))):
                bad += 1
            ne.append(len(diag))
            sm.append(float(Ms.sum()))
            ls.append(float(np.sum(Ms - ms)))
        print('  n=%5d  L_i == i violations: %d   epochs/n: min %.3f mean %.3f   '
              'sum M/n^2: min %.3f mean %.3f'
              % (n, bad, min(ne) / n, np.mean(ne) / n,
                 min(sm) / n ** 2, np.mean(sm) / n ** 2))


# ---------------------------------------------------------------- I
def exact_variance_identity(Ms=(8, 10, 12)):
    """Exhaustive: for every word with #0=#2, E[D_M(theta)^2] = s^2 * theta * M/(M-1)
    with s^2 = (1/M) sum b_i^2 (per-word normalisation).  Also checks the exact-zero
    cross-covariance for theta' >= 2 theta."""
    import sys
    sys.path.insert(0, '../probe_ologn')
    from p0_tapes import enum_tapes
    print('\n== I. exact variance identity (exhaustive over the ensemble) ==')
    for M in Ms:
        th = [0.25, 0.5, 0.75, 1.0]
        acc = np.zeros((len(th), len(th)))
        cnt = 0
        for a in enum_tapes(M):
            b = a - 1.0
            s2 = np.mean(b * b)
            if s2 == 0:
                continue
            X = np.concatenate(([0.0], np.cumsum(b)))
            g = np.array([X[int(round(t * M))] - 2 * X[int(round(t * M)) // 2] for t in th])
            acc += np.outer(g, g) / (M * s2)
            cnt += 1
        acc /= cnt
        print('  M=%2d  words=%6d   E[D(th)D(th\')]/(M s^2) observed:' % (M, cnt))
        for i, t in enumerate(th):
            print('        th=%.2f :' % t, ' '.join('%8.4f' % x for x in acc[i]))
        print('        exact diagonal theta*M/(M-1):',
              ' '.join('%8.4f' % (t * M / (M - 1.0)) for t in th))
        print('        limit (2 th< - th>)_+       :',
              ' '.join('%8.4f' % max(2 * min(t, u) - max(t, u), 0) for t in th for u in [t]))


if __name__ == '__main__':
    which = sys.argv[1] if len(sys.argv) > 1 else 'all'
    if which in ('A', 'all'):
        cov_check()
    if which in ('C', 'all'):
        first_epoch()
    if which in ('D', 'all'):
        brownian()
        brownian_grid_check()
    if which in ('E', 'all'):
        rigor_p()
    if which in ('F', 'all'):
        ensemble_independence()
    if which in ('G', 'all'):
        # p(c) from the Brownian run in section D (N=4096, 40000 paths)
        validate(pp=(0.3810, 0.3282, 0.2504, 0.1914, 0.0913))
        window_geometry()
    if which in ('I', 'all'):
        exact_variance_identity()
