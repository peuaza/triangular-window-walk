# -*- coding: utf-8 -*-
"""Route 5 probes: the dyadic-contrast structure of the b2 rule.

Nothing here modifies any .tex.  Everything is checked against the *same* walk
as `_arr_rules.walk` with b1 forced to 0 (`walk_b0` of p2_b0.py), which
Theorem B of REPORT.md identifies with the largest-concavity scan on all
306,965 ensemble tapes for n = 6..13.

Checks
  V1  Haar / half-difference normal form, both parities (exhaustive windows).
  V2  shift invariance: b2 and hence delta depend only on z = a-1.
  V3  epoch decomposition: one descending scan per anchor, anchors strictly
      increasing, halt iff no negative contrast at the anchor.
  V4  the epoch answer is a one-pass function of the window (3 counters).
  V5  model-R accounting: sum_e (l_e - 2r) for r = 0,1,2,4,...
  V6  lossiness: #distinct delta values vs #distinct first-anchor sign vectors.

Usage:
  $P lb_route5.py verify 6 7 8 9 10 11 12
  $P lb_route5.py measure 6 8 10 12 14 16 24 32 48 64 96 128
  $P lb_route5.py lossy 6 8 10 12
"""
import sys
import math

import numpy as np

sys.path.insert(0, '..')
from p0_tapes import enum_tapes


# --------------------------------------------------------------------------- #
# the reference walk (b1 forced to 0) -- copied verbatim from p2_b0.py
# --------------------------------------------------------------------------- #
def bit_b2(z):
    m = len(z)
    k = max(1, m // 2)
    return int(z.sum()) * k > int(z[:k].sum()) * m


def walk_b0(a, n):
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
        b2 = bit_b2(a[ell:r + 1])
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


# --------------------------------------------------------------------------- #
# the epoch form:  delta = the first anchor L at which no contrast is negative
# --------------------------------------------------------------------------- #
def epoch_walk(z, keep=False):
    """delta of the b1=0 rule, computed from the epoch decomposition.

    Returns (delta, epochs) with epochs = [(L, R, M, mstar_or_0), ...];
    mstar = 0 means "no negative contrast at this anchor" (the halt).
    """
    n = len(z)
    Z = np.concatenate(([0], np.cumsum(z)))
    L, R = 0, n - 1
    E = []
    while True:
        M = R - L + 1
        if M < 2:
            # walk_b0 tests `ell >= r` before `R-L+1 < 3`, so M=1 returns L
            return L / (n - 1), E
        if M < 3:
            return None, E
        ms = np.arange(2, M + 1)
        k = ms // 2
        ZL = Z[L]
        neg = (Z[L + ms] - ZL) * k <= (Z[L + k] - ZL) * ms
        if not neg.any():
            E.append((L, R, M, 0))
            return L / (n - 1), E
        mstar = int(ms[neg][-1])
        E.append((L, R, M, mstar))
        if mstar == M:
            L, R = L + 1, R - 1
        else:
            L += 1


# --------------------------------------------------------------------------- #
# V4: the epoch answer as a one-pass function of the window with 3 counters
# --------------------------------------------------------------------------- #
def epoch_onepass(z, L, R):
    """m* = max{m in [2,M] : b2(L,m)=0}, or 0.  One left-to-right pass,
    counters T (whole prefix sum), P (first-half prefix sum), best (m*)."""
    M = R - L + 1
    T = 0
    P = 0
    half = 0
    best = 0
    for m in range(1, M + 1):
        T += int(z[L + m - 1])
        h = m // 2
        while half < h:
            P += int(z[L + half])
            half += 1
        if m >= 2:
            if (T - P) * h <= P * (m - h):     # b2 = 0
                best = m
    return best


# --------------------------------------------------------------------------- #
def enum_z(n):
    for a in enum_tapes(n):
        yield a - 1


def rand_z(n, nw, seed=2026):
    rng = np.random.default_rng(seed)
    h = n // 2
    for _ in range(nw):
        B = np.zeros((n, 2), dtype=np.int16)
        for c in range(2):
            B[rng.choice(n, h, replace=False), c] = 1
        yield B.sum(1).astype(np.int64) - 1


# --------------------------------------------------------------------------- #
def cmd_verify(ns):
    print('== V1  Haar / half-difference normal form ==')
    bad1 = bad2 = bad3 = tot = 0
    rng = np.random.default_rng(7)
    for trial in range(200000):
        n = int(rng.integers(2, 60))
        a = rng.integers(0, 3, size=n)
        L = int(rng.integers(0, n - 1))
        r = int(rng.integers(L + 1, n))
        m = r - L + 1
        k = m // 2
        z = a - 1
        tot += 1
        b2 = bit_b2(a[L:r + 1])
        # even/odd closed forms on z: mean(second block) > mean(first block)
        Pz = int(z[L:L + k].sum())
        Qz = int(z[L + k:L + m].sum())
        pred = int(Qz * k > Pz * (m - k))
        bad1 += (pred != b2)
        if m % 2 == 0:
            # Haar coefficient at scale k: sum(second half) - sum(first half) > 0
            h = int(z[L + k:L + 2 * k].sum()) - int(z[L:L + k].sum())
            bad2 += (int(h > 0) != b2)
        else:
            # k*Q - (k+1)*P > 0 on z
            bad3 += (int(k * Qz - (k + 1) * Pz > 0) != b2)
    print('   %d random windows: block-mean form mismatches %d '
          '(even Haar %d, odd form %d)' % (tot, bad1, bad2, bad3))

    print('== V2  shift invariance: bit_b2(a) == bit_b2(a+c) ==')
    bad = 0
    checked = 0
    for _ in range(20000):
        n = int(rng.integers(2, 40))
        a = rng.integers(0, 3, size=n)
        c = int(rng.integers(-1, 3))
        if a.min() + c < 0:
            continue
        L = int(rng.integers(0, n - 1))
        r = int(rng.integers(L + 1, n))
        checked += 1
        bad += (bit_b2(a[L:r + 1]) != bit_b2(a[L:r + 1] + c))
    print('   %d windows, mismatches %d' % (checked, bad))

    print('== V3/V4  epoch decomposition vs walk_b0, one-pass epoch solver ==')
    for n in ns:
        tapes = list(enum_tapes(n))
        agree = onepass_bad = 0
        anchor_strict = True
        halt_cond_ok = True
        for a in tapes:
            z = a - 1
            d0, s0 = walk_b0(a, n)
            d1, E = epoch_walk(z)
            ok = (d0 is None and d1 is None) or (
                d0 is not None and d1 is not None and abs(d0 - d1) < 1e-12)
            agree += ok
            Ls = [e[0] for e in E]
            if any(Ls[i] >= Ls[i + 1] for i in range(len(Ls) - 1)):
                anchor_strict = False
            for (L, R, M, ms) in E:
                got = epoch_onepass(z, L, R)
                if got != ms:
                    onepass_bad += 1
            # halt condition: last epoch has mstar=0 and, at that anchor, all
            # m in [2,M] have b2 = 1
            L, R, M, ms = E[-1]
            if ms == 0:
                ZL = int(z[:L].sum())
                Z = np.concatenate(([0], np.cumsum(z)))
                for m in range(2, M + 1):
                    k = m // 2
                    if not ((Z[L + m] - Z[L]) * k > (Z[L + k] - Z[L]) * m):
                        halt_cond_ok = False
        print('   n=%-3d tapes=%-7d delta agreement %d/%d  anchors strictly '
              'increasing %s  one-pass epoch mismatches %d  halt cond %s'
              % (n, len(tapes), agree, len(tapes), anchor_strict,
                 onepass_bad, halt_cond_ok))


# --------------------------------------------------------------------------- #
def cmd_measure(ns):
    print('== V5  epoch accounting and the model-R quantity ==')
    print('   n     tapes  E/n    sum l / n^2  mean l/ n  max l/ n   '
          'r=0      r=1      r=2      r=4     (sum_e max(0,l_e-2r)/n^2)')
    for n in ns:
        if n <= 13:
            tapes = list(enum_z(n))
        else:
            tapes = list(rand_z(n, 200 if n <= 64 else 60))
        Esum = 0.0
        Lsum = 0.0
        Lmax = 0
        eq = {0: 0.0, 1: 0.0, 2: 0.0, 4: 0.0}
        N = len(tapes)
        for z in tapes:
            d, E = epoch_walk(z)
            Esum += len(E)
            s = 0
            for (L, R, M, ms) in E:
                if ms == 0:
                    le = M - 1        # probes m = M..2
                else:
                    le = M - ms + 1
                s += le
                Lmax = max(Lmax, le)
                for r in eq:
                    eq[r] += max(0, le - 2 * r)
            Lsum += s
        print('   %-5d %-6d %-6.3f %-12.4f %-10.3f %-10.3f %-8.5f %-8.5f %-8.5f %-8.5f'
              % (n, N, Esum / N / n, Lsum / N / n ** 2, Lsum / Esum / n,
                 Lmax / n, eq[0] / N / n ** 2, eq[1] / N / n ** 2,
                 eq[2] / N / n ** 2, eq[4] / N / n ** 2))


# --------------------------------------------------------------------------- #
def _H(counts, tot):
    return -sum((c / tot) * math.log2(c / tot) for c in counts if c)


def cmd_rounds(ns):
    """Round-additivity: H(A_1..A_k) vs sum_{e<=k} H(A_e) for the epoch answers."""
    print('== V7  the round-additivity gap (epoch answers A_e = m*_e) ==')
    print('   n   tapes   k   #tapes>=k   sum_e H(A_e)   H(A_1..A_k)   '
          'ratio   naive reads / log2(3)   true reads <= H(z)/log2(3)')
    for n in ns:
        tapes = list(enum_z(n))
        seq = []
        for z in tapes:
            d, E = epoch_walk(z)
            seq.append(tuple(e[3] for e in E))
        for k in (1, 2, 3, 4, 6):
            sub = [s[:k] for s in seq if len(s) >= k]
            if len(sub) < 50:
                continue
            tot = len(sub)
            marg = []
            for j in range(k):
                c = {}
                for s in sub:
                    c[s[j]] = c.get(s[j], 0) + 1
                marg.append(_H(c.values(), tot))
            cj = {}
            for s in sub:
                cj[s] = cj.get(s, 0) + 1
            Hj = _H(cj.values(), tot)
            print('   %-3d %-7d %-3d %-11d %-14.3f %-13.3f %-7.3f %-23.2f %.2f'
                  % (n, len(tapes), k, tot, sum(marg), Hj, Hj / sum(marg) if sum(marg) else 0,
                     sum(marg) / math.log2(3), n * math.log2(3) / math.log2(3)))
        # the *full* sum of marginals, padding with a 'stop' symbol
        Emax = max(len(s) for s in seq)
        tot = len(seq)
        full = 0.0
        for j in range(Emax):
            c = {}
            for s in seq:
                v = s[j] if j < len(s) else -1
                c[v] = c.get(v, 0) + 1
            full += _H(c.values(), tot)
        print('   %-3d %-7d ALL %-11d %-14.3f %-13s %-7s %-23.2f %.2f'
              % (n, len(tapes), Emax, full, '(joint n/a)', '-',
                 full / math.log2(3), n))


def cmd_gaps(ns):
    """How often a contrast can be decided without reading both blocks fully."""
    print('== V8  contrast gap distribution (can partial reads decide?) ==')
    print('   n    queries   frac m>=M/2   frac |gap|<=sqrt(k)   '
          'mean |gap|/k   frac |gap|=0')
    for n in ns:
        tapes = list(enum_z(n)) if n <= 12 else list(rand_z(n, 100))
        tot = 0
        long_ = 0
        small = 0
        zero = 0
        rel = []
        for z in tapes:
            Z = np.concatenate(([0], np.cumsum(z)))
            d, E = epoch_walk(z)
            for (L, R, M, ms) in E:
                hi = M if ms == 0 else ms
                for m in range(M, max(hi, 1) - 1, -1):
                    if m < 2:
                        continue
                    k = m // 2
                    gap = int((Z[L + m] - Z[L + k]) - (Z[L + k] - Z[L]))
                    tot += 1
                    long_ += (2 * m >= M)
                    small += (abs(gap) <= math.sqrt(k))
                    zero += (gap == 0)
                    rel.append(abs(gap) / k)
        print('   %-4d %-9d %-13.3f %-21.3f %-14.3f %.4f'
              % (n, tot, long_ / tot, small / tot, np.mean(rel), zero / tot))


def cmd_hist(ns):
    """Where does the epoch scan stop?  Distribution of m*/M at the first epoch."""
    print('== V9  first-epoch answer distribution (why H(A_e) is small) ==')
    print('   n    tapes   H(m*/M, 20 bins)   mean m*/M   frac m*=M   frac m*<=M/4'
          '   E/n')
    for n in ns:
        tapes = list(enum_z(n)) if n <= 12 else list(rand_z(n, 300))
        vals = []
        mstarM = []
        fracM = 0
        fracQ = 0
        Es = 0
        for z in tapes:
            d, E = epoch_walk(z)
            Es += len(E)
            L, R, M, ms = E[0]
            vals.append(round(ms / M, 2))
            mstarM.append(ms / M)
            fracM += (ms == M)
            fracQ += (ms <= M / 4)
        c = {}
        for v in vals:
            c[v] = c.get(v, 0) + 1
        tot = len(vals)
        print('   %-4d %-7d %-19.3f %-11.3f %-11.3f %-15.3f %.3f'
              % (n, tot, _H(c.values(), tot), np.mean(mstarM), fracM / tot,
                 fracQ / tot, Es / tot / n))


# --------------------------------------------------------------------------- #
def cmd_lossy(ns):
    print('== V6  lossiness: delta vs the first-anchor contrast sign vector ==')
    print('   n   tapes   #delta   #signatures   mean #signatures per delta   '
          'H(delta)  H(sig)')
    for n in ns:
        tapes = list(enum_z(n))
        dmap = {}
        sigs = set()
        bydelta = {}
        for z in tapes:
            d, E = epoch_walk(z)
            key = None if d is None else round(d * (n - 1))
            Z = np.concatenate(([0], np.cumsum(z)))
            sig = []
            for m in range(2, n + 1):
                k = m // 2
                sig.append(1 if (Z[m] * k > Z[k] * m) else 0)
            sig = tuple(sig)
            sigs.add(sig)
            dmap[key] = dmap.get(key, 0) + 1
            bydelta.setdefault(key, set()).add(sig)
        tot = len(tapes)
        Hd = -sum((c / tot) * math.log2(c / tot) for c in dmap.values())
        Hs = math.log2(len(sigs))
        mean_sig = np.mean([len(v) for v in bydelta.values()])
        print('   %-4d %-7d %-8d %-13d %-28.1f %-9.2f %-9.2f'
              % (n, tot, len(dmap), len(sigs), mean_sig, Hd, Hs))


# --------------------------------------------------------------------------- #
if __name__ == '__main__':
    cmd = sys.argv[1] if len(sys.argv) > 1 else 'verify'
    ns = [int(x) for x in sys.argv[2:]]
    if cmd == 'verify':
        cmd_verify(ns)
    elif cmd == 'measure':
        cmd_measure(ns)
    elif cmd == 'lossy':
        cmd_lossy(ns)
    elif cmd == 'rounds':
        cmd_rounds(ns)
    elif cmd == 'gaps':
        cmd_gaps(ns)
    elif cmd == 'hist':
        cmd_hist(ns)
    else:
        print(__doc__)
