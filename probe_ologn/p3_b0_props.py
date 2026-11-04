# -*- coding: utf-8 -*-
"""p3: (a) test simplifications of the b1=0 process, (b) where does b1 fire, (c) full vs b0."""
import sys, itertools
import numpy as np
sys.path.insert(0, '..')
from _arr_longrun_fast import walk_naive, bit_b2_ref, bit_b1_ref, longest_run_vec
from p2_b0 import walk_b0, scan_for_L


def enum_tapes(n):
    for t in range(n // 2 + 1):
        c1 = n - 2 * t
        if c1 < 0:
            continue
        for zeros in itertools.combinations(range(n), t):
            zs = set(zeros)
            rest = [i for i in range(n) if i not in zs]
            for twos in itertools.combinations(rest, t):
                a = np.ones(n, dtype=np.int64)
                a[list(zeros)] = 0
                a[list(twos)] = 2
                yield a


def delta_b0_simplified(a, n):
    """min L such that b2([L,r])=1 for all r in [L+1, n-1]  (R fixed at n-1)."""
    for L in range(n - 1):
        ok = True
        for r in range(L + 1, n):
            if bit_b2_ref(a[L:r + 1]) == 0:
                ok = False
                break
        if ok:
            return L / (n - 1)
    return None


def feature_test(n):
    tapes = list(enum_tapes(n))
    simple_hit = 0
    tot = 0
    for a in tapes:
        d0, _ = walk_b0(a, n)
        ds = delta_b0_simplified(a, n)
        tot += 1
        if (d0 is None and ds is None) or (d0 is not None and ds is not None and abs(d0 - ds) < 1e-12):
            simple_hit += 1
    print('n=%d: "R frozen at n-1" closed form matches b1=0 walk on %d/%d' % (n, simple_hit, tot))


def b1_stats(n, nw, seed):
    """how often / at which window sizes does b1 fire on random ensemble tapes?"""
    rng = np.random.default_rng(seed)
    h = n // 2 if n % 2 == 0 else (n + 1) // 2
    tot_steps = tot_true = 0
    ms = []
    per_tape = []
    maxm = 0
    for _ in range(nw):
        B = np.zeros((n, 2), dtype=np.int16)
        for c in range(2):
            B[rng.choice(n, h, replace=False), c] = 1
        a = B.sum(1).astype(np.int64)
        L, R = 0, n - 1
        ell, r = 0, n - 1
        seen = set()
        steps = tru = 0
        while True:
            if ell >= r:
                break
            if R - L + 1 < 3:
                break
            steps += 1
            if steps > 20 * n * n:
                break
            if (ell, r) in seen:
                break
            seen.add((ell, r))
            m = r - ell + 1
            b2 = bit_b2_ref(a[ell:r + 1])
            b1 = bit_b1_ref(a[ell:r + 1])
            if b1:
                tru += 1
                ms.append(m)
                maxm = max(maxm, m)
            if b2 and b1:
                dy, dx = 1, 0
            elif b2 and not b1:
                dy, dx = 0, -1
            elif b1 and not b2:
                dy, dx = 0, 1
            else:
                dy, dx = -1, 0
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
                break
            if r + dx < L or ell + dy > R:
                L, R = L + 1, R - 1
                ell, r = L, R
                seen.clear()
                continue
            ell, r = ell + dy, r + dx
        tot_steps += steps
        tot_true += tru
        per_tape.append(tru)
    ms = np.array(ms) if ms else np.array([0])
    print('n=%d tapes=%d: steps=%d b1true=%d (%.5f of steps)  per-tape mean %.2f max %d  m at b1true: median %.2f q90 %.2f max %d'
          % (n, nw, tot_steps, tot_true, tot_true / max(1, tot_steps), np.mean(per_tape), max(per_tape),
             np.median(ms), np.quantile(ms, .9), maxm))


def full_vs_b0(n, nw, seed):
    rng = np.random.default_rng(seed)
    h = n // 2 if n % 2 == 0 else (n + 1) // 2
    same = 0
    for _ in range(nw):
        B = np.zeros((n, 2), dtype=np.int16)
        for c in range(2):
            B[rng.choice(n, h, replace=False), c] = 1
        a = B.sum(1).astype(np.int64)
        d1, s1, _ = walk_naive(a, n)
        d2, s2 = walk_b0(a, n)
        if (d1 is None and d2 is None) or (d1 is not None and d2 is not None and abs(d1 - d2) < 1e-12):
            same += 1
    print('n=%d tapes=%d: full delta == b1=0 delta on %d/%d' % (n, nw, same, nw))


if __name__ == '__main__':
    for n in (6, 8, 10, 11, 12):
        feature_test(n)
    for n in (64, 128, 256, 512):
        b1_stats(n, 40, 2026)
        full_vs_b0(n, 40, 2026)
