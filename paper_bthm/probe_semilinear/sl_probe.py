# -*- coding: utf-8 -*-
"""sl_probe.py -- quick shape probe for the threshold-run rule.

b1 = [ longest 2-run in the window > t(n) ],  t(n) = floor(log2 n) + c.
If rho (global longest 2-run) <= t, b1 == 0 on every window and the b0 shortcut
applies.  This probe charges:
    init  n + n/2 reads
    per step 1.5 reads
    per test with rho > t and m > t: m reads (full scan, upper bound)
No literal-machine bookkeeping yet; this is only to see the shape of e and g.
"""
import sys
import numpy as np
from collections import Counter

sys.path.insert(0, '../probe_ologn2')
sys.path.insert(0, '.')
from rs_common import Capsule, sample_tapes, b0_shortcut  # noqa: E402


def rho_of(a):
    best = cur = 0
    for v in a:
        if v == 2:
            cur += 1
            if cur > best:
                best = cur
        else:
            cur = 0
    return best


def walk_capsule(a, t, charge=True):
    """Capsule walk: exact b1 = [longest 2-run > t], with read accounting."""
    n = len(a)
    C = Capsule(a)
    rho = rho_of(a)
    reads = n + max(1, n // 2)
    L, R, ell, r = 0, n - 1, 0, n - 1
    seen = set()
    steps = 0
    cap = 40 * n * n
    while True:
        if ell >= r:
            return ell / (n - 1), reads, steps, rho
        if R - L + 1 < 3:
            return None, reads, steps, rho
        steps += 1
        if steps > cap:
            return None, reads, steps, rho
        if (ell, r) in seen:
            return None, reads, steps, rho
        seen.add((ell, r))
        m = r - ell + 1
        b2 = C.b2(ell, m)
        if rho <= t or m <= t:
            b1 = False
            reads += 0
        else:
            b1 = C.longest_run(ell, r) > t
            reads += m if charge else 1
        reads += 1.5
        if b2:
            dy, dx = (1, 0) if b1 else (0, -1)
        else:
            dy, dx = (0, 1) if b1 else (-1, 0)
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
            return (ell + dy) / (n - 1), reads, steps, rho
        if r + dx < L or ell + dy > R:
            L, R = L + 1, R - 1
            ell, r = L, R
            seen.clear()
            continue
        ell, r = ell + dy, r + dx


def main():
    ns = [int(x) for x in sys.argv[1].split(',')]
    N = int(sys.argv[2]) if len(sys.argv) > 2 else 2000
    seed = int(sys.argv[3]) if len(sys.argv) > 3 else 2026
    c = int(sys.argv[4]) if len(sys.argv) > 4 else 0
    print('%5s %4s %6s %10s %8s %10s %8s %10s %8s' %
          ('n', 't', 'rare', 'mean_rd', 'e/n', 'max_rd', 'g/n', 'exc_mean', 'exc/n2'))
    for n in ns:
        t = int(np.floor(np.log2(n))) + c
        rs = []
        exc = []
        rhos = []
        for a in sample_tapes(n, N, seed):
            d, rd, st, rho = walk_capsule(a, t)
            rs.append(rd)
            rhos.append(rho)
            if rho > t:
                exc.append(rd)
        rs = np.array(rs, float)
        print('%5d %4d %6d %10.1f %8.2f %10.1f %8.1f %10.1f %8.3f' %
              (n, t, len(exc), rs.mean(), rs.mean() / n, rs.max(), rs.max() / n,
               np.mean(exc) if exc else float('nan'),
               (np.mean(exc) / n ** 2) if exc else float('nan')), flush=True)


if __name__ == '__main__':
    main()
