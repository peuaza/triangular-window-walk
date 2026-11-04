# -*- coding: utf-8 -*-
"""Epoch structure of the four rules: how many cuts, how long is each epoch's
scan?  Supports (or refutes) the picture 'Theta(n) queries, each of read
complexity Theta(n)' that is the honest reason for the n^2."""
import sys
import numpy as np

PAT = np.array([2, 1, 1])


def longest_run(z):
    b = (z > 1).astype(np.int8)
    if not b.any():
        return 0
    d = np.diff(np.concatenate((np.zeros(1, np.int8), b, np.zeros(1, np.int8))))
    return int((np.flatnonzero(d == -1) - np.flatnonzero(d == 1)).max())


def bit_b1(z, rule):
    m = len(z)
    if rule == 'longrun':
        return longest_run(z) * 4 > m
    if rule == 'runs':
        return int((np.diff(z) != 0).sum()) * 2 > m
    if rule == 'pattern':
        if m < 3:
            return False
        return bool((np.lib.stride_tricks.sliding_window_view(z, 3) == PAT).all(1).any())
    if rule == 'pal':
        return bool((z == z[::-1]).all())
    raise ValueError(rule)


def walk(a, n, rule):
    L, R = 0, n - 1
    ell, r = 0, n - 1
    seen = set()
    steps = 0
    epochs = []
    cur = 0
    while True:
        if ell >= r:
            epochs.append((cur, R - L + 1))
            return ell / (n - 1), steps, epochs
        if R - L + 1 < 3:
            epochs.append((cur, R - L + 1))
            return None, steps, epochs
        steps += 1
        if steps > 20 * n * n:
            epochs.append((cur, R - L + 1))
            return None, steps, epochs
        if (ell, r) in seen:
            epochs.append((cur, R - L + 1))
            return None, steps, epochs
        seen.add((ell, r))
        cur += 1
        z = a[ell:r + 1]
        m = r - ell + 1
        k = max(1, m // 2)
        b2 = int(z.sum()) * k > int(z[:k].sum()) * m
        b1 = bit_b1(z, rule)
        if b2 and b1:
            dy, dx = 1, 0
        elif b2 and not b1:
            dy, dx = 0, -1
        elif b1 and not b2:
            dy, dx = 0, 1
        else:
            dy, dx = -1, 0
        if ell + dy < L or r + dx > R:
            epochs.append((cur, R - L + 1))
            cur = 0
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
            epochs.append((cur, R - L + 1))
            return (ell + dy) / (n - 1), steps, epochs
        if r + dx < L or ell + dy > R:
            epochs.append((cur, R - L + 1))
            cur = 0
            L, R = L + 1, R - 1
            ell, r = L, R
            seen.clear()
            continue
        ell, r = ell + dy, r + dx


if __name__ == '__main__':
    for N in [int(x) for x in sys.argv[1:]] or [128, 256]:
        NW, H = max(30, 3000 // N), N // 2
        print('n=%d tapes=%d' % (N, NW))
        for rule in ['longrun', 'runs', 'pattern', 'pal']:
            rng = np.random.default_rng(2026)
            ep, sl, tot = [], [], []
            for _ in range(NW):
                B = np.zeros((N, 2), dtype=np.int16)
                for c in range(2):
                    B[rng.choice(N, H, replace=False), c] = 1
                a = B.sum(1).astype(np.int64)
                d, s, epochs = walk(a, N, rule)
                ep.append(len(epochs))
                if epochs:
                    sl.append(np.mean([e[0] for e in epochs]))
                tot.append(s)
            print('  %-8s epochs=%6.1f (%.2f n)  mean scan=%6.1f (%.2f n)  steps=%7.1f'
                  % (rule, np.mean(ep), np.mean(ep) / N, np.mean(sl), np.mean(sl) / N,
                     np.mean(tot)))
            sys.stdout.flush()
