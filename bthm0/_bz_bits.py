# -*- coding: utf-8 -*-
"""How many bits of workspace does the exact incremental simulator actually hold?

Tracks the largest magnitude reached by every state variable on real runs, and
reports ceil(log2) bits for each, for several n.

Usage: python _bz_bits.py <n> <tapes> [seed]
"""
import sys

import numpy as np

from _bz_inc import make_tape, Tap, Win


def run(T, n, acc):
    L, R = 0, n - 1
    w = Win(T, 0, n - 1)
    steps = 0
    while True:
        if w.y >= w.x:
            break
        if R - L + 1 < 3:
            break
        steps += 1
        if steps > 20 * n:
            break
        for name in ('y', 'x', 'S', 'Q', 'Rr', 'Sh', 'hy', 'hx', 'cy', 'cx'):
            acc[name] = max(acc.get(name, 0), abs(getattr(w, name)))
        acc['L'] = max(acc.get('L', 0), abs(L))
        acc['R'] = max(acc.get('R', 0), abs(R))
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
            b1 = False if (np_ * C_ - A * A <= 0 or np_ * D_ - B * B <= 0) else (np_ * w.Rr - A * B) > 0
        o = ("B" if b1 else "D") if b2 else ("C" if b1 else "A")
        dy, dx = {"A": (-1, 0), "B": (1, 0), "C": (0, 1), "D": (0, -1)}[o]
        ny, nx = w.y + dy, w.x + dx
        if ny < L or nx > R:
            if ny < L:
                R = w.x - 1
                w.del_right()
            else:
                L = w.y + 1
                w.del_left()
            w.sync_head()
            continue
        if ny >= nx:
            break
        if nx < L or ny > R:
            L, R = L + 1, R - 1
            w = Win(T, L, R)
            continue
        w.move(dy, dx)


def main():
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 256
    NW = int(sys.argv[2]) if len(sys.argv) > 2 else 60
    seed = int(sys.argv[3]) if len(sys.argv) > 3 else 1
    rng = np.random.default_rng(seed)
    acc = {}
    for _ in range(NW):
        a = make_tape(rng, n)
        run(Tap(np.round(a * 2).astype(np.int32)), n, acc)
    bits = 0
    print('n=%d  tapes=%d' % (n, NW))
    for k in sorted(acc):
        b = int(np.ceil(np.log2(max(1, acc[k]) + 1)))
        bits += b
        print('   %-3s max %6d  -> %2d bits' % (k, acc[k], b))
    print('   total: %d bits of state  =  %.1f * log2(n)   (log2 n = %.1f)'
          % (bits, bits / np.log2(n), np.log2(n)))


if __name__ == '__main__':
    main()
