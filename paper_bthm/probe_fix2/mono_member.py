"""mono_member -- a content-dependent member of the monotone family, and its gate.

The measured member (`lin_engine.mono_walk`) moves the headers by (ell+r) parity alone: its output is a
function of n, so a zero-read machine answers it exactly and its rule-level cost is 0. This probe
implements the intended family member -- f reads the window's sum -- and checks:

  gate (tape dependence)  the output must vary with the tape; a member whose output is constant is
                          degenerate and is reported as content-independent
  cost                    the collapse machine reads the opening window (n) plus one cell per step
                          (n-1 steps), so 2n-1 on every tape: G = e = g = 2n-1 at every n
  lower bound (sketch)    every step needs the cell that leaves the window, so any exact machine reads
                          at least n-1 cells on every tape: G_rho >= n-1, hence Theta(n)

Usage: python mono_member.py [sizes] [K]
"""
import sys
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
PAPER = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(PAPER, 'probe_en_linear'))

import el_lib as E     # noqa: E402


def mono_member(a, n):
    """f: move the left header iff the window's sum exceeds its length (mean > 1); else the right."""
    reads = n
    ell, r = 0, n - 1
    s = int(np.sum(a))
    steps = 0
    while ell < r:
        m = r - ell + 1
        if s > m:
            s -= int(a[ell]); ell += 1
        else:
            s -= int(a[r]); r -= 1
        reads += 1
        steps += 1
    return dict(out=ell / (n - 1), reads=reads, steps=steps)


def main():
    sizes = [int(v) for v in (sys.argv[1].split(',') if len(sys.argv) > 1 else ['512', '1024', '2048', '4096'])]
    K = int(sys.argv[2]) if len(sys.argv) > 2 else 200
    print("%6s %6s %8s %8s %8s %8s %8s" % ("n", "K", "distinct", "sd(out)", "min reads", "max reads", "e/n"))
    for n in sizes:
        outs, reads = [], []
        for i in range(K):
            X = E.tape(n, i)
            a = (np.asarray(X)[:, 0] + np.asarray(X)[:, 1]).astype(np.int64)
            r = mono_member(a, n)
            outs.append(r['out']); reads.append(r['reads'])
        outs = np.asarray(outs, float); reads = np.asarray(reads, float)
        print("%6d %6d %8d %8.4f %8.0f %8.0f %8.3f" %
              (n, K, len(set(np.round(outs, 12))), outs.std(), reads.min(), reads.max(), reads.mean() / n))
        assert reads.min() == 2 * n - 1 and reads.max() == 2 * n - 1


if __name__ == '__main__':
    main()
