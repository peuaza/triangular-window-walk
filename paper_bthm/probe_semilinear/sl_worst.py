# -*- coding: utf-8 -*-
"""sl_worst.py -- targeted worst-case tapes for the threshold-run rule.

Constructs valid ensemble tapes (two columns, h = n/2 ones each, i.e. #2 = #0,
#1 = n - 2#2) carrying a forced 2-run of length t+1 at position p, then measures
the in-class machine.  This gives a deterministic lower bound on g(n) instead of
waiting for a 1/n event in a random sample.
"""
import sys
import numpy as np

from sl_machine import solve, t_of


def make_tape(n, t, p, seed, n2=None):
    rng = np.random.default_rng(seed)
    if n2 is None:
        n2 = n // 4
    a = [1] * n
    L = t + 1
    assert p + L <= n
    for i in range(p, p + L):
        a[i] = 2
    rem2 = n2 - L
    rem0 = n2
    free = [i for i in range(n) if a[i] == 1]
    rng.shuffle(free)
    assert rem2 >= 0 and rem2 + rem0 <= len(free)
    for i in free[:rem2]:
        a[i] = 2
    for i in free[rem2:rem2 + rem0]:
        a[i] = 0
    return a


def main():
    ns = [int(x) for x in sys.argv[1].split(',')]
    step = int(sys.argv[2]) if len(sys.argv) > 2 else 8
    seeds = [int(x) for x in sys.argv[3].split(',')] if len(sys.argv) > 3 else [1]
    frac = float(sys.argv[4]) if len(sys.argv) > 4 else 0.25
    print('%5s %4s %4s %8s %12s %10s %8s %10s' %
          ('n', 't', 'p', 'seed', 'steps', 'steps/n', 'reads', 'reads/n2'))
    for n in ns:
        t = t_of(n)
        best = None
        for seed in seeds:
            for p in range(0, n - t, step):
                a = make_tape(n, t, p, seed, n2=int(n * frac))
                r = solve(a, t, store=True)
                if best is None or r['steps'] > best['steps']:
                    best = r
                    best['p'] = p
                    best['seed'] = seed
        r = best
        print('%5d %4d %4d %8d %12d %10.1f %8d %10.3f' %
              (n, t, r['p'], r['seed'], r['steps'], r['steps'] / n,
               r['reads'], r['reads'] / n ** 2), flush=True)


if __name__ == '__main__':
    main()
