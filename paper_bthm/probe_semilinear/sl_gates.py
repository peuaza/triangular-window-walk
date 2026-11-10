# -*- coding: utf-8 -*-
"""sl_gates.py -- the gate measurements.

G1  output distribution (from the run files, printed separately).
G2  size of the exact witness structure: number of maximal 2-runs longer than t,
    on random tapes and on a constructed tape that packs as many as the ensemble
    allows; does the fixed-size register list overflow?
G3  availability of the tape-wide certificate rho <= t.
G4  (support) the naive scan machine on a constructed worst-case tape.
"""
import sys
import numpy as np

from sl_machine import solve, t_of


def packed_valid(n, t):
    """Pack runs of length t+1 with #2 = #0 = n/2 exactly (n even)."""
    L = t + 1
    a = [1] * n
    n2 = n // 2
    cnt = n2 // L
    p = 0
    for _ in range(cnt):
        if p + L > n:
            break
        for i in range(p, p + L):
            a[i] = 2
        p += L + 1
    n2_now = sum(1 for v in a if v == 2)
    # remaining 2s anywhere
    free = [i for i in range(n) if a[i] == 1]
    for i in free[:n2 - n2_now]:
        a[i] = 2
    # zeros: need n2
    free = [i for i in range(n) if a[i] == 1]
    for i in free[:n2]:
        a[i] = 0
    assert sum(1 for v in a if v == 2) == n2
    assert sum(1 for v in a if v == 0) == n2
    return a


def nlong_of(a, t):
    best = 0
    cur = 0
    for v in a:
        if v == 2:
            cur += 1
        else:
            if cur > t:
                best += 1
            cur = 0
    if cur > t:
        best += 1
    return best


def main():
    ns = [int(x) for x in sys.argv[1].split(',')]
    print('%5s %4s %8s %8s %10s %8s %10s %12s' %
          ('n', 't', 'max_runs', 'packed', 'overflow', 'steps', 'reads', 'reads/n2'))
    for n in ns:
        t = t_of(n)
        a = packed_valid(n, t)
        nl = nlong_of(a, t)
        r = solve(a, t, store=True)
        print('%5d %4d %8d %8d %10s %8d %10d %12.3f' %
              (n, t, n // (2 * (t + 1)), nl,
               str(r['overflow']), r['steps'], r['reads'],
               r['reads'] / n ** 2), flush=True)
    # G4: naive scan machine on the constructed worst-case tape
    print('\nG4 naive scan machine on a constructed worst-case tape')
    print('%5s %8s %12s %12s %10s' % ('n', 'steps', 'reads', 'reads/n2', 'reads/n3'))
    for n in [int(x) for x in sys.argv[2].split(',')] if len(sys.argv) > 2 else []:
        t = t_of(n)
        # best p from sl_worst for this n (approx 0.68 n)
        from sl_worst import make_tape
        a = make_tape(n, t, int(0.68 * n), 1)
        r = solve(a, t, store=False)
        print('%5d %8d %12d %12.3f %10.4f' %
              (n, r['steps'], r['reads'], r['reads'] / n ** 2,
               r['reads'] / n ** 3), flush=True)


if __name__ == '__main__':
    main()
