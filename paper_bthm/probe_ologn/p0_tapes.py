# -*- coding: utf-8 -*-
"""p0: correct ensemble enumeration.

a[i] = x_i + y_i, x,y binary columns each with ceil(n/2) ones.
n even: #0 = #2 = t, #1 = n-2t
n odd : #0 = #2+1  (c0=t+1, c2=t, c1=n-2t-1)
"""
import itertools
import numpy as np


def enum_tapes(n):
    if n % 2 == 0:
        for t in range(n // 2 + 1):
            c0 = c2 = t
            c1 = n - 2 * t
            if c1 < 0:
                continue
            for zeros in itertools.combinations(range(n), c0):
                zs = set(zeros)
                rest = [i for i in range(n) if i not in zs]
                for twos in itertools.combinations(rest, c2):
                    a = np.ones(n, dtype=np.int64)
                    a[list(zeros)] = 0
                    a[list(twos)] = 2
                    yield a
    else:
        for t in range(n // 2 + 1):
            c0, c2 = t + 1, t
            c1 = n - 2 * t - 1
            if c1 < 0:
                continue
            for zeros in itertools.combinations(range(n), c0):
                zs = set(zeros)
                rest = [i for i in range(n) if i not in zs]
                for twos in itertools.combinations(rest, c2):
                    a = np.ones(n, dtype=np.int64)
                    a[list(zeros)] = 0
                    a[list(twos)] = 2
                    yield a


def count_tapes(n):
    import math
    tot = 0
    if n % 2 == 0:
        for t in range(n // 2 + 1):
            c1 = n - 2 * t
            if c1 >= 0:
                tot += math.factorial(n) // (math.factorial(t) * math.factorial(c1) * math.factorial(t))
    else:
        for t in range(n // 2 + 1):
            c0, c2, c1 = t + 1, t, n - 2 * t - 1
            if c1 >= 0:
                tot += math.factorial(n) // (math.factorial(c0) * math.factorial(c1) * math.factorial(c2))
    return tot


if __name__ == '__main__':
    import sys
    for n in range(3, 16):
        print(n, count_tapes(n))
    for n in (6, 7, 8):
        print('enum', n, sum(1 for _ in enum_tapes(n)))
