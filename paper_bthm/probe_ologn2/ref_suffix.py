# -*- coding: utf-8 -*-
"""Referee check 3:
 (i) cost is translation invariant; an unbalanced multiset has the same cost law
     as its centred version (finite-M check).
 (ii) a concrete tape whose suffix is all-1s: cost(i,n-1)=1, below c0*(n-i).
 (iii) conditional mean cost of the suffix given the prefix, n=12.
"""
import numpy as np
from itertools import combinations
from q0_walk import prefix, epoch


def cost_of(word):
    P = prefix(np.asarray(word))
    return epoch(P, 0, len(word) - 1)[1]


def sample_arrangements(counts, ntapes, seed):
    """uniform arrangements of the multiset given by counts (0,1,2)."""
    rng = np.random.default_rng(seed)
    base = np.concatenate([np.zeros(counts[0], int), np.ones(counts[1], int),
                           np.full(counts[2], 2, int)])
    for _ in range(ntapes):
        yield rng.permutation(base)


print("(i) cost law for balanced vs strongly unbalanced multisets, M=256, 4000 arrangements")
for counts, label in [((64, 128, 64), "balanced 1/4,1/2,1/4"),
                      ((26, 128, 102), "unbalanced .10/.50/.40"),
                      ((102, 128, 26), "unbalanced .40/.50/.10")]:
    M = sum(counts)
    cs = np.array([cost_of(w) for w in sample_arrangements(counts, 4000, 7)])
    print("   %-26s E[cost]/M=%.3f  P(cost>=.05M)=%.3f" % (label, cs.mean() / M, (cs >= .05 * M).mean()))

# translation invariance, brute force
w = np.array([0, 2, 0, 1, 2, 1, 1, 2, 0, 1])
print("\n(ii) translation invariance, brute-force check on one word:",
      cost_of(w), cost_of(w + 1), cost_of(w + 2))

# (ii) the all-1s suffix tape
for n in (64, 512):
    k = 3
    a = np.concatenate([np.zeros(k, int), np.full(k, 2, int), np.ones(n - 2 * k, int)])
    print("   n=%4d  a=[0]*%d+[2]*%d+[1]*%d : #2-#0=%d  cost(6,n-1)=%d  c0*(n-6)=%.2f (c0=.008) / %.2f (c0=.004)"
          % (n, k, k, n - 2 * k, (a == 2).sum() - (a == 0).sum(), cost_of(a[2 * k:]),
             .008 * (n - 2 * k), .004 * (n - 2 * k)))

# (iii) conditional mean given prefix, n=12 exhaustive
print("\n(iii) n=12 exhaustive: E[cost(i,n-1) | prefix] over prefixes, i=6")
from collections import defaultdict


def tapes_even(n):
    for k in range(n // 2 + 1):
        for pos2 in combinations(range(n), k):
            rest = [j for j in range(n) if j not in pos2]
            for pos0 in combinations(rest, k):
                a = np.ones(n, dtype=np.int64)
                a[list(pos2)] = 2
                a[list(pos0)] = 0
                yield a


i = 6
groups = defaultdict(list)
for a in tapes_even(12):
    groups[tuple(a[:i])].append(cost_of(a[i:]))
means = np.array([np.mean(v) for v in groups.values()])
print("   prefixes=%d  E over prefixes: %.3f   min=%.3f  max=%.3f   fraction with mean<2: %.2f"
      % (len(means), means.mean(), means.min(), means.max(), (means < 2).mean()))
allc = np.array([c for v in groups.values() for c in v])
print("   unconditional E[cost(i,n-1)]/(n-i) = %.3f" % (allc.mean() / (12 - i)))
