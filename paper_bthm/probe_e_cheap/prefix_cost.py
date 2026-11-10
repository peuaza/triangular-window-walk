"""Where does a hard tape's cost sit?  First epoch, or spread out?

The user's question: can a slightly larger budget suppress the residual from the hard tapes?  A
predictor can only stop early if stopping early saves cost, and it can only save cost if the early
epochs are cheap.  Each epoch costs about (its length) x (its window width), so we read the epoch
trace and report, for the hard tapes (total cost > 10x the median), the share of the cost carried by
the first epoch, the first two, and the last one.

If the first epoch already dominates, no early-stop shortcut exists and the hard tapes must be paid
in full.

Usage: python prefix_cost.py [nt]
"""
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
PAPER = os.path.dirname(HERE)
for d in (os.path.join(PAPER, 'probe_en_linear'), os.path.join(PAPER, 'probe_triple_nat'),
          os.path.join(PAPER, 'probe_e_uniform')):
    if d not in sys.path:
        sys.path.insert(0, d)

import el_lib as E                                              # noqa: E402

RULES = ('AMAX', 'AMIN', 'BAL>m/2', 'BAL>m/4', 'KCHG', 'LAS', 'MAXP', 'NBLK=1', 'NBLK>2',
         'NBLK>3', 'dense2', 'oddrun', 'trip')
NT = int(sys.argv[1]) if len(sys.argv) > 1 else 60
N = 1024


def main():
    print(f"{'rule':<10}{'tapes':>6}{'hard':>5}{'cost1/tot':>10}{'cost2/tot':>10}"
          f"{'last/tot':>9}{'nep_hard':>9}")
    for rule in RULES:
        shares1, shares2, sharesL, neps = [], [], [], []
        costs = []
        traces = []
        for i in range(NT * 20):
            if len(traces) >= NT:
                break
            try:
                X = E.tape(N, 5000 + N + 100003 * i)
                M = E.Mach(X, N)
                t = M.trace(rule)
            except Exception:                                   # noqa: BLE001
                continue
            if t['end'] != 0:                                   # meeting tapes only
                continue
            ep = np.asarray(t['eplen'], float) * np.asarray(t['epM'], float)
            if ep.size == 0 or ep.sum() <= 0:
                continue
            costs.append(float(ep.sum()))
            traces.append(ep)
        if not traces:
            print(f'{rule:<10}{0:>6}')
            continue
        med = float(np.median(costs))
        for ep in traces:
            if ep.sum() <= 10 * med:
                continue
            tot = ep.sum()
            shares1.append(ep[0] / tot)
            shares2.append(ep[:2].sum() / tot)
            sharesL.append(ep[-1] / tot)
            neps.append(ep.size)
        if not shares1:
            print(f'{rule:<10}{len(traces):>6}{0:>5}')
            continue
        print(f'{rule:<10}{len(traces):>6}{len(shares1):>5}{np.mean(shares1):>10.2f}'
              f'{np.mean(shares2):>10.2f}{np.mean(sharesL):>9.2f}{np.mean(neps):>9.0f}')


if __name__ == '__main__':
    main()
