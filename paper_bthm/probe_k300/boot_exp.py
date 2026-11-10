"""Bootstrap confidence intervals for the per-rule exponents of e and G (review item 7).

For each rule and each fit, the tapes are resampled with replacement at every size (independent
resampling per size, 400 draws), the statistic is recomputed and the log-log exponent refitted; the
2.5th and 97.5th percentiles of the 400 exponents are the interval.  This is the honest interval for
a fit over four sizes: it answers "how much of the exponent is the tape sample", not "how much is the
functional form".

Usage: python boot_exp.py [draws]
"""
import glob
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
PAPER = os.path.dirname(HERE)
NS = (512, 1024, 2048, 4096)
DRAWS = int(sys.argv[1]) if len(sys.argv) > 1 else 400


def expo(ns, ys):
    ns = np.asarray(ns, float)
    ys = np.asarray(ys, float)
    return float(np.polyfit(np.log(ns), np.log(ys), 1)[0])


def main():
    cells = {}
    for f in glob.glob(os.path.join(PAPER, 'probe_k300', 'k3_*.json')):
        j = json.load(open(f))
        cells.setdefault(j['rule'], {})[j['n']] = np.array(j['reads'], float)
    rng = np.random.default_rng(20261110)
    print(f"{'rule':<10}{'e exp':>7}{'e 95% CI':>16}{'G exp':>7}{'G 95% CI':>16}{'tapes@4k':>9}")
    for rule in sorted(cells):
        ns = [n for n in NS if n in cells[rule]]
        if len(ns) < 3:
            continue
        E = [float(np.mean(cells[rule][n])) for n in ns]
        G = [float(np.exp(np.mean(np.log(cells[rule][n])))) for n in ns]
        be, bg = expo(ns, E), expo(ns, G)
        de, dg = [], []
        for _ in range(DRAWS):
            Es, Gs = [], []
            for n in ns:
                a = cells[rule][n]
                s = a[rng.integers(0, a.size, a.size)]
                Es.append(float(np.mean(s)))
                Gs.append(float(np.exp(np.mean(np.log(s)))))
            de.append(expo(ns, Es))
            dg.append(expo(ns, Gs))
        qe = np.percentile(de, [2.5, 97.5])
        qg = np.percentile(dg, [2.5, 97.5])
        print(f"{rule:<10}{be:>7.2f}[{qe[0]:>6.2f},{qe[1]:>6.2f}]"
              f"{bg:>7.2f}[{qg[0]:>6.2f},{qg[1]:>6.2f}]{cells[rule][ns[-1]].size:>9}")


if __name__ == '__main__':
    main()
