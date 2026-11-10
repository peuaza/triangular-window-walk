"""Fit the witness panel's g and compare with e and the sample maximum.

Outputs, per rule:
  g_wit(n)   max over the structured witnesses that meet
  g_K(n)     the sample maximum from probe_master/master_table.json
  e(n)       the mean from the same table
  exponents  log-log OLS over the sizes present
  g_wit/e    at the largest size
"""
import json
import math
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
PAPER = os.path.dirname(HERE)


def exp_of(xs, ys):
    xs = np.asarray(xs, float)
    ys = np.asarray(ys, float)
    if xs.size < 2 or np.any(ys <= 0):
        return float('nan')
    return float(np.polyfit(np.log(xs), np.log(ys), 1)[0])


def main():
    cells = json.load(open(os.path.join(HERE, 'gw_cells.json')))
    rows = json.load(open(os.path.join(PAPER, 'probe_master', 'master_table.json')))['rows']
    tab = {}
    for r in rows:
        tab[(r['rule'], r['n'])] = r
    ns = sorted(int(k) for k in cells)
    rules = sorted({k.split('|')[0] for k in cells[str(ns[0])]})
    out = []
    print(f"{'rule':<10}{'g_wit exp':>10}{'g_K exp':>9}{'e exp':>8}{'g_wit/e':>9}{'n':>6}")
    for rule in rules:
        gw, gk, ev, nn = [], [], [], []
        for n in ns:
            d = cells[str(n)]
            vals = [v['reads'] for k, v in d.items()
                    if k.split('|')[0] == rule and 'err' not in v and v.get('delta') is not None]
            if not vals:
                continue
            gw.append(max(vals))
            r = tab.get((rule, n))
            if r:
                gk.append(r.get('g_meet') or float('nan'))
                ev.append(r.get('e') or float('nan'))
                nn.append(n)
        if len(nn) < 2:
            continue
        rec = dict(rule=rule, ns=nn, g_wit=gw[-len(nn):], g_K=gk, e=ev,
                   g_wit_exp=exp_of(nn, gw[-len(nn):]), g_K_exp=exp_of(nn, gk),
                   e_exp=exp_of(nn, ev), ratio=gw[-1] / ev[-1] if ev[-1] else float('nan'))
        out.append(rec)
        print(f"{rule:<10}{rec['g_wit_exp']:>10.2f}{rec['g_K_exp']:>9.2f}{rec['e_exp']:>8.2f}"
              f"{rec['ratio']:>9.1f}{nn[-1]:>6}")
    with open(os.path.join(HERE, 'gw_fit.json'), 'w') as f:
        json.dump(out, f, indent=1)


if __name__ == '__main__':
    main()
