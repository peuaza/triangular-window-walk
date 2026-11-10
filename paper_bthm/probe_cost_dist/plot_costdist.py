"""Figure: per-tape cost distributions, one ridge per rule, fully-linear above the band family.

x = log10(reads / n) on the meeting tapes at n = 4096.  The four fully-linear rules concentrate in a
narrow band; the band rules spread over two to four decades and carry their mass in the right tail.
Writes probe_cost_dist/costdist.pdf and .png.
"""
import glob
import json
import os

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde

HERE = os.path.dirname(os.path.abspath(__file__))
PAPER = os.path.dirname(HERE)
N = 4096

FL = ['monotone', 'pattern', 'runs', 'peak']
BAND = ['KCHG', 'BAL>m/2', 'BAL>m/4', 'NBLK>2', 'NBLK>3', 'NBLK=1',
        'AMIN', 'AMAX', 'LAS', 'MAXP', 'oddrun']

data, meta = {}, {}
p = os.path.join(HERE, 'fl_costs.json')
if os.path.exists(p):
    d = json.load(open(p))
    for k, v in d.items():
        if v:
            data[k] = np.array(v, float)
            meta[k] = 'fully linear'
for f in glob.glob(os.path.join(PAPER, 'probe_k300', 'k3_*_%d.json' % N)):
    j = json.load(open(f))
    r = j['rule']
    if r in BAND:
        data[r] = np.array(j['reads'], float)
        meta[r] = 'band'

rows = [r for r in FL + BAND if r in data]
fig, ax = plt.subplots(figsize=(7.2, 0.30 * len(rows) + 0.9))
xs_all = []
for k, r in enumerate(rows):
    c = data[r] / N
    z = np.log10(c[c > 0])
    if z.size < 5:
        continue
    lo, hi = np.percentile(z, [0.5, 99.5])
    grid = np.linspace(min(lo - 0.3, z.min()), max(hi + 0.3, z.max()), 400)
    if z.std() < 1e-6:
        dens = np.exp(-((grid - z[0]) / 0.02) ** 2)
    else:
        dens = gaussian_kde(z, bw_method=0.35)(grid)
    dens = dens / dens.max() * 0.85
    col = '#1f4e79' if meta[r] == 'fully linear' else '#b45309'
    ax.fill_between(grid, k, k + dens, color=col, alpha=0.55, lw=0)
    ax.plot([np.median(z)] * 2, [k, k + 0.55], color='k', lw=1.0, alpha=0.8)
    ax.text(grid.min() - 0.06, k + 0.28, r, ha='right', va='center', fontsize=8,
            family='DejaVu Sans Mono')
    xs_all.append(grid)

lo = min(g.min() for g in xs_all) if xs_all else -1
hi = max(g.max() for g in xs_all) if xs_all else 3
ax.set_xlim(lo - 0.05, hi + 0.1)
ax.set_ylim(-0.4, len(rows) + 0.2)
ax.set_yticks([])
ax.set_xlabel(r'$\log_{10}(\mathrm{reads}/n)$ on the meeting tapes, $n=4096$')
for s in ('top', 'right', 'left'):
    ax.spines[s].set_visible(False)
ax.grid(axis='x', ls=':', lw=0.5, alpha=0.6)
ax.set_title('Cost per cell, one density per rule (vertical tick: median)', fontsize=9, loc='left')
fig.tight_layout()
fig.savefig(os.path.join(HERE, 'costdist.pdf'))
fig.savefig(os.path.join(HERE, 'costdist.png'), dpi=170)
print('rows', len(rows), 'wrote', os.path.join(HERE, 'costdist.pdf'))
for r in rows:
    c = data[r] / N
    print('%-9s n=%3d  min %.3g  med %.3g  p99 %.3g  max %.3g' %
          (r, c.size, c.min(), np.median(c), np.percentile(c, 99), c.max()))
