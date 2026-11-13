"""mp_plot -- one figure per rule, and a 21-panel atlas.

For a rule, every machine in mp_cells.json contributes a density of log10(C/n) over the meeting
tapes at the largest n available, and three markers: ln G_M (the geometric mean), e_M (the mean),
g_M (the maximum).  The rule's own G_rho, e_rho, g_rho are the three minima over machines, drawn as
three vertical lines at the top of the panel with the machine that attains each.

Usage:
  python mp_plot.py rule RULE [N]
  python mp_plot.py atlas [N]
Writes probe_minima/figs/<rule>.pdf and probe_minima/mp_atlas.pdf.
"""
import json
import os
import sys

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde

HERE = os.path.dirname(os.path.abspath(__file__))
CELLS = os.path.join(HERE, 'mp_cells.json')
FIGS = os.path.join(HERE, 'figs')
FIG_DIR = FIGS


def load():
    return json.load(open(CELLS)) if os.path.exists(CELLS) else []


def stats(recs, n):
    """{(machine): (G, e, g, reads)} at size n."""
    out = {}
    for r in recs:
        if r['n'] != n or not r['reads']:
            continue
        c = np.asarray(r['reads'], float)
        out[r['machine']] = (float(np.exp(np.mean(np.log(c)))), float(np.mean(c)), float(c.max()), c)
    return out


def panel(ax, rule, recs, n, title=True):
    """Only the machines that attain min G, min e or min g are drawn.

    The scan and pcr engines cost two to three decades more than the predictive ones; drawing them
    compresses every curve that matters into a sliver, so they are excluded by construction and named
    in the caption instead.
    """
    st = stats(recs, n)
    if not st:
        ax.text(0.5, 0.5, 'no data', ha='center', va='center', fontsize=8)
        return
    Gs = {m: v[0] for m, v in st.items()}
    es = {m: v[1] for m, v in st.items()}
    gs = {m: v[2] for m, v in st.items()}
    mG, me, mg = min(Gs, key=Gs.get), min(es, key=es.get), min(gs, key=gs.get)
    keep = []
    for m in (mG, me, mg):
        if m not in keep:
            keep.append(m)
    xs = np.concatenate([np.log10(st[m][3] / n) for m in keep])
    lo, hi = np.percentile(xs, [0.5, 99.5])
    grid = np.linspace(min(lo - 0.15, xs.min()), max(hi + 0.15, xs.max()), 300)
    colors = {mG: '#1f4e79', me: '#b45309', mg: '#7f1d1d'}
    LSQ = {'G': ':', 'e': '--', 'g': '-'}      # by quantity, not by machine
    tags = []
    for q, m in (('G', mG), ('e', me), ('g', mg)):
        tags.append('%s%s' % (q, '' if m == mG else ''))
    for mn in keep:
        G, e, g, c = st[mn]
        z = np.log10(c / n)
        # the density exists only on the observed range: a KDE drawn beyond min/max would put mass
        # where no tape was measured, and g (the maximum) would no longer sit at the curve's end
        gmn = np.linspace(z.min(), z.max(), 300)
        if z.std() < 1e-9:
            d = np.exp(-((gmn - z[0]) / 0.02) ** 2)
        else:
            d = gaussian_kde(z, bw_method=0.35)(gmn)
        d = d / d.max() * 0.7
        ax.fill_between(gmn, 0, d, color=colors[mn], alpha=0.18, lw=0)
        got = [q for q, m in (('G', mG), ('e', me), ('g', mg)) if m == mn]
        lab = '%s (%s)' % (mn, ', '.join(got))
        ax.plot(gmn, d, color=colors[mn], lw=1.0, label=lab)
    # the three values: merge those that coincide, one horizontal rule, labels above
    raw = [(np.log10(Gs[mG] / n), 'v', mG, 'G'), (np.log10(es[me] / n), 'o', me, 'e'),
           (np.log10(gs[mg] / n), '^', mg, 'g')]
    groups = []
    for x, mk, mm, q in raw:
        placed = False
        for g in groups:
            if abs(x - g['x']) < 0.01:
                g['qs'].append(q)
                g['machines'].add(mm)
                placed = True
                break
        if not placed:
            groups.append({'x': x, 'qs': [q], 'mk': mk, 'machines': {mm}})
    order = r'$\Theta(n)$' if rule == 'monotone' else r'$O(n)$'
    for g in groups:
        mm = sorted(g['machines'])[0]
        merged = len(g['qs']) > 1
        ls = LSQ[g['qs'][0]] if not merged else '-.'
        ax.axvline(g['x'], color=colors[mm], ls=ls, lw=1.0, alpha=0.9)
        ax.plot([g['x']], [0.99], g['mk'] if not merged else 'D', ms=5, color=colors[mm],
                markeredgecolor='k', markeredgewidth=0.3, zorder=5)
        ax.annotate('%s = %.2f$n$, %s' % ('='.join(g['qs']), 10 ** g['x'], order),
                    xy=(g['x'], 1.01), xytext=(0, 1), textcoords='offset points',
                    fontsize=5.5, rotation=0, va='bottom', ha='center', color=colors[mm],
                    bbox=dict(boxstyle='round,pad=0.15', fc='white', ec='none', alpha=0.85))
    if title:
        NAME = {'BAL>m/2': 'BAL_m/2', 'BAL>m/4': 'BAL_m/4', 'NBLK=1': 'NBLK_1',
                'NBLK>2': 'NBLK_gt2', 'NBLK>3': 'NBLK_gt3'}
        ax.set_title(NAME.get(rule, rule), fontsize=8, loc='left')
    ax.set_yticks([])
    ax.set_ylim(0, 1.05)
    for s in ('top', 'right', 'left'):
        ax.spines[s].set_visible(False)
    ax.tick_params(labelsize=6)


def one(rule, n):
    os.makedirs(FIGS, exist_ok=True)
    recs = [r for r in load() if r['rule'] == rule and r['n'] == n]
    fig, ax = plt.subplots(figsize=(7.0, 2.4))
    panel(ax, rule, recs, n, title=False)
    ax.set_xlabel(r'$\log_{10}(\mathrm{reads}/n)$, $n=%d$' % n, fontsize=8)
    ax.legend(fontsize=6, frameon=False, ncol=1, loc='upper left',
              bbox_to_anchor=(1.005, 1.02), borderaxespad=0.0)
    fig.suptitle(rule, fontsize=9, x=0.01, ha='left')
    fig.tight_layout()
    out = os.path.join(FIGS, '%s.pdf' % rule.replace('>', 'gt').replace('=', 'eq').replace('/', 's'))
    fig.savefig(out)
    fig.savefig(out.replace('.pdf', '.png'), dpi=170)
    plt.close(fig)
    return out


def atlas(n=None):
    recs = load()
    rules = sorted({r['rule'] for r in recs})
    if not rules:
        return None
    cols = 4
    rows = (len(rules) + cols - 1) // cols
    fig, axes = plt.subplots(rows, cols, figsize=(11, 2.1 * rows))
    axes = np.atleast_1d(axes).ravel()
    for i, rule in enumerate(rules):
        ns = [r['n'] for r in recs if r['rule'] == rule and r['reads']]
        nn = max(ns) if ns else (n or 4096)
        panel(axes[i], rule, [r for r in recs if r['rule'] == rule], nn)
    for j in range(len(rules), len(axes)):
        axes[j].axis('off')
    fig.tight_layout()
    out = os.path.join(HERE, 'mp_atlas.pdf')
    fig.savefig(out)
    fig.savefig(out.replace('.pdf', '.png'), dpi=170)
    plt.close(fig)
    return out


def quad(rules, n, out):
    """A 2x2 page: one panel per rule, the four mechanisms in one float."""
    recs = load()
    fig, axes = plt.subplots(2, 2, figsize=(11, 4.6))
    axes = axes.ravel()
    for i, rule in enumerate(rules[:4]):
        panel(axes[i], rule, [r for r in recs if r['rule'] == rule], n)
        axes[i].set_xlabel(r'$\log_{10}(\mathrm{reads}/n)$', fontsize=7)
    handles, labels = [], []
    for ax in axes:
        h, l = ax.get_legend_handles_labels()
        for hh, ll in zip(h, l):
            if ll not in labels:
                handles.append(hh); labels.append(ll)
    fig.legend(handles, labels, fontsize=6, frameon=False, ncol=4,
               loc='lower center', bbox_to_anchor=(0.5, -0.02))
    fig.tight_layout(rect=(0, 0.06, 1, 1))
    fig.savefig(out)
    fig.savefig(out.replace('.pdf', '.png'), dpi=170)
    plt.close(fig)
    return out


if __name__ == '__main__':
    if sys.argv[1] == 'quad':
        rules = sys.argv[2].split(',')
        n = int(sys.argv[3]) if len(sys.argv) > 3 else 4096
        out = os.path.join(HERE, 'figs', 'quad.pdf')
        os.makedirs(FIG_DIR, exist_ok=True)
        print(quad(rules, n, out))
    elif sys.argv[1] == 'rule':
        rule = sys.argv[2]
        n = int(sys.argv[3]) if len(sys.argv) > 3 else 4096
        print(one(rule, n))
    else:
        n = int(sys.argv[2]) if len(sys.argv) > 2 else 4096
        print(atlas(n))
