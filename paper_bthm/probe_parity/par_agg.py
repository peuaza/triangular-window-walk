# -*- coding: utf-8 -*-
"""par_agg -- par_cells.json and the per-rule table from the per-cell JSONs.

Every number is read off the per-tape records; nothing is modelled.  The
monotone-flag side is NOT recomputed: it is taken from probe_trim/cut_cells.json
(the stored BTHM-C measurement), as the comparison protocol requires.

Two tape-level checks are made and recorded:
  * every BTHM-0 field of every one of the 18,000 tapes must equal the stored
    probe_trim/cut_<rule>_<n>.json row exactly (the no-memory walk is untouched);
  * every parity dynamical field must equal the stored monotone-flag field:
    the two models must agree tape for tape, not merely on the aggregates.

writes par_cells.json, par_table.md
"""
import json
import os
import statistics as st

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
TRIM = os.path.join(os.path.dirname(HERE), 'probe_trim')
NS = (512, 1024, 2048)
B0_FIELDS = ('end0', 'steps0', 'reads0', 'cyc_start', 'cyc_len', 'width',
             'cuts0')
# parity field -> stored monotone-flag field (cut_run.py row layout)
PF = (('endP', 'endC'), ('stepsP', 'stepsC'), ('readsP', 'readsC'),
      ('ntrim', 'ntrim'), ('cutsP', 'cutsC'), ('nl', 'nl'), ('nr', 'nr'),
      ('nused', 'nused'))


def safe(name):
    return (name.replace('>', 'g').replace('=', 'e').replace('/', 'd')
            .replace('.', 'p'))


def med(xs):
    xs = [x for x in xs if x is not None]
    return float(st.median(xs)) if xs else None


def fits(ns, ys):
    pts = [(n, y) for n, y in zip(ns, ys) if y is not None and y > 0]
    if len(pts) < 3:
        return None
    a = np.polyfit(np.log([p[0] for p in pts]), np.log([p[1] for p in pts]), 1)
    return round(float(a[0]), 3)


def cell(rule, n, bad0, badp):
    f = os.path.join(HERE, 'par_%s_%d.json' % (safe(rule), n))
    if not os.path.exists(f):
        return None
    J = json.load(open(f))
    rows = J['rows']
    nt = len(rows)
    old = json.load(open(os.path.join(TRIM, 'cut_%s_%d.json'
                                      % (safe(rule), n))))['rows']
    # ---- the two tape-level checks
    for k, x in enumerate(rows):
        for fl in B0_FIELDS:
            if int(old[k][fl]) != int(x[fl]):
                bad0.append((rule, n, k, fl, old[k][fl], x[fl]))
        for fa, fb in PF:
            if int(old[k][fb]) != int(x[fa]):
                badp.append((rule, n, k, fa, old[k][fb], x[fa]))
    e0 = [x['end0'] for x in rows]
    eP = [x['endP'] for x in rows]
    cyc = [x for x in rows if x['end0'] == 4]           # BTHM-0 cycle deaths
    met0 = [x for x in rows if x['end0'] == 0]
    metP = [x for x in rows if x['endP'] == 0]
    resc = [x for x in cyc if x['endP'] == 0]
    term = [x for x in cyc if x['endP'] in (0, 1)]
    cap = 40 * n * n + 100
    return dict(
        n=n, tapes=nt, repair=J['repair'], cap=cap, hlog=J['hlog'], mode=5,
        rrid=J['rrid'], variant=J['variant'], conv=J['conv'],
        # (a) BTHM-0, no memory -- the same walk as probe_trim's
        b0_meet=len(met0), b0_meet_frac=round(len(met0) / nt, 4),
        b0_cyc=len(cyc), b0_cyc_frac=round(len(cyc) / nt, 4),
        b0_empty=e0.count(1), b0_empty_frac=round(e0.count(1) / nt, 4),
        b0_capcap=e0.count(3), b0_full=e0.count(5),
        # (b) BTHM-C_2, the parity colour with the same trim repair
        pP_meet=len(metP), pP_meet_frac=round(len(metP) / nt, 4),
        pP_empty=eP.count(1), pP_empty_frac=round(eP.count(1) / nt, 4),
        pP_cyc=eP.count(4), pP_halt=eP.count(7),
        pP_cap=eP.count(3) + eP.count(6), pP_cap_after_repair=eP.count(6),
        pP_full=eP.count(5),
        pP_rescue=len(resc),
        pP_rescue_frac=(round(len(resc) / len(cyc), 4) if cyc else None),
        pP_rescue_any=len(term),
        pP_rescue_any_frac=(round(len(term) / len(cyc), 4) if cyc else None),
        pP_trim_med=med([x['ntrim'] for x in resc]),
        pP_trim_max=max([x['ntrim'] for x in resc], default=None),
        pP_trim_med_all=med([x['ntrim'] for x in rows]),
        pP_trim_max_all=max([x['ntrim'] for x in rows], default=None),
        pP_trim_med_met=med([x['ntrim'] for x in metP]),
        pP_trim_over_n=sum(1 for x in rows if x['ntrim'] > n),
        pP_left=sum(int(x['nl']) for x in rows),
        pP_right=sum(int(x['nr']) for x in rows),
        pP_split_viol=sum(1 for x in rows
                          if int(x['nl']) + int(x['nr']) != int(x['ntrim'])),
        pP_steps_med=med([x['stepsP'] for x in rows]),
        pP_steps_max=max([x['stepsP'] for x in rows], default=None),
        pP_reads_med_all=med([x['readsP'] for x in rows]),
        pP_reads_max=max([x['readsP'] for x in rows], default=None),
        nused_max=max([x['nused'] for x in rows], default=None),
        secs=round(J['secs'], 1),
    )


def main():
    import sys
    sys.path.insert(0, HERE)
    import par_engine                                             # noqa: F401
    from tr_lib import RULENAME
    rules = [r for r in RULENAME if r != 'monotone']
    cutJ = json.load(open(os.path.join(TRIM, 'cut_cells.json'),
                          encoding='utf-8'))['rules']
    bad0, badp = [], []
    cells = {}
    cmp_ = {}
    for rule in rules:
        cells[rule] = {}
        cmp_[rule] = {}
        for n in NS:
            c = cell(rule, n, bad0, badp)
            if c is None:
                continue
            cells[rule][str(n)] = c
            f = cutJ[rule][str(n)]                    # the stored flag cell
            par_meet = c['pP_meet_frac']
            flag_meet = f['bC_meet_frac']
            cmp_[rule][str(n)] = dict(
                par_meet=c['pP_meet'], par_meet_frac=par_meet,
                flag_meet=f['bC_meet'], flag_meet_frac=flag_meet,
                gap=round(par_meet - flag_meet, 4),
                par_rescue=c['pP_rescue'], flag_rescue=f['bC_rescue'],
                par_rescue_frac=c['pP_rescue_frac'],
                flag_rescue_frac=f['bC_rescue_frac'],
                d_rescue=c['pP_rescue'] - f['bC_rescue'],
                par_cap=c['pP_cap'], flag_cap=f['bC_cap'],
                par_full=c['pP_full'], flag_full=f['bC_full'],
                par_trim_med_all=c['pP_trim_med_all'],
                flag_trim_med_all=f['bC_trim_med_all'],
                par_trim_max_all=c['pP_trim_max_all'],
                flag_trim_max_all=f['bC_trim_max_all'],
                par_trim_med=c['pP_trim_med'], flag_trim_med=f['bC_trim_med'],
                par_trim_max=c['pP_trim_max'], flag_trim_max=f['bC_trim_max'],
            )
    # ---- pooled
    pooled = {}
    for n in NS:
        s = dict(tapes=0, b0_meet=0, b0_cyc=0, par_meet=0, flag_meet=0,
                 par_rescue=0, flag_rescue=0, par_cap=0, flag_cap=0,
                 par_full=0, flag_full=0)
        trimsP, trimsF = [], []
        for rule in rules:
            c = cells[rule].get(str(n))
            if c is None:
                continue
            f = cutJ[rule][str(n)]
            s['tapes'] += c['tapes']
            s['b0_meet'] += c['b0_meet']
            s['b0_cyc'] += c['b0_cyc']
            s['par_meet'] += c['pP_meet']
            s['flag_meet'] += f['bC_meet']
            s['par_rescue'] += c['pP_rescue']
            s['flag_rescue'] += f['bC_rescue']
            s['par_cap'] += c['pP_cap']
            s['flag_cap'] += f['bC_cap']
            s['par_full'] += c['pP_full']
            s['flag_full'] += f['bC_full']
            trimsP.append(c['pP_trim_med_all'])
            trimsF.append(f['bC_trim_med_all'])
        p = dict(n=n, **s)
        p['par_meet_frac'] = round(s['par_meet'] / s['tapes'], 4)
        p['flag_meet_frac'] = round(s['flag_meet'] / s['tapes'], 4)
        p['gap'] = round(p['par_meet_frac'] - p['flag_meet_frac'], 4)
        p['par_rescue_frac'] = round(s['par_rescue'] / s['b0_cyc'], 4)
        p['flag_rescue_frac'] = round(s['flag_rescue'] / s['b0_cyc'], 4)
        p['trim_med_of_cell_meds'] = med(trimsP)
        p['flag_trim_med_of_cell_meds'] = med(trimsF)
        pooled[str(n)] = p
    T = sum(pooled[str(n)]['tapes'] for n in NS)
    K = sum(pooled[str(n)]['b0_cyc'] for n in NS)
    M0 = sum(pooled[str(n)]['b0_meet'] for n in NS)
    A = sum(pooled[str(n)]['par_meet'] for n in NS)
    Fm = sum(pooled[str(n)]['flag_meet'] for n in NS)
    R = sum(pooled[str(n)]['par_rescue'] for n in NS)
    RF = sum(pooled[str(n)]['flag_rescue'] for n in NS)
    pooled['all'] = dict(
        tapes=T, b0_meet=M0, b0_cyc=K, par_meet=A, flag_meet=Fm,
        par_meet_frac=round(A / T, 4), flag_meet_frac=round(Fm / T, 4),
        gap=round((A - Fm) / T, 4), par_rescue=R, flag_rescue=RF,
        par_rescue_frac=round(R / K, 4), flag_rescue_frac=round(RF / K, 4),
        n_reported='512-2048 pooled')
    # ---- the gap across n: fit exponent if any gap is nonzero
    gaps = [pooled[str(n)]['gap'] for n in NS]
    if any(g != 0 for g in gaps):
        gap_fit = dict(gaps=gaps, exponent=fits(NS, [abs(g) for g in gaps]))
    else:
        gap_fit = dict(gaps=gaps, exponent=None,
                       note='gap identically 0 at every n; nothing to fit')
    ref = None
    fr = os.path.join(HERE, 'par_ref.json')
    if os.path.exists(fr):
        R_ = json.load(open(fr))
        ref = dict(n=R_['n'], tapes=R_['tapes'],
                   mismatches=len(R_['bad']),
                   fire_hist={str(k): v for k, v in R_['hist_all'].items()},
                   fresh_odd=R_['fresh_odd_repeats'],
                   fresh_even=R_['fresh_even_repeats'])
    out = dict(
        meta=dict(
            model='BTHM-C_2: the parity colour per window state (ell, r, L, R); '
                  'black iff an odd number of visits so far; the bit toggles on '
                  'every visit; the repair (mode-4 tape trim) fires on a BLACK '
                  'window read; action tables unchanged (variant (a), '
                  'colour-blind rules)',
            protocol='tr_lib.tapes(n), 300 random tapes, n in [512, 1024, 2048]',
            engine='probe_parity/par_engine.py: par_engine is cut_engine.py '
                   'with walk_par added (mode 5); modes 0-4 byte-identical',
            flag_source='probe_trim/cut_cells.json (not recomputed)',
            cap='40 n^2 + 100, unchanged',
        ),
        rules=cells, comparison=cmp_, pooled=pooled, gap_fit=gap_fit,
        verification=dict(
            b0_tapes_checked=sum(c['tapes'] for r in rules
                                 for c in cells[r].values()),
            b0_field_mismatches=len(bad0),
            b0_first_mismatch=(bad0[0] if bad0 else None),
            par_vs_flag_field_mismatches=len(badp),
            par_vs_flag_first_mismatch=(badp[0] if badp else None),
            ref=(ref if ref else 'par_ref.json missing'),
            cell_aggregate_mismatches=sum(
                1 for r in rules for k, c in cmp_[r].items()
                for a, b in (('par_meet', 'flag_meet'),
                             ('par_rescue', 'flag_rescue'),
                             ('par_cap', 'flag_cap'), ('par_full', 'flag_full'),
                             ('par_trim_max_all', 'flag_trim_max_all'),)
                if c[a] != c[b]),
        ),
    )
    with open(os.path.join(HERE, 'par_cells.json'), 'w') as fh:
        json.dump(out, fh, indent=1)
    # ---- the table
    L = []
    L.append('| rule | n | par meet | flag meet | gap (pp) | par resc | '
             'flag resc | trim med | trim max | cap P | cap F |')
    L.append('|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|')

    def fres(x):
        return ('%.1f%%' % (100 * x)) if x is not None else '\u2013'

    for rule in rules:
        for n in NS:
            c = cmp_[rule].get(str(n))
            if c is None:
                L.append('| %s | %d | - | - | - | - | - | - | - | - | - |'
                         % (rule, n))
                continue
            L.append('| %s | %d | %.1f%% | %.1f%% | %+.1f | %s | %s | '
                     '%g | %g | %d | %d |'
                     % (rule, n, 100 * c['par_meet_frac'],
                        100 * c['flag_meet_frac'], 100 * c['gap'],
                        fres(c['par_rescue_frac']), fres(c['flag_rescue_frac']),
                        c['par_trim_med_all'] if c['par_trim_med_all']
                        is not None else -1,
                        c['par_trim_max_all'] if c['par_trim_max_all']
                        is not None else -1,
                        c['par_cap'], c['flag_cap']))
    L.append('| **pooled** | %s | %.1f%% | %.1f%% | %+.1f | %.1f%% | %.1f%% | '
             '- | - | %d | %d |'
             % ('all',
                100 * pooled['all']['par_meet_frac'],
                100 * pooled['all']['flag_meet_frac'],
                100 * pooled['all']['gap'],
                100 * pooled['all']['par_rescue_frac'],
                100 * pooled['all']['flag_rescue_frac'],
                sum(c['pP_cap'] for r in rules for c in cells[r].values()),
                sum(f['bC_cap'] for r in rules for f in cutJ[r].values())))
    open(os.path.join(HERE, 'par_table.md'), 'w').write('\n'.join(L) + '\n')
    print('\n'.join(L))
    print()
    print('b0 tapes checked %d, mismatches %d | parity-vs-flag mismatches %d'
          % (out['verification']['b0_tapes_checked'],
             len(bad0), len(badp)))
    for n in NS:
        p = pooled[str(n)]
        print('pooled n=%-5d tapes=%d meet0=%d cap0=%d | par meet=%d (%.4f) '
              'flag meet=%d (%.4f) gap=%+.4f | rescue %d/%d vs %d/%d | '
              'cap P/F=%d/%d full P/F=%d/%d'
              % (n, p['tapes'], p['b0_meet'], p['b0_cyc'], p['par_meet'],
                 p['par_meet_frac'], p['flag_meet'], p['flag_meet_frac'],
                 p['gap'], p['par_rescue'], p['b0_cyc'], p['flag_rescue'],
                 p['b0_cyc'], p['par_cap'], p['flag_cap'], p['par_full'],
                 p['flag_full']))
    p = pooled['all']
    print('pooled all tapes=%d meet0=%d (%.4f) cap0=%d meet P=%d (%.4f) '
          'F=%d (%.4f) gap=%+.4f rescue P=%d F=%d of %d (%.4f / %.4f)'
          % (p['tapes'], p['b0_meet'], p['b0_meet'] / p['tapes'], p['b0_cyc'],
             p['par_meet'], p['par_meet_frac'], p['flag_meet'],
             p['flag_meet_frac'], p['gap'], p['par_rescue'], p['flag_rescue'],
             p['b0_cyc'], p['par_rescue_frac'], p['flag_rescue_frac']))
    print('gap by n:', out['gap_fit'])
    print('ref:', ref)


if __name__ == '__main__':
    main()
