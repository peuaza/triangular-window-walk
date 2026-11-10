# -*- coding: utf-8 -*-
"""ce_step1 -- per-epoch failure rate of the CES outcome rule.

For every tape and every cut epoch of the TRUE walk, compare
    pred = 'D' if b2@start == 0 else ('R' if b1@start else 'L')
with the true epoch outcome.  Report the rate per epoch, the distribution of
the first failing epoch, and the rate restricted to M0 > W (the part CES
shortcuts).

usage: ce_step1.py <rule> <n> <nt> [seed] [W]
Writes ce_fail_<rule>_<n>.json
"""
import json
import os
import sys
from collections import Counter

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
for _p in (HERE, os.path.join(ROOT, 'probe_e_uniform'),
           os.path.join(ROOT, 'probe_kappa_pred'),
           os.path.join(ROOT, 'probe_ologn2')):
    if _p not in sys.path:
        sys.path.insert(0, _p)

import ce_lib as CL                                              # noqa: E402
import eu_lib as EL                                              # noqa: E402


def buckets():
    return [0, 2, 3, 4, 6, 9, 13, 17, 25, 33, 49, 65, 97, 129, 193,
            257, 385, 513, 1025, 2049, 4097, 10 ** 9]


def run(rule, n, nt, seed=0, W=16):
    tapes = []
    for i in range(nt):
        a, s, X = EL.make_tape(n, 5000 + n + 100003 * i)
        tapes.append(np.asarray(a, np.int64))
    recs = []
    tot_ep = tot_fail = 0
    ep_b2_1 = fail_b2_1 = 0
    ep_b2_0 = fail_b2_0 = 0
    short_ep = short_fail = 0
    n_cut_mis = n_term = 0
    term_M0 = Counter()
    term_ln = []
    tapes_any = tapes_any_short = 0
    meet_tapes = meet_fail = 0
    first_idx = []          # first failing epoch index
    first_M0 = []
    first_kind = Counter()
    fail_hist_M0 = Counter()
    ep_hist_M0 = Counter()
    fail_hist_pos = Counter()
    ep_hist_pos = Counter()
    fail_ln = []
    b2_1_lens = []
    for a in tapes:
        r = CL.measure_tape(a, rule, W)
        recs.append(dict(why=r['why'], delta=r['delta'], neps=r['neps'],
                         nfail=r['nfail'], first_fail=r['first_fail'],
                         nfail_short=r['nfail_short'], n_short=r['n_short'],
                         ncut_mismatch=r['ncut_mismatch'], nterm=r['nterm']))
        tot_ep += r['neps']; tot_fail += r['nfail']
        n_cut_mis += r['ncut_mismatch']; n_term += r['nterm']
        if r['neps']:
            if r['nfail']:
                tapes_any += 1
                first_idx.append(r['first_fail']['idx'] / r['neps'])
                first_M0.append(r['first_fail']['M0'])
                first_kind[r['first_fail']['kind']] += 1
            if r['nfail_short']:
                tapes_any_short += 1
            if r['why'] == 'meet':
                meet_tapes += 1
                if r['nfail']:
                    meet_fail += 1
        for e in r['epochs']:
            M0 = e['M0']; pos = e['j'] / max(1, r['neps'])
            if e['b2']:
                ep_b2_1 += 1
                b2_1_lens.append(e['ln'])
                if e['fail']:
                    fail_b2_1 += 1
            else:
                ep_b2_0 += 1
                if e['fail']:
                    fail_b2_0 += 1
            if e['terminal'] and M0 > W:
                term_M0[M0] += 1
                term_ln.append(e['ln'])
            if M0 > W:
                short_ep += 1
                if e['fail']:
                    short_fail += 1
            if e['fail']:
                fail_hist_M0[M0] += 1
                fail_hist_pos[int(pos * 10)] += 1
                fail_ln.append(e['ln'])
            ep_hist_M0[M0] += 1
            ep_hist_pos[int(pos * 10)] += 1
    out = dict(rule=rule, n=n, nt=nt, W=W, seed=seed,
               tapes=nt, meet_tapes=meet_tapes, meet_with_fail=meet_fail,
               tot_epochs=tot_ep, tot_fail=tot_fail,
               rate=tot_fail / max(1, tot_ep),
               rate_x_n=tot_fail / max(1, tot_ep) * n,
               ep_b2_1=ep_b2_1, fail_b2_1=fail_b2_1,
               rate_b2_1=fail_b2_1 / max(1, ep_b2_1),
               ep_b2_0=ep_b2_0, fail_b2_0=fail_b2_0,
               short_epochs=short_ep, short_fail=short_fail,
               rate_short=short_fail / max(1, short_ep),
               n_cut_mismatch=n_cut_mis, n_term=n_term,
               term_M0=sorted(term_M0.items()),
               term_ln_mean=float(np.mean(term_ln)) if term_ln else None,
               term_ln_max=int(np.max(term_ln)) if term_ln else None,
               first_kind=dict(first_kind),
               tapes_any_fail=tapes_any, tapes_any_fail_short=tapes_any_short,
               first_idx_frac=first_idx, first_M0=first_M0,
               fail_hist_M0=sorted(fail_hist_M0.items()),
               ep_hist_M0=sorted(ep_hist_M0.items()),
               fail_hist_pos=sorted(fail_hist_pos.items()),
               ep_hist_pos=sorted(ep_hist_pos.items()),
               fail_ln_mean=float(np.mean(fail_ln)) if fail_ln else None,
               fail_ln_max=int(np.max(fail_ln)) if fail_ln else None,
               b2_1_ln_mean=float(np.mean(b2_1_lens)) if b2_1_lens else None,
               recs=recs)
    fn = os.path.join(HERE, 'ce_fail_%s_%d.json' % (rule, n))
    with open(fn, 'w') as f:
        json.dump(out, f, indent=1)
    print('%-8s n=%-5d nt=%-4d meet=%-4d epochs=%-8d fail=%-6d rate=%.2e (rate*n=%.3g)'
          % (rule, n, nt, meet_tapes, tot_ep, tot_fail, out['rate'], out['rate_x_n']))
    print('   b2=1 epochs %d fail %d rate %.2e | b2=0 epochs %d fail %d'
          % (ep_b2_1, fail_b2_1, out['rate_b2_1'], ep_b2_0, fail_b2_0))
    print('   failure modes: cut_mismatch %d (all M) / early_term %d ; of the latter M0>W %d'
          % (n_cut_mis, n_term, short_fail))
    print('   early-terminal M0: %s  len mean %s max %s'
          % (dict(Counter(term_M0).most_common(8)), out['term_ln_mean'], out['term_ln_max']))
    print('   shortcut part (M0>%d): epochs %d fail %d rate %.2e ; tapes any %d/%d'
          % (W, short_ep, short_fail, out['rate_short'], tapes_any_short, nt))
    print('   tapes with any fail %d/%d ; meeting tapes with fail %d/%d'
          % (tapes_any, nt, meet_fail, meet_tapes))
    print('   first-fail kind %s ; M0: %s'
          % (dict(first_kind), dict(Counter(first_M0).most_common(10))))


if __name__ == '__main__':
    rule = sys.argv[1]
    n = int(sys.argv[2])
    nt = int(sys.argv[3]) if len(sys.argv) > 3 else 100
    seed = int(sys.argv[4]) if len(sys.argv) > 4 else 0
    W = int(sys.argv[5]) if len(sys.argv) > 5 else 16
    run(rule, n, nt, seed, W)
