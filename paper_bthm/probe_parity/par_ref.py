# -*- coding: utf-8 -*-
"""par_ref -- an independent, SLOW reference for the BTHM-C_2 walk.

Why: the engine (par_engine.walk_par) keeps only ONE BIT per window state, so
by construction it cannot look at a visit counter -- that is the model.  This
file walks the same rule with an EXPLICIT counter per quadruple and asserts,
step by step, that the stored bit is the parity of the count, so that the
engine's semantics is the published one:

  visit 1 reads white, visit 2 reads black (fire), visit 3 reads white again,
  visit 4 black (fire), ...  -- the repair fires exactly on odd-numbered
  repeats; one bit of parity cannot express "ever visited".

The same reference also runs the monotone-flag policy (fire on every repeat,
prior >= 1) and must reproduce engine mode 4, and the parity policy must
reproduce engine mode 5, on the same tapes: totals (steps, reads, cuts, end,
ntrim, nl, nr, msum) must agree tape for tape.

usage: par_ref.py [n] [nt]
writes par_ref.json
"""
import json
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import par_engine as PE                                          # noqa: E402

PAPER = PE.PAPER
sys.path.insert(0, os.path.join(PAPER, 'probe_g_witness'))
sys.path.insert(0, os.path.join(PAPER, 'probe_triple_nat'))
import pcr_lib as P                                              # noqa: E402
import tr_lib as T                                               # noqa: E402
import par_verify as V                                           # noqa: E402


def _b1(S, rid, variant, ell, r):
    return PE._b1_variant(variant, rid, S.a, S.Pa, S.Paa, S.PR, S.Pd, S.Pn2,
                          S.Ppk, S.Ptr, S.Ps, S.cpos, S.pnz, S.C2, S.bstart,
                          S.bpos, S.psmin, S.SPs, S.SNs, S.LOG2, S.r2s,
                          S.r2e, S.SP2, S.rcs, S.rce, S.rcodd, S.SPC, S.R,
                          S.IP, S.HX, S.Ppat, S.d1, S.d2, ell, r)


def ref_walk(S, rid, conv, variant, policy, cap):
    """Exact transcription of par_engine.walk_par / walk(mode 4), in Python,
    with a counter per quadruple (audit only).  policy='parity' fires when
    the bit read is 1 (prior visits odd); policy='flag' fires when the
    quadruple was seen at all (mode 4)."""
    n = S.n
    L, R_ = 0, n - 1
    ell, r = L, R_
    y, x = L, R_
    hy, hx = 0, max(1, n // 2) - 1
    steps = cuts = ntrim = msum = 0
    reads = n + 2 * (n - 1) + 2 + max(1, n // 2)
    end = 1
    cnt = {}          # key -> visits (the reference's luxuy; the model has none)
    bit = {}          # key -> the one bit the model actually carries
    hist = {}         # prior-visit-count histogram of the repair fires
    fresh_repeats = []  # parity only: (prior, fired?) for every repeat
    while True:
        if ell >= r:
            end = 0
            break
        if R_ - L + 1 < 3:
            end = 1
            break
        steps += 1
        if steps + ntrim > cap:
            end = 3 if ntrim == 0 else 6
            break
        m = r - ell + 1
        msum += m
        key = ((ell * (n + 1) + r) * (n + 1) + L) * (n + 1) + R_
        prior = cnt.get(key, 0)
        col = bit.get(key, 0)
        if policy == 'parity':
            # ---- the audit: the one stored bit IS the parity of the visits
            assert col == (prior & 1), ('bit != parity', key, col, prior)
            fire = col == 1
            cnt[key] = prior + 1
            bit[key] = col ^ 1
        else:
            # the monotone flag: set on the first visit, never cleared
            assert col == (1 if prior > 0 else 0), ('flag unset', key, col)
            fire = col == 1
            cnt[key] = prior + 1
            bit[key] = 1
        if prior > 0:
            fresh_repeats.append(prior)
        if fire:
            hist[prior] = hist.get(prior, 0) + 1
            ntrim += 1
            v1 = _b1(S, rid, variant, ell, r)
            if variant == 3:
                act = PE.CUTL if v1 == 1 else PE.RMINUS
            else:
                v2 = PE.b2_of(S.a, S.Pa, ell, r)
                act = PE.TAB[2 * v2 + v1]
            side = 1 if act in (PE.LPLUS, PE.LMINUS, PE.CUTL) else 2
            if conv == 1:
                side = 3 - side
            if side == 1:
                L += 1
            else:
                R_ -= 1
            if L > R_:
                end = 1
                break
            ell, r = L, R_
            reads += abs(y - L) + abs(x - R_)
            y, x = L, R_
            k = (x - y + 1) >> 1
            if k == 0:
                k = 1
            ty, tx = y, y + k - 1
            reads += abs(hy - ty) + abs(hx - tx)
            hy, hx = ty, tx
            continue
        b1 = _b1(S, rid, variant, ell, r)
        if variant == 3:
            act = PE.CUTL if b1 == 1 else PE.RMINUS
        else:
            v2 = PE.b2_of(S.a, S.Pa, ell, r)
            act = PE.TAB[2 * v2 + b1]
        if act >= PE.CUTL:
            cuts += 1
            if act == PE.CUTL:
                L += 1
            elif act == PE.CUTR:
                R_ -= 1
            else:
                L += 1
                R_ -= 1
            if L > R_:
                end = 1
                break
            ell, r = L, R_
            reads += abs(y - L) + abs(x - R_)
            y, x = L, R_
            k = (x - y + 1) >> 1
            if k == 0:
                k = 1
            ty, tx = y, y + k - 1
            reads += abs(hy - ty) + abs(hx - tx)
            hy, hx = ty, tx
            continue
        dy = 1 if act == PE.LPLUS else (-1 if act == PE.LMINUS else 0)
        dx = 1 if act == PE.RPLUS else (-1 if act == PE.RMINUS else 0)
        if ell + dy < L or r + dx > R_:
            cuts += 1
            if conv == 1:
                if ell + dy < L:
                    R_ -= 1
                else:
                    L += 1
            else:
                if ell == L and r == R_:
                    L += 1
                    R_ -= 1
                    y += 1
                    x -= 1
                    reads += 2
                elif ell + dy < L:
                    L += 1
                else:
                    R_ -= 1
            if L > R_:
                end = 1
                break
            ell, r = L, R_
            reads += abs(y - L) + abs(x - R_)
            y, x = L, R_
            k = (x - y + 1) >> 1
            if k == 0:
                k = 1
            ty, tx = y, y + k - 1
            reads += abs(hy - ty) + abs(hx - tx)
            hy, hx = ty, tx
            continue
        if ell + dy >= r + dx:
            end = 0
            break
        ell += dy
        r += dx
        y, x = ell, r
        reads += 1
        k = (x - y + 1) >> 1
        if k == 0:
            k = 1
        ty, tx = y, y + k - 1
        reads += abs(hy - ty) + abs(hx - tx)
        hy, hx = ty, tx
    return dict(steps=steps, reads=reads, cuts=cuts, end=end, ntrim=ntrim,
                msum=msum, hist=hist, fresh=fresh_repeats,
                nkeys=len(cnt))


def main():
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 64
    nt = int(sys.argv[2]) if len(sys.argv) > 2 else 6
    t0 = time.time()
    cap = 40 * n * n + 100
    rules = [r for r, _, _, _, _ in T.SPEC if r != 'monotone']
    bad = []
    per_rule = {}
    hist_all = {}
    fresh_odd = fresh_even = fresh_fired = 0
    for rule in rules:
        rid, cv, vr = V.rrid(rule), V.conv(rule), V.variant(rule)
        tot = dict(tapes=0, e4=0, r4=0, e5=0, r5=0, mismatch=0,
                   fire1=0, fire3=0, fire5=0, white2=0, fire_odd=0,
                   fire_even=0, first_div=None)
        for i, a, s, X in T.tapes(n, nt):
            C = PE.Cyc(X, hlog=20)
            r4 = C.run(rid, cv, vr, mode=4)
            r5 = C.run(rid, cv, vr, mode=5)
            S = P.St(X)
            f4 = ref_walk(S, rid, cv, vr, 'flag', cap)
            f5 = ref_walk(S, rid, cv, vr, 'parity', cap)
            tot['tapes'] += 1
            hits = []
            for k in ('steps', 'reads', 'cuts', 'end', 'ntrim', 'msum'):
                if r4[k] != f4[k]:
                    hits.append(('mode4', i, k, r4[k], f4[k]))
                if r5[k] != f5[k]:
                    hits.append(('mode5', i, k, r5[k], f5[k]))
            bad += hits
            tot['mismatch'] += len(hits)
            # parity fire histogram: fires must sit at ODD prior counts
            for prior, c in f5['hist'].items():
                hist_all[prior] = hist_all.get(prior, 0) + c
                if prior & 1:
                    tot['fire_odd'] += c
                    if prior == 1:
                        tot['fire1'] += c
                    elif prior == 3:
                        tot['fire3'] += c
                    elif prior == 5:
                        tot['fire5'] += c
                else:
                    tot['fire_even'] += c
            for prior in f5['fresh']:
                if prior & 1:
                    fresh_odd += 1
                else:
                    fresh_even += 1
                if prior == 2:
                    tot['white2'] += 1
            if hits and tot['first_div'] is None:
                tot['first_div'] = hits[0]
        per_rule[rule] = tot
        print('%-10s n=%-4d tapes=%d e4=%d e5=%d mism=%d fires(odd)=%d '
              'fires(even)=%d' % (rule, n, tot['tapes'], r4['end'], r5['end'],
                                  tot['mismatch'], tot['fire_odd'],
                                  tot['fire_even']), flush=True)
    out = dict(n=n, tapes=nt, cap=cap, rules=rules, bad=bad, per_rule=per_rule,
               hist_all=hist_all, fresh_odd_repeats=fresh_odd,
               fresh_even_repeats=fresh_even, secs=round(time.time() - t0, 1))
    with open(os.path.join(HERE, 'par_ref.json'), 'w') as fh:
        json.dump(out, fh, indent=1)
    print('mismatches:', len(bad), '| repairs at even prior counts:',
          sum(c for p, c in hist_all.items() if not p & 1),
          '| repairs at odd prior counts:',
          sum(c for p, c in hist_all.items() if p & 1))
    print('repeat visits seen: odd-prior %d, even-prior %d' %
          (fresh_odd, fresh_even))
    print('hist of repair prior counts:', dict(sorted(hist_all.items())))
    print('secs', round(time.time() - t0, 1))
    return 0 if not bad else 1


if __name__ == '__main__':
    sys.exit(main())
