# -*- coding: utf-8 -*-
"""el_cell -- one (rule, n) cell: the epoch trace and the jump machine.

For each drawn tape the walk is run once (el_lib.trace_walk, gated against
pcr_lib.Runner by el_gate.py) and every epoch is recorded.  On the same tape the
jump machine (below) is run for each W2 and charged reads.

THE JUMP MACHINE  (CES's shape, generalised)
--------------------------------------------
State: the region [L,R] and one statistic per rule, O(log n) bits.
  while M > W2:  evaluate (b2,b1) at the region, take the epoch's cut from
                 b2 == 0  -> double cut
                 b2 == 1  -> cutR if b1 else cutL
                 charge: 1 read per cell the cut gives up, + `eval` reads for
                 the predicate (0 if the statistic slides, M if it has to be
                 re-scanned).
  while M <= W2: run the epoch literally, one step at a time, charging 1 read
                 per step (the cell entering the window) + `eval` per step.
The machine never consults the truth; its (delta, end) is compared with the
trace's.  A wrong cut at any epoch sends it down a different branch, so the
exactness count is a property of the machine, not of the trace.

READ MODEL (the record's)
-------------------------
One cell touch costs 1, reads counted with multiplicity, opening pass included.
`eval` is the honest per-query price of the rule's predicate under the
O(log n)-bit state rule: 0 for a window statistic that slides (one counter, one
read per cell given up), M for a range query (max/min over the window, or a
level index) that has to be re-scanned because O(log n) bits cannot hold it.

Usage: python el_cell.py RULE N [NMEET NHALT DRAW_CAP]
Writes el_<rule>_<n>.json
"""
import json
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import el_lib as E
import pcr_lib as P
import tr_bal4

TAB = P.TABLE_NP
CONV = dict(P.CONV)
CONV.setdefault('BAL>m/4', 0)
W2S = (8, 16, 32, 64, 128, 256)
BINS = (1, 2, 3, 4, 5, 9, 17, 33, 65, 257, 1025, 1 << 40)
BNAMES = ('1', '2', '3', '4', '5-8', '9-16', '17-32', '33-64', '65-256',
          '257-1024', '>1024')


def safe(name):
    return (name.replace('>', 'gt').replace('=', 'e').replace('/', 'd')
            .replace('.', 'p'))


def pred_of(rule, b2, b1):
    """The epoch outcome the machine predicts from the two bits at its start.

    conv 'a' (fourteen rules): b2 = 0 -> the move leaves the region at once, and
    the window is the whole region, so it is a double cut.  b2 = 1 -> the window
    shrinks, and the epoch exits on the side the first move pushed from:
    b1 = 0 (r -= 1) -> cut on the left, b1 = 1 (ell += 1) -> cut on the right.
    conv 'c1a' (canonical) cuts the OPPOSITE side, so the side is b1's alone.
    """
    if CONV.get(rule, 0) == 1:
        return 1 if b1 else 2
    if b2 == 0:
        return 3
    return 2 if b1 else 1


def b1_at(S, rule, ell, r):
    if rule in E.NBLSAT:
        return int(E.nblk_b1(S.Pa, S.Ps, S.psmin, S.bstart, S.bpos,
                             E.NBLSAT[rule], ell, r))
    if rule in E.BALTHR:
        return int(tr_bal4.bal_b1(S.a, S.Pa, S.Ps, S.psmin, S.bstart, S.bpos,
                                  E.BALTHR[rule], ell, r))
    return int(P.b1_of(P.RID[rule], S.a, S.Pa, S.Paa, S.PR, S.Pd, S.Pn2, S.Ppk,
                       S.Ptr, S.Ps, S.cpos, S.pnz, S.C2, S.bstart, S.bpos,
                       S.psmin, S.SPs, S.SNs, S.LOG2, S.r2s, S.r2e, S.SP2,
                       S.rcs, S.rce, S.rcodd, S.SPC, S.R, S.IP, S.HX, S.Ppat,
                       S.d1, S.d2, ell, r))


def b2_at(S, ell, r):
    return int(P.b2_of(S.a, S.Pa, ell, r))


def jm_walk(rule, S, n, W2, mode):
    """The jump machine.  Returns (delta, end, reads, epochs, cuts, eg, nscan)."""
    a = S.a
    sv = S.s
    L, R = 0, n - 1
    reads = n                      # opening pass: region sum (+ argmax)
    epochs = 0
    cuts = 0
    eg = 0
    nscan = 0
    lazy = (mode == 'lazy')
    ismax = rule in ('MAXP', 'AMAX')
    # the range queries are over Ps, the prefix walk of the LABEL tape s, not
    # over the a-sum (sum(s) = 0, so Ps[n] = 0 -- using a here would be wrong)
    Ssum = 0
    for j in range(n):
        Ssum += int(sv[j])
    p = -1
    pv = 0
    if lazy:
        qhi = (R + 1) if rule == 'MAXP' else R
        acc = 0
        best = 0
        p = L
        for j in range(L, qhi + 1):
            if j > L:
                acc += int(sv[j - 1])
            if (acc > best) if ismax else (acc < best):
                best = acc
                p = j
        pv = best
    while True:
        M = R - L + 1
        if M < 3:
            return None, 1, reads, epochs, cuts, eg, nscan
        if M <= W2:
            ell, r = L, R
            seen = set()
            while True:
                if ell >= r:
                    return ell / (n - 1.0), 0, reads, epochs, cuts, eg, nscan
                if R - L + 1 < 3:
                    return None, 1, reads, epochs, cuts, eg, nscan
                if (ell, r) in seen:
                    return None, 2, reads, epochs, cuts, eg, nscan
                seen.add((ell, r))
                m = r - ell + 1
                b1 = b1_at(S, rule, ell, r)
                b2 = b2_at(S, ell, r)
                # one read per step (the cell entering the window); the eval is
                # free when the statistic slides, m otherwise
                reads += 1 + (0 if mode == 'sliding' else m)
                eg += 1
                act = TAB[2 * b2 + b1]
                dy = 1 if act == P.LPLUS else (-1 if act == P.LMINUS else 0)
                dx = 1 if act == P.RPLUS else (-1 if act == P.RMINUS else 0)
                if ell + dy < L or r + dx > R:
                    if CONV.get(rule, 0) == 1:
                        if ell + dy < L:
                            R -= 1
                        else:
                            L += 1
                    else:
                        if ell == L and r == R:
                            L += 1
                            R -= 1
                        elif ell + dy < L:
                            L += 1
                        else:
                            R -= 1
                    if L > R:
                        return None, 1, reads, epochs, cuts, eg, nscan
                    ell, r = L, R
                    seen.clear()
                    continue
                if ell + dy >= r + dx:
                    return (ell + dy) / (n - 1.0), 0, reads, epochs, cuts, eg, nscan
                if r + dx < L or ell + dy > R:
                    L += 1
                    R -= 1
                    if L > R:
                        return None, 1, reads, epochs, cuts, eg, nscan
                    ell, r = L, R
                    seen.clear()
                    continue
                ell += dy
                r += dx
        else:
            if lazy:
                qhi = (R + 1) if rule == 'MAXP' else R
                if p < L or p > qhi:
                    reads += qhi - L + 1
                    nscan += 1
                    acc = 0
                    best = 0
                    p = L
                    for j in range(L, qhi + 1):
                        if j > L:
                            acc += int(sv[j - 1])
                        if (acc > best) if ismax else (acc < best):
                            best = acc
                            p = j
                    pv = best
                if rule == 'MAXP':
                    b1 = 1 if pv * pv > M else 0
                elif rule == 'AMAX':
                    b1 = 1 if Ssum > pv else 0
                else:
                    b1 = 1 if Ssum < pv else 0
            else:
                b1 = b1_at(S, rule, L, R)
                if mode == 'scan':
                    reads += M
            b2 = b2_at(S, L, R)
            epochs += 1
            pr = pred_of(rule, b2, b1)
            cuts += 1
            if pr == 3:
                reads += 2
                dw = int(sv[L])
                Ssum -= dw + int(sv[R])
                if lazy:
                    pv -= dw
                L += 1
                R -= 1
            elif pr == 2:
                reads += 1
                Ssum -= int(sv[R])
                R -= 1
            else:
                reads += 1
                dw = int(sv[L])
                Ssum -= dw
                if lazy:
                    pv -= dw
                L += 1


def hist(vals):
    h = [0] * len(BNAMES)
    for v in vals:
        for j in range(len(BINS) - 1):
            if BINS[j] <= v < BINS[j + 1]:
                h[j] += 1
                break
    return h


def run(rule, n, nmeet, nhalt, cap_draw):
    t0 = time.time()
    mode = ('sliding' if E.EVAL[rule] == 'sliding'
            else ('lazy' if E.EVAL[rule] == 'rangemax' else 'scan'))
    taps = []
    i = 0
    nm = nh = 0
    while (nm < nmeet or nh < nhalt) and i < cap_draw:
        X = E.tape(n, i)
        M = E.Mach(X, n)
        t = M.trace(rule)
        meet = (t['end'] == 0)
        keep = (meet and nm < nmeet) or ((not meet) and nh < nhalt)
        if not keep:
            i += 1
            continue
        if meet:
            nm += 1
        else:
            nh += 1
        # ---- per-epoch prediction audit
        epM = t['epM']
        epb2 = t['epb2']
        epb1 = t['epb1']
        epout = t['epout']
        eplen = t['eplen']
        nep = t['nep']
        bad = []
        for k in range(nep):
            if epb2[k] < 0:
                continue
            pred = pred_of(rule, int(epb2[k]), int(epb1[k]))
            if int(epout[k]) != pred:
                bad.append([int(epM[k]), int(pred), int(epout[k]), k,
                            1 if k == nep - 1 else 0, int(eplen[k])])
        # ---- the machine, per W2
        mach = {}
        for W2 in W2S:
            d, en, rd, epc, cu, egs, nsc = jm_walk(rule, M.S, n, W2, mode)
            ok = (int(en) == int(t['end'])) and (
                (d is None and t['delta'] < 0) or
                (d is not None and t['delta'] >= 0 and abs(d - t['delta']) < 1e-12))
            mach[str(W2)] = dict(delta=(-1.0 if d is None else float(d)),
                                 end=int(en), reads=int(rd), epochs=int(epc),
                                 cuts=int(cu), eg=int(egs), nscan=int(nsc),
                                 exact=int(bool(ok)))
        taps.append(dict(i=i, end=int(t['end']),
                         delta=float(t['delta']), steps=int(t['steps']),
                         cuts=int(t['cuts']), inc=int(t['inc']),
                         msum=int(t['msum']), nep=int(nep),
                         termM0=int(epM[-1]) if nep else -1,
                         termout=int(epout[-1]) if nep else -1,
                         termlen=int(eplen[-1]) if nep else -1,
                         mach=mach,
                         bad=bad))
        i += 1
    mean_ep_len = [float(np.mean([x['steps'] / max(1, x['nep']) for x in taps]))]
    d = dict(rule=rule, n=n, mode=mode, nmeet=nm, nhalt=nh, ndraw=i,
             secs=round(time.time() - t0, 1), w2s=list(W2S),
             taps=taps,
             mean_epoch_len=mean_ep_len[0],
             ep_len_hist=hist([x['termlen'] for x in taps]),
             binnames=list(BNAMES))
    p = os.path.join(HERE, 'el_%s_%d.json' % (safe(rule), n))
    with open(p, 'w') as fh:
        json.dump(d, fh)
    return d


def main():
    rule = sys.argv[1]
    n = int(sys.argv[2])
    nmeet = int(sys.argv[3]) if len(sys.argv) > 3 else 200
    nhalt = int(sys.argv[4]) if len(sys.argv) > 4 else 100
    cap = int(sys.argv[5]) if len(sys.argv) > 5 else 60000
    d = run(rule, n, nmeet, nhalt, cap)
    print('%-9s n=%4d mode=%-7s meet=%d halt=%d draws=%d secs=%.0f  '
          'e_old/n(meet)=%.1f' % (
              rule, n, d['mode'], d['nmeet'], d['nhalt'], d['ndraw'],
              d['secs'],
              np.mean([x['inc'] for x in d['taps'] if x['end'] == 0]) / n))


if __name__ == '__main__':
    main()
