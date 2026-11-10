# -*- coding: utf-8 -*-
"""k3_run -- K = 300 meeting tapes per (rule, n) cell, band family.

WHY
---
The band family's e exponents (KCHG 1.29, BAL>m/2 1.57, BAL>m/4 1.61,
NBLK>3 1.69) rest on K = 60-200 meeting tapes at the top sizes.  If the cost
distribution over meeting tapes is heavy-tailed, a small K understates the mean
and the bias grows with n -- which is exactly the shape that fakes a
sub-quadratic exponent.  This probe holds K fixed at 300 on all four sizes.

PROTOCOL -- identical to MASTER_FINAL (probe_table_eg / probe_triple_nat)
-------------------------------------------------------------------------
tapes   T_n: two independent uniform h-subsets, h = n//2, one per column,
        rng = default_rng(5000 + n + 100003 i).  Same bytes as
        eu_lib.make_tape / teg_lib.tape / tr_lib.tapes.
machine tr_lib.run(rule, a, s), unmodified: hmach3.H3 (B=64) for the band
        rules, pcr_lib.Runner for the rest; charge INC.  The machine label and
        the charge string are copied out of tr_lib.spec_of into every cell.
meet    end_code == 0 (a numeric delta in [0,1]).
e       mean reads over the meeting tapes.  This is the MASTER_FINAL
        definition: KCHG@512 reproduces 6191.7 -> 6192 from the same cell.
g_K     max reads over the FIRST 300 meeting tapes of the cell (fewer if the
        cell is low-K, and then the cell is flagged).

Caps: 300 meetings, 12000 draws, 2400 s wall -- whichever comes first.  A cell
that stops on a cap writes what it has and is flagged lowK.

USAGE
  python k3_run.py cell RULE N
  python k3_run.py grid [--jobs 3] [--only a,b] [--ns 512,1024]
  python k3_run.py bench N [K]
"""
import json
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
PAPER = os.path.dirname(HERE)
for _d in ('probe_table_eg', 'probe_triple_nat'):
    _p = os.path.join(PAPER, _d)
    if _p not in sys.path:
        sys.path.insert(0, _p)
if HERE not in sys.path:
    sys.path.insert(0, HERE)

import teg_lib as G                                              # noqa: E402
import tr_lib as T                                               # noqa: E402

RULES = ['KCHG', 'BAL>m/2', 'BAL>m/4', 'NBLK=1', 'NBLK>2', 'NBLK>3',
         'AMAX', 'AMIN', 'MAXP', 'LAS', 'oddrun']
NS = (512, 1024, 2048, 4096)

TARGET = 300
CAP = 12000
SECONDS = 2400.0
BATCH = 500


def safe(name):
    return (name.replace('>', 'g').replace('=', 'e').replace('/', 'd')
            .replace('.', 'p'))


def logpath(rule, n):
    return os.path.join(HERE, 'k3_%s_%d.log' % (safe(rule), n))


def jsonpath(rule, n):
    return os.path.join(HERE, 'k3_%s_%d.json' % (safe(rule), n))


def log(rule, n, msg):
    line = '[%s] %-8s n=%-5d %s' % (time.strftime('%H:%M:%S'), rule, n, msg)
    print(line, flush=True)
    with open(logpath(rule, n), 'a') as fh:
        fh.write(line + '\n')


def write_cell(rule, n, reads, meet_i, ndraw, t0, stop):
    e = float(np.mean(reads)) if reads else float('nan')
    k_used = min(TARGET, len(reads))
    g = float(max(reads[:k_used])) if reads else float('nan')
    tag, sub, mach, chg = T.spec_of(rule)
    d = dict(rule=rule, n=n, machine=mach, charge=chg, spec_tag=tag,
             spec_sub=sub, nmeet=len(reads), ndraw=ndraw, K_used=k_used,
             e=e, g_K=g, reads=[int(x) for x in reads],
             meet_i=[int(x) for x in meet_i],
             meet_rate=(len(reads) / ndraw if ndraw else 0.0),
             lowK=bool(len(reads) < TARGET),
             stop=stop, secs=round(time.time() - t0, 1),
             seed_formula='5000+n+100003*i')
    with open(jsonpath(rule, n), 'w') as fh:
        json.dump(d, fh)
    d.pop('reads')
    d.pop('meet_i')
    return d


def run_cell(rule, n, target=TARGET, cap=CAP, seconds=SECONDS):
    t0 = time.time()
    reads, meet_i = [], []
    i = 0
    stop = 'target'
    while True:
        if len(reads) >= target:
            stop = 'target'
            break
        if i >= cap:
            stop = 'draws_cap'
            break
        if time.time() - t0 > seconds:
            stop = 'wall_cap'
            break
        nb = min(BATCH, cap - i)
        for _ in range(nb):
            if len(reads) >= target:
                break
            a, s = G.tape(n, i)
            r = T.run(rule, a, s)
            if G.end_code(r.get('end')) == 0:
                reads.append(int(round(float(r['reads']))))
                meet_i.append(i)
            i += 1
            if time.time() - t0 > seconds:
                break
        d = write_cell(rule, n, reads, meet_i, i, t0, stop)
        log(rule, n, '%5d draws %4d meet  e=%.0f  g_K=%.0f  %.0fs  %s'
            % (i, len(reads), d['e'], d['g_K'], time.time() - t0, stop))
        if stop != 'target' or len(reads) >= target:
            break
    d = write_cell(rule, n, reads, meet_i, i, t0, stop)
    log(rule, n, 'DONE %d draws %d meet e=%.0f g_K=%.0f %s [%.0fs]'
        % (i, len(reads), d['e'], d['g_K'],
           'GATE-OK' if len(reads) >= target else 'low-K', d['secs']))
    return d


# ------------------------------------------------------------------- bench
def bench(n, k=6, rules=None):
    rules = RULES if rules is None else rules
    out = []
    for rule in rules:
        t0 = time.time()
        mt = 0
        readsum = 0
        for i in range(k):
            a, s = G.tape(n, i)
            r = T.run(rule, a, s)
            if G.end_code(r.get('end')) == 0:
                mt += 1
            readsum += int(round(float(r['reads'])))
        dt = (time.time() - t0) / k
        out.append((rule, n, dt, mt, readsum / k))
        print('%-8s n=%-5d %.3f s/draw  meet %d/%d  meanreads %.0f'
              % (rule, n, dt, mt, k, readsum / k), flush=True)
    return out


# -------------------------------------------------------------------- grid
def _task(a):
    rule, n = a
    try:
        d = run_cell(rule, n)
        return (rule, n, d['nmeet'], d['ndraw'], d['e'], d['g_K'],
                d['lowK'], d['secs'], None)
    except Exception:                                            # noqa: BLE001
        import traceback
        return (rule, n, 0, 0, 0, 0, True, 0, traceback.format_exc())


def grid(rules, ns, jobs):
    from multiprocessing import Pool
    tasks = [(r, n) for r in rules for n in ns]
    # longest first: big n and the band rules decide the makespan
    tasks.sort(key=lambda t: (-t[1], t[0]))
    print('grid: %d cells, jobs=%d, target=%d, cap=%d, wall=%ds'
          % (len(tasks), jobs, TARGET, CAP, SECONDS), flush=True)
    t0 = time.time()
    done = 0
    with Pool(jobs, maxtasksperchild=1) as pool:
        for rule, n, meet, ndraw, e, g, lowK, secs, err in \
                pool.imap_unordered(_task, tasks):
            done += 1
            if err:
                print('!! %s n=%d FAILED\n%s' % (rule, n, err), flush=True)
            else:
                print('[%2d/%2d] %-8s n=%-5d meet=%-4d draws=%-6d e=%-12.0f '
                      'g_K=%-14.0f %s  cell %.0fs  total %.0fs'
                      % (done, len(tasks), rule, n, meet, ndraw, e, g,
                         'low-K' if lowK else 'OK', secs, time.time() - t0),
                      flush=True)
    print('grid done %.0fs' % (time.time() - t0), flush=True)


if __name__ == '__main__':
    cmd = sys.argv[1] if len(sys.argv) > 1 else 'grid'
    if cmd == 'cell':
        run_cell(sys.argv[2], int(sys.argv[3]))
    elif cmd == 'bench':
        n = int(sys.argv[2])
        k = int(sys.argv[3]) if len(sys.argv) > 3 else 6
        rs = sys.argv[4].split(',') if len(sys.argv) > 4 else None
        bench(n, k, rs)
    else:
        jobs = 3
        rules = RULES
        ns = list(NS)
        for a in sys.argv[2:]:
            if a.startswith('--jobs'):
                jobs = int(a.split('=')[1])
            elif a.startswith('--only'):
                rules = a.split('=')[1].split(',')
            elif a.startswith('--ns'):
                ns = [int(x) for x in a.split('=')[1].split(',')]
        grid(rules, ns, jobs)
