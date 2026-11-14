# -*- coding: utf-8 -*-
"""col_check -- the two gates before any colour number is quoted.

A  tape identity / BTHM-0: every mode-0 field (end0, steps0, reads0, cyc_start,
   cyc_len, width, cuts0) of every tape of col_<rule>_<n>.json must equal the
   stored probe_trim/cut_<rule>_<n>.json row (cut_run.py's BTHM-0 run on the
   same tr_lib.tapes ensemble).  This is the same-tape check of COLOUR.md.

B  the kmax=1 control: col_k1_<rule>_<n>.json must reproduce the stored cut
   row on (endC, stepsC, readsC, ntrim, cutsC, nl, nr, nused) tape for tape --
   the colour machinery at kmax=1 with colour pinned to 0 is BTHM-C, and the
   stored file is the independent run of it.

C  structural invariants of the colour runs: every repair fires with the
   stored count exactly kmax (repok == ntrim); nl + nr == ntrim; a tape that
   never closed a period in BTHM-0 ends and reads identically under the colour
   (invariant); no run exceeds n repairs.

writes col_check.json
"""
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import col_engine as E                                          # noqa: E402

TRIM = os.path.join(E.PAPER, 'probe_trim')
RULES = ('pattern', 'runs', 'peak')
NS = (512, 1024, 2048)
B0F = ('end0', 'steps0', 'reads0', 'cyc_start', 'cyc_len', 'width', 'cuts0')
C1F = ('endC', 'stepsC', 'readsC', 'ntrim', 'cutsC', 'nl', 'nr', 'nused')


def safe(name):
    return (name.replace('>', 'g').replace('=', 'e').replace('/', 'd')
            .replace('.', 'p'))


def load(f):
    p = os.path.join(HERE, f)
    return json.load(open(p)) if os.path.exists(p) else None


def stored(rule, n):
    p = os.path.join(TRIM, 'cut_%s_%d.json' % (safe(rule), n))
    if not os.path.exists(p):
        return None
    return {int(r['i']): r for r in json.load(open(p))['rows']}


def main():
    out = {}
    for rule in RULES:
        for n in NS:
            J = load('col_%s_%d.json' % (safe(rule), n))
            if J is None:
                continue
            K = load('col_k1_%s_%d.json' % (safe(rule), n))
            S = stored(rule, n)
            cell = dict(rule=rule, n=n, tapes=J['tapes'], kmax=J['kmax'],
                        secs=round(J['secs'], 1))
            # A: mode 0 vs the stored file
            a_bad = []
            b_bad = []
            c_bad = []
            inv_bad = []
            trim_over_n = 0
            for r in J['rows']:
                i = r['i']
                if S is not None and i in S:
                    s = S[i]
                    d = [f for f in B0F if int(s[f]) != int(r[f])]
                    if d:
                        a_bad.append((i, d, [s[f] for f in d],
                                      [r[f] for f in d]))
                if int(r['repok']) != int(r['ntrim']):
                    c_bad.append((i, 'repok', r['repok'], r['ntrim']))
                if int(r['nl']) + int(r['nr']) != int(r['ntrim']):
                    c_bad.append((i, 'nl+nr', r['nl'], r['nr']))
                if int(r['ntrim']) > n:
                    trim_over_n += 1
                if int(r['end0']) in (0, 1, 3) and (
                        int(r['endC']) != int(r['end0'])
                        or int(r['stepsC']) != int(r['steps0'])
                        or int(r['readsC']) != int(r['reads0'])
                        or int(r['ntrim']) != 0):
                    inv_bad.append((i, r['end0'], r['endC'], r['steps0'],
                                    r['stepsC'], r['reads0'], r['readsC'],
                                    r['ntrim']))
            cell['A_b0_tapes'] = len(J['rows'])
            cell['A_b0_mismatch'] = len(a_bad)
            cell['A_first'] = a_bad[0] if a_bad else None
            # B: the kmax=1 control
            if K is not None and S is not None:
                kb = []
                for r in K['rows']:
                    i = r['i']
                    if i not in S:
                        continue
                    s = S[i]
                    d = [f for f in C1F if int(s[f]) != int(r[f])]
                    if d:
                        kb.append((i, d, [s[f] for f in d],
                                   [r[f] for f in d]))
                cell['B_k1_tapes'] = len(K['rows'])
                cell['B_k1_mismatch'] = len(kb)
                cell['B_first'] = kb[0] if kb else None
            cell['C_repair_invariants'] = len(c_bad)
            cell['C_first'] = c_bad[0] if c_bad else None
            cell['C_no_period_invariant'] = len(inv_bad)
            cell['C_first_inv'] = inv_bad[0] if inv_bad else None
            cell['C_trim_over_n'] = trim_over_n
            out['%s|%d' % (rule, n)] = cell
            print('%-8s n=%-5d A b0bad=%d  B k1bad=%s  C repinv=%d noinv=%d '
                  'trim>n=%d' % (rule, n, cell['A_b0_mismatch'],
                                 cell.get('B_k1_mismatch', 'n/a'),
                                 cell['C_repair_invariants'],
                                 cell['C_no_period_invariant'],
                                 trim_over_n), flush=True)
    json.dump(out, open(os.path.join(HERE, 'col_check.json'), 'w'), indent=1)
    print('wrote col_check.json')
    return 0


if __name__ == '__main__':
    sys.exit(main())
