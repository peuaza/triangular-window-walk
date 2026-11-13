"""mp_tex -- the machine-profile table as LaTeX, from mp_cells.json.

For each rule, at the largest n with data: min_M G, min_M e, min_M g, the machine attaining each, and
whether one machine attains all three.  Values are in units of n.
"""
import json
import os
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
D = json.load(open(os.path.join(HERE, 'mp_cells.json')))
RULES = sorted({r['rule'] for r in D})

NAME = {'BAL>m/2': 'BAL_m/2', 'BAL>m/4': 'BAL_m/4', 'NBLK=1': 'NBLK_1',
        'NBLK>2': 'NBLK_gt2', 'NBLK>3': 'NBLK_gt3'}


def esc(s):
    s = NAME.get(s, s)
    return s.replace('_', '\\_')

rows = []
for rule in RULES:
    ns = sorted({r['n'] for r in D if r['rule'] == rule and r['reads']})
    if not ns:
        continue
    n = ns[-1]
    st = {}
    for r in D:
        if r['rule'] != rule or r['n'] != n or not r['reads']:
            continue
        c = np.asarray(r['reads'], float)
        st[r['machine']] = (float(np.exp(np.mean(np.log(c)))), float(np.mean(c)), float(c.max()))
    mG = min(st, key=lambda m: st[m][0]); me = min(st, key=lambda m: st[m][1]); mg = min(st, key=lambda m: st[m][2])
    same = len({mG, me, mg}) == 1
    rows.append((rule, n, st[mG][0] / n, st[me][1] / n, st[mg][2] / n,
                 mG, me, mg, same))

L = []
L.append("\\begin{table}[p]")
L.append("\\centering")
L.append("\\footnotesize")
L.append("\\begin{tabular}{@{}lrlrrrll@{}}")
L.append("\\toprule")
L.append("rule & $n$ & $\\min G/n$ & $\\min e/n$ & $\\min g/n$ & for $G$ & for $e$ & for $g$ \\\\")
L.append("\\midrule")
for rule, n, g, e, gg, mG, me, mg, same in rows:
    L.append("%s & %d & %.2f & %.2f & %.1f & %s & %s & %s \\\\" %
             (esc(rule).replace('\\\\_', '\\_'), n, g, e, gg, esc(mG), esc(me), esc(mg)))
L.append("\\bottomrule")
L.append("\\end{tabular}")
L.append("\\caption{The machine profile at the largest measured size: the three minima over the machines "
         "tested for each rule, in units of $n$, and the machine that attains each.  The minima are over "
         "the machines tested (\\texttt{pcr}, \\texttt{band\\_B64}, \\texttt{ces\\_kp\\_W64}, \\texttt{scan}), "
         "so each is an \\emph{upper} bound on the rule's own minimum, not the minimum over all machines.  "
         "Where the three columns name one machine, that machine attains all three minima at this size.}")
L.append("\\label{tab:profile}")
L.append("\\end{table}")

out = os.path.join(HERE, '_mp_profile.tex')
open(out, 'w', encoding='utf-8').write("\n".join(L) + "\n")
print("rows:", len(rows), "same-machine:", sum(1 for r in rows if r[8]), "->", out)
