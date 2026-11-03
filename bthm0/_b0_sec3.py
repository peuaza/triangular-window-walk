# -*- coding: utf-8 -*-
"""论文 4 的 §3 扩写：策略带、pos:0.3 的反常、slide 退化、n 依赖"""
import io

p = 'bthm0.tex'
lines = io.open(p, encoding='utf-8', newline='').read().split('\n')
B = chr(92)

i = [k for k, l in enumerate(lines) if l.startswith(B + 'section{How much a single pass explains}')][0]
j = [k for k, l in enumerate(lines) if l.startswith(B + 'section{Two scanners that predict perfectly}')][0]

new = [
    B + 'section{How much a single pass explains}',
    '',
    'Table~' + B + 'ref{tab:main} gives the headline value: the lag-one autocorrelation of the tape explains',
    '$R^2=.608$ of the variance of $' + B + 'delta$, with a bootstrap interval $.565$--$.645$ over $472$ meetings',
    'from $4000$ tapes. The rest of this section asks whether that number is a property of the machine',
    'or of the particular rule.',
    '',
    B + 'begin{table}[h]',
    B + 'centering', B + 'small',
    B + 'begin{tabular}{lcccc}',
    B + 'toprule',
    'rule & meetings & $' + B + 'sd(' + B + 'delta)$ & support & $R^2$ ' + B + B,
    B + 'midrule',
    'canonical four-way table & 150 & $.394$ & 91 & $.653$ ' + B + B,
    'jump of length two & 153 & $.394$ & 89 & $.633$ ' + B + B,
    'centre rule, $c=0.7$ & 151 & $.304$ & 84 & $.585$ ' + B + B,
    'centre rule, $c=0.6$ & 145 & $.281$ & 73 & $.509$ ' + B + B,
    'centre rule, $c=0.5$ & 161 & $.091$ & 47 & $.059$ ' + B + B,
    'position rule, $' + B + 'theta=0.7$ & 432 & $.132$ & 114 & $.064$ ' + B + B,
    'position rule, $' + B + 'theta=0.5$ & 684 & $.252$ & 209 & $.208$ ' + B + B,
    B + 'midrule',
    'position rule, $' + B + 'theta=0.3$ & $1200$ & $.350$ & 248 & $.579$ ' + B + B,
    'fixed slide & $1200$ & $0$ & 1 & --- ' + B + B,
    B + 'bottomrule', B + 'end{tabular}',
    B + 'caption{The canonical rule and seven variants at $n=256$, $1200$ tapes each. The two rows below',
    'the rule are the ones the constraint excludes or that have no dispersion: the position rule at',
    '$' + B + 'theta=0.3$ meets on every tape and therefore never cycles, and the fixed slide gives the same',
    'output on every tape.}',
    B + 'label{tab:band}', B + 'end{table}',
    '',
    'Two features stand out. First, the spread is wide: from $.059$ for the centre rule with a narrow',
    'band to $.653$ for the canonical table, a factor of eleven. So the headline number is not a',
    'constant of the machine; it is a property of the rule, and the rules differ in how much of their',
    'output the tape''s correlation tracks.',
    '',
    'Second, the rule that meets on every tape is the one that does not cycle, and its value $.579$ is',
    'below the canonical $.653$. That is the pattern the constraint in Section~5 formalises: a rule that',
    'always terminates is a rule whose trajectory never revisits a state, and such rules are the',
    'scanners, not the walks. The fixed slide at the other extreme produces one value for every tape',
    'and has no dispersion at all, so $R^2$ is undefined for it rather than one.',
    '',
    B + 'paragraph{Dependence on $n$.} Three sizes give $R^2=.655$, $.653$ and $.594$ at $n=128$, $256$',
    'and $512$, with $150$, $150$ and $124$ meetings. The meeting rate falls with $n$ ($.136$, $.125$,',
    '$.103$), which is expected: a longer tape offers more opportunities to close a cycle before the',
    'headers meet. The explained variance falls with it, by about $.06$ between $n=256$ and $n=512$;',
    'with samples of this size the drop is about one standard error, so we read it as suggestive',
    'rather than established.', '', '']

lines[i:j] = new
io.open(p, 'w', encoding='utf-8', newline='').write('\n'.join(lines))
print('section 3 expanded')
