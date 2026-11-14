"""Build bthm_ck.tex: the merged paper (the flag and the count) from bthm_c.tex.

Every substitution asserts its anchor, so a silent no-op is impossible.
"""
import io
import os
import re

HERE = os.path.dirname(os.path.abspath(__file__))
src = io.open(os.path.join(HERE, "bthm_c.tex"), encoding="utf-8").read()
T = src

TITLE_NEW = (r"\title{One query per state: a flag buys termination, a count buys behaviour}")

ABSTRACT_NEW = r"""\begin{abstract}
A tape machine with two headers, a window and three trims answers a query at each step; in the companion
paper its cost is measured on the tapes whose walk \emph{meets}. This paper asks what one query per state
buys the machine, and measures two readings of the same query over $18{,}000$ tape-runs in sixty (rule, $n$)
cells. The first reading is a \emph{flag}: has this window occurred before? It buys termination. Without it
$14.9\%$ of the runs meet, $61.9\%$ close a period and then repeat it for ever --- the run is stopped at its first repeat and recorded under
the cap category --- and
$23.2\%$ empty the region without ever meeting; the cycle is short, a median of two steps. With the flag the
meeting rate is $37.8\%$, $36.9\%$ of the runs that would die are rescued, and a run halts in at most
$(n+1)(n^{2}+1)$ steps because a repair consumes a cell. The flag is not a memory: the visited set is the
oracle's --- a run reached $2{,}945{,}071$ quadruples at $n=2048$ --- and the machine keeps $O(\log n)$ bits;
no lower bound on that memory is claimed. The
second reading is a \emph{count}: how often? The rule reads the count modulo $k$, and it buys behaviour on top of
termination: the meeting fraction rises by up to $19$ points as $k$ goes from $1$ to $4$, front-loaded at
$k=2$ for two of the three rules (for the third the largest step is $k=2\to3$), and between a third ($k=2$) and nearly half ($k=4$) of the rescues needs no repair at all, because the
rule's variant sends the walk the other way and the period never closes. Neither reading buys cheapness: every fitted exponent stays
on the rungs the catalogue already has, the read bound stays $O(n^{3})$ under the paper's own accounting (one charge per cell read), and the
distribution stays
heavy-tailed. Two failures are recorded for the same reason: a cheaper repair that buys meeting rate at an unbounded
cost, and one bit of writable ink per cell --- a mark read costs what a cell read costs, so a mark can only pay by answering a
question about several cells at once, and one bit cannot.
\end{abstract}"""

INTRO_NEITHER = r"""\section{What neither reading buys}\label{sec:neither}

Both readings buy termination, and the count buys behaviour on top of it; neither buys a cheaper order or a
lighter tail. The reason is structural rather than incidental. The charge is paid by the rule, one read per
action, so an oracle that answers a query without acting cannot make any rule read less: what it can do is
stop a walk from cycling, and stopping a walk does not change how many cells a rule must look at before it
decides. The failure recorded next is in the same spirit --- a repair, and a piece of ink, that buy nothing
on the count --- and it is reported because the negative is the same negative."""

OPEN_NEW = r"""\paragraph{The residue, beyond this design.} Section~\ref{sec:count} measures one design of the colour
variants, at $k\le4$, on three rules. What it leaves open is whether another design, another $k$, or a
variant that acts through the boundary cut --- these three tables carry no cut entry, so that half is
untested --- can move a cost onto a rung the ladder does not have, or lighten the tail. And the
colour-blind equivalence at $k=2$ says where the residue's value lies: entirely in the rule's reading of
it, not in the repair's."""

# 1. title
anchor = r"\title{One query bit stops the walk: cycle detection with a membership oracle}"
assert anchor in T, "title anchor"
T = T.replace(anchor, TITLE_NEW)

# 2. abstract
m = re.search(r"\\begin\{abstract\}.*?\\end\{abstract\}", T, re.S)
assert m, "abstract anchor"
T = T[:m.start()] + ABSTRACT_NEW + T[m.end():]

# 3. insert the count section, and open the "neither" section
mid = r"\section{Why an oracle, and not ink}\label{sec:not}"
assert mid in T, "ink section anchor"
sec5 = io.open(os.path.join(HERE, "_sec5_count.tex"), encoding="utf-8").read().strip()
T = T.replace(mid, sec5 + "\n\n" + INTRO_NEITHER + "\n\n" + r"\subsection{An oracle is not ink}\label{sec:not}")

# 4. the open paragraph about counting
old_open = re.search(r"\\paragraph\{Counting instead of remembering\.\}.*?different paper\.\n", T, re.S)
assert old_open, "open paragraph anchor"
T = T[:old_open.start()] + OPEN_NEW + "\n" + T[old_open.end():]

# 5. engines: name the colour engines
eng_anchor = r"that the two models are comparable; the reference engines charge differently"
assert eng_anchor in T, "engines anchor"
T = T.replace(
    eng_anchor,
    "that the two models are comparable; the colour model of Section~\\ref{sec:count} uses the same engine "
    "with the oracle answering a count instead of a flag, and its own independent engine agrees on $384$ "
    "sampled runs at $k=1\\dots4$; the reference engines charge differently",
)

# 6. the two forward references: the counting variant is Section 5 of the merged paper
for _old, _new in [
    ("the rule cannot recover a count from it (Section~\\ref{sec:not})",
     "the rule cannot recover a count from it (Section~\\ref{sec:count})"),
    ("The counting variant is a different model, and Section~\\ref{sec:not} says why.",
     "The counting variant is a different model; Section~\\ref{sec:count} measures it."),
]:
    assert _old in T, "forward-ref anchor missing"
    T = T.replace(_old, _new)

io.open(os.path.join(HERE, "bthm_ck.tex"), "w", encoding="utf-8", newline="\n").write(T)
print("wrote bthm_ck.tex  %d lines" % T.count("\n"))
# quick structural report
for pat in [r"\title\{", r"\section\{", r"\subsection\{"]:
    print(pat, len(re.findall(pat, T)))
