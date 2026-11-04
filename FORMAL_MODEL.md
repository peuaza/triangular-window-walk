# A formal model for the pointer/walk classification

*Drop-in section. Notation follows the manuscript: $L,R$ are the region boundaries, $\ell,r$ the
headers, $\delta=\ell/(n-1)$ the output, $\varphi$ the window function, $f$ the action map, and
$\mathcal R$ the set of rules. Definition and theorem numbers are local to this file and should be
renumbered on insertion; §7 lists where each piece replaces or supplements the current text.*

---

## 1. The ensemble

Throughout $n\ge 3$ and $[n]=\{0,\dots,n-1\}$.

\begin{definition}[Tape]\label{def:tape}
Let $\Gamma=\{0,\tfrac12,1\}$. A \emph{tape of length $n$} is a word $\tau\in\Gamma^n$; $\tau_i$ is the
value of cell $i$. The tape is the whole input: no machine below is given anything else.
\end{definition}

\begin{definition}[Ensemble]\label{def:ensemble}
Put $h=h(n)=\lceil n/2\rceil$. Let $X,Y$ be independent random subsets of $[n]$, each uniform over the
$h$-element subsets, and set
\[
\tau_i \;=\; \tfrac12\bigl(\mathbf 1[i\in X]+\mathbf 1[i\in Y]\bigr),\qquad 0\le i\le n-1 .
\]
Write $\mu_n$ for the law of $\tau$ on $\Gamma^n$ and $\mathcal T_n=\operatorname{supp}\mu_n$. Tapes at
different $n$ are independent; no claim below uses more than the marginal $\mu_n$.
\end{definition}

\begin{lemma}[Cell law]\label{lem:cell}
Under $\mu_n$, for even $n$
\[
P(\tau_i=0)=P(\tau_i=1)=\tfrac14,\qquad P(\tau_i=\tfrac12)=\tfrac12,
\]
and for odd $n$, with $h=(n+1)/2$,
\[
P(\tau_i=0)=\Bigl(\frac{n-1}{2n}\Bigr)^{2},\qquad
P(\tau_i=\tfrac12)=\frac{n^2-1}{2n^2},\qquad
P(\tau_i=1)=\Bigl(\frac{n+1}{2n}\Bigr)^{2}.
\]
For $i\ne j$, $\operatorname{Cov}(\mathbf 1[\tau_i=0],\mathbf 1[\tau_j=0])<0$.
\end{lemma}

\begin{proof}
$\tau_i=0$ iff both columns miss $i$, of probability $((n-h)/n)^2$; $\tau_i=1$ iff both hit it, of
probability $(h/n)^2$; $\tau_i=\tfrac12$ is the complement. The negative covariance is sampling without
replacement.
\end{proof}

\begin{lemma}[Support]\label{lem:support}
With $\#_v\tau=|\{i:\tau_i=v\}|$,
\[
\mathcal T_n \;=\; \bigl\{\tau\in\Gamma^n \;:\; \#_1\tau-\#_0\tau=2h-n
\ \text{ and }\ \#_{\frac12}\tau \text{ is even}\bigr\}.
\]
Hence on every tape of the ensemble $\#_1\tau=\#_0\tau$ for even $n$, and $\#_1\tau=\#_0\tau+1$ for odd
$n$. The cardinality is $\sum_a\binom na\binom{n-a}{2h-2a}$, so $|\mathcal T_6|=141$ and
$|\mathcal T_{13}|=201\,643$.
\end{lemma}

\begin{proof}
A tape determines $A=\{i:\tau_i=1\}=X\cap Y$ and $B=\{i:\tau_i\ne\tfrac12\}=X\cup Y$, whence
$|A|+|B|=2h$. Conversely, given $A\subseteq B$ with $|A|+|B|=2h$, choose $S\subseteq B\setminus A$ with
$|S|=h-|A|$ and put $X=A\cup S$, $Y=A\cup(B\setminus S)$; the resulting tape is $\tau$. The map is a
bijection, and $\#_1=|A|$, $\#_0=n-|B|$, $\#_{\frac12}=|B|-|A|$.
\end{proof}

\begin{remark}[Three consequences]\label{rem:support}
(i) $\mathcal T_n$ is thin: $141$ of the $3^6=729$ words at $n=6$. (ii) Every worst-case quantifier in
this paper ranges over $\mathcal T_n$ and not over $\Gamma^n$; a rule may be arbitrarily bad off the
support and still lie in every class defined below. (iii) The manuscript's phrase "every admissible
sequence" already denotes $\mathcal T_n$ --- the counts $141$ and $201\,643$ are exactly the
cardinalities in Lemma~\ref{lem:support} --- but the set itself is never defined there.
\end{remark}

---

## 2. The device

The manuscript's only definition of a read is one sentence ("A read is one access to one cell"), and it
has no machine at all. This subsection supplies both. The device is deliberately weak --- a finite-state
transducer with heads, not a register machine --- so that no arithmetic on cell values can be smuggled
into a step.

\begin{definition}[Device]\label{def:device}
A \emph{device} is a tuple $M=(Q,q_0,q_h,k,s,I)$ where
\begin{itemize}
\item $Q$ is a finite set of states, with start state $q_0$ and halting state $q_h$;
\item $k\in\{1,2\}$ is the number of input heads;
\item $s:\mathbb N\to\mathbb N$ bounds the work tape, which is binary and carries one head;
\item $I$ is a finite instruction set with instructions of exactly two kinds:
\emph{sense}, $\mathrm{sense}(q,a_1,\dots,a_k,b)=(q',\allowbreak b',e_1,\dots,e_k,e_w)$, which observes the symbols
$a_j\in\Gamma$ under the $k$ input heads and the bit $b$ under the work head, then displaces every head
by $e_j,e_w\in\{-1,0,+1\}$ and writes $b'$; and \emph{move},
$\mathrm{move}(q,b)=(q',b',e_1,\dots,e_k,e_w)$, which is the same but does not consult the input symbols.
\end{itemize}
The input tape is read-only and holds $\tau$. The output tape is write-once, left to right, binary, and
carries the answer; it is not charged to $s$.
\end{definition}

\begin{definition}[Step and read]\label{def:read}
A \emph{step} is one execution of one instruction of $I$. A \emph{read} is one cell examination: a
\emph{sense} instruction executed with input head $j$ at cell $i$ is one read of cell $i$,
\emph{counted with multiplicity} --- sensing the same cell twice is two reads. A \emph{move} instruction
costs one step and \emph{no} read. A pass that senses each of $m$ consecutive cells therefore costs $m$
reads, and the $m-1$ displacements between them are free in reads.
\end{definition}

\begin{definition}[Cost profile]\label{def:cost}
For a device $M$ and $\tau\in\Gamma^n$, $\mathrm{reads}_M(\tau)$ is the number of \emph{sense}
instructions executed before halting, $\mathrm{steps}_M(\tau)$ the number of instructions executed, and
$\mathrm{space}_M(\tau)$ the number of work-tape cells visited. $M$ is \emph{$s(n)$-bounded} if
$\mathrm{space}_M(\tau)\le s(n)$ for every $\tau\in\Gamma^n$.
\end{definition}

\begin{remark}[What is not charged, and why it is harmless]\label{rem:uncharged}
Head positions, the input length, the output tape and the number of heads are not charged; only the work
tape is charged to $s(n)$, which is the standard convention and the one "$O(\log n)$ bits of state"
intends. Free head positions create no loophole, because a finite control cannot read a head position as
a number --- it can only move the head --- so a device that must remember a position has to spend
work-tape bits on it. This is exactly the difference between a sensing-head automaton, which can detect
that two heads coincide, and one that can compute with their coordinates.
\end{remark}

\begin{remark}[Indistinguishability]\label{rem:indist}
If two tapes agree on every cell that $M$ senses during a run, the two runs are identical, step for step.
Hence if a set $S$ of cells determines $\delta_{\rho,n}(\tau)$ for every $\tau\in\mathcal T_n$, every
device computing $\delta_{\rho,n}$ reads at least $|S|$ cells on some tape --- and since $|S|\le n$, this
form of argument can never prove more than a linear bound. This is the exchange argument of §6 of the
manuscript, stated once and in general.
\end{remark}

\begin{remark}[No arithmetic on cell values]\label{rem:arithmetic}
Instructions branch on symbols; they cannot add, compare or average cell values. Every arithmetic
operation the simulators perform --- the integer comparisons $Sk>S_hm$ and $(m-1)R_r-AB>0$, the updates
of the window sums --- must be realised as a finite program on the $s$-bit work tape, at one step per
work-tape operation, and costs \emph{no reads}. This reproduces the manuscript's accounting exactly:
reads are charged for cells sensed and never for register operations, which is why the register count is
invisible in Table~5.1 but visible in $\mathrm{steps}$.
\end{remark}

\begin{remark}[Head count]\label{rem:heads}
$k=1$ and $k=2$ give different classes (Remark~\ref{rem:relations}). The manuscript's machine has two
heads and its only lower bound is for one; the classes below carry the head count as a subscript so that
the two cannot be confused.
\end{remark}

---

## 3. The rule, and the two readings of $\varphi$

\begin{definition}[Rule]\label{def:rule}
A \emph{rule} is a pair $\rho=(\varphi,f)$ where
\[
\varphi:W\to A,\qquad W=\bigcup_{m\ge1}\Gamma^m,\qquad A \text{ finite},
\]
\[
f:A\times\mathbb N\times\mathbb N^4\longrightarrow
\{\ell{+}1,\ \ell{-}1,\ r{+}1,\ r{-}1,\ \text{left cut},\ \text{right cut},\ \text{double cut},\
\text{halt}\},
\]
evaluated at $\bigl(\varphi(\tau[\ell..r]),\,n,\,\ell,\,r,\,L,\,R\bigr)$. Two conditions hold by
definition and are worth naming, because the manuscript uses them without stating them: $\varphi$ depends
on the current window only --- not on the history, so $(\varphi,f)$ makes the state $(L,R,\ell,r)$
Markov --- and all conventions by which a rule resolves an inadmissible move, a cut or a halt are part of
$f$. No computability, complexity or measurability condition is placed on $\varphi$: it is a mathematical
function and it is the \emph{specification} of the rule, not an object the rule computes.
\end{definition}

\begin{definition}[Run]\label{def:run}
Fix $n$ and $\tau\in\Gamma^n$. The \emph{state} is $(L,R,\ell,r)$ with
$0\le L\le\ell\le r\le R\le n-1$, and the initial state is $(0,n-1,0,n-1)$. Given a state, let
$a=f(\varphi(\tau[\ell..r]),n,\ell,r,L,R)$.
\begin{enumerate}
\item If $\ell=r$, the run halts and outputs $\delta=\ell/(n-1)$.
\item If $a$ is a move taking a header to $p$, it is \emph{admissible} iff $L\le p\le R$. If admissible
and $p$ meets or passes the other header, the run halts and outputs $\delta=p/(n-1)$. If inadmissible,
the run aborts with no output.
\item If $a$ is a cut, it is admissible iff the new region is non-empty; then $(L,R)\leftarrow(L',R')$ and
$(\ell,r)\leftarrow(L',R')$.
\item If $a$ is \emph{halt}, the run halts with no output.
\item If the new state has occurred earlier in the run, the run halts with no output.
\end{enumerate}
\end{definition}

\begin{definition}[Output function]\label{def:output}
$\delta_{\rho,n}:\Gamma^n\to D_n\cup\{\bot\}$, where $D_n=\{i/(n-1):0\le i\le n-1\}$ and $\bot$ denotes
"no output", is the map defined by Definition~\ref{def:run}. It is total: the state space is finite, so
every run halts by clause~5.
\end{definition}

\begin{remark}[The repeat test costs no space]\label{rem:repeat}
Clause~5 is a property of the orbit and not a resource of the device: $f$ and the geometry determine the
next state, so a device can implement the test by Brent's algorithm in $O(1)$ extra work bits. It costs
reads and steps, never space. The clause therefore does not smuggle an $\Omega(n^2)$-bit table into the
model --- but a device that merely simulates must pay for cycle detection, and the manuscript's read
counts already include that payment.
\end{remark}

\begin{definition}[Oracle semantics, (O)]\label{def:oracle}
$\varphi$ is an oracle: evaluating $\varphi(w)$ costs $0$ reads and $0$ steps, for every window $w$.
\end{definition}

\begin{proposition}[The oracle reading is vacuous]\label{prop:collapse}
Under (O), every rule $\rho$ satisfies $\rho\in\mathsf{Pointer}_2(0,\allowbreak O(\log n))$, and
$\mathsf{Walk}_2=\varnothing$.
\end{proposition}

\begin{proof}
The device stores $(L,R,\ell,r)$ --- four integers below $n$, so $4\lceil\log_2 n\rceil$ bits --- queries
the oracle for $\varphi(\tau[\ell..r])$, looks $f$ up in its finite control, and updates. It never senses
a cell, so $\mathrm{reads}=0$. Clauses 1--5 of Definition~\ref{def:run} are decided from the stored state
and the oracle's answer.
\end{proof}

\noindent\emph{Consequence.} The manuscript cannot be using (O) for its classification, on pain of making
it trivial. Every read count in §3 and §5 of the manuscript is therefore a statement under the next
definition.

\begin{definition}[Computed semantics, (C)]\label{def:computed}
The device has no oracle. On input $\tau$ it must write $\delta_{\rho,n}(\tau)$ on its output tape,
sensing only cells of $\tau$. $\varphi$ enters only through $\delta_{\rho,n}$; the device may compute
$\varphi$ however it likes, precompute it, cache it, or bypass it entirely.
\end{definition}

\begin{definition}[$\varphi$-state; incremental computability]\label{def:phistate}
A \emph{$\varphi$-state} for $\rho$ on a window $w$ is a value $\sigma$ in a set $\Sigma$ together with an
evaluation map $\mathrm{ev}:\Sigma\times\mathbb N\to A$ such that $\mathrm{ev}(\sigma,|w|)=\varphi(w)$. A
device \emph{carries} $\varphi$ if its work tape holds a $\varphi$-state and $\mathrm{ev}$ costs $O(1)$
steps and $0$ reads. $\varphi$ is \emph{incrementally computable} if there is
$\mathrm{upd}:\Sigma\times\Gamma\to\Sigma$, costing $O(1)$ steps and $0$ reads, with
$\mathrm{upd}(\sigma_w,c)=\sigma_{w\setminus c}$ for every window $w$ and every cell $c$ at either end of
$w$.
\end{definition}

\begin{remark}[Which reading the manuscript uses where]\label{rem:which}
All three readings are present, and the manuscript never says which is in force.
\begin{center}
\small
\begin{tabular}{@{}p{.27\textwidth}p{.20\textwidth}p{.44\textwidth}@{}}
\toprule
Location & Reading & What it actually establishes\\
\midrule
§2, Definition of strategy (l.122--129) & (O) & the \emph{semantics} of the rule: $f$ is handed
$\varphi(\tau[\ell..r])$. Not a cost statement.\\
§3.1 read counts, §5 simulators & (C) & read counts of devices that compute $\varphi$ themselves.\\
Theorem 3.1 (\texttt{thm:monotone}) & (C) with a carried $\varphi$-state & the only place the assumption
is written down; strictly stronger than the Definition of strategy.\\
§7 ceiling and its predictors & neither & a statistical estimator with a read budget, not required to be
exact. Not a device in the sense of §2.\\
\bottomrule
\end{tabular}
\end{center}
\end{remark}

---

## 4. The classes

\begin{definition}[Pointer]\label{def:pointer}
Let $T,s:\mathbb N\to\mathbb N$ and $k\in\{1,2\}$. A rule $\rho$ belongs to
$\mathsf{Pointer}_k(T,s)$ iff there is a \emph{single} device $M$ with $k$ heads and $s(n)$ work-tape bits
such that
\[
\forall n\ge3,\ \forall\tau\in\mathcal T_n:\quad
\mathrm{reads}_M(\tau)\le T(n),\quad
\mathrm{space}_M(\tau)\le s(n),\quad
M(1^n,\tau)=\delta_{\rho,n}(\tau).
\]
The device is uniform in $n$ --- one program, handed $1^n$. The quantifier over tapes is universal and
worst case: one tape on which the device errs or exceeds $T(n)$ puts $\rho$ outside the class. $\bot$ is
an output like any other, so a device must also certify that the run halts or cycles.
\end{definition}

\begin{definition}[Walk]\label{def:walk}
$\mathsf{Walk}_k(T,s)=\mathcal R\setminus\mathsf{Pointer}_k(T,s)$, where $\mathcal R$ is the set of rules
of Definition~\ref{def:rule}.
\end{definition}

\begin{definition}[The named classes]\label{def:named}
With $O(n)$ reads and $O(\log n)$ bits,
\[
\mathsf{Pointer}_k=\bigcup_{c\ge1}\mathsf{Pointer}_k(c\,n,\,c\log n),
\qquad
\mathsf{Walk}_k=\mathcal R\setminus\mathsf{Pointer}_k .
\]
Unwinding: $\rho\in\mathsf{Walk}_k$ iff for every $c$ and every $k$-head device $M$ with $c\log n$ bits
there are infinitely many $n$ and a tape $\tau\in\mathcal T_n$ on which $M$ errs or reads more than
$c\,n$ cells. Membership in $\mathsf{Walk}_k$ is a statement about \emph{all} devices, and it is what a
lower bound has to prove.
\end{definition}

\begin{definition}[Average-case and high-probability variants]\label{def:avg}
The same clauses with the tape quantifier weakened:
\[
\mathsf{Pointer}^{\mathrm{avg}}_k(T,s):\quad \exists M\ \forall n:\ \
\mathbb E_{\tau\sim\mu_n}\bigl[\mathrm{reads}_M(\tau)\bigr]\le T(n)
\ \text{ and }\ M \text{ correct on every } \tau\in\mathcal T_n;
\]
\[
\mathsf{Pointer}^{1-\varepsilon}_k(T,s):\quad \exists M\ \forall n:\ \
P\bigl(M \text{ correct and } \mathrm{reads}_M\le T(n)\bigr)\ge 1-\varepsilon(n);
\]
\[
\mathsf{Pointer}^{\mathrm{avg}|\mathrm{out}}_k(T,s):\quad \text{the first, with } \mu_n
\text{ conditioned on } \delta_{\rho,n}(\tau)\ne\bot .
\]
\end{definition}

\begin{remark}[Elementary relations]\label{rem:relations}
\[
\mathsf{Pointer}_k(T,s)\subseteq\mathsf{Pointer}^{\mathrm{avg}}_k(T,s)
\subseteq\mathsf{Pointer}^{\mathrm{avg}|\mathrm{out}}_k(T,s),
\]
\[
\mathsf{Pointer}_1\subseteq\mathsf{Pointer}_2,
\qquad
\mathsf{Pointer}_k(T,s)\subseteq\mathsf{Pointer}_k(T',s')\ \ \text{if } T\le T',\ s\le s',
\]
\[
\mathsf{Pointer}^{1-\varepsilon}_k(T,s)\subseteq
\mathsf{Pointer}^{\mathrm{avg}}_k\bigl(T+\varepsilon R_{\max},\,s\bigr),
\]
where $R_{\max}=\max_{\tau}\mathrm{reads}_M(\tau)$. The manuscript's evidence does not sit in a single
one of these: its upper bounds are worst case, its tables are
$\mathsf{Pointer}^{\mathrm{avg}|\mathrm{out}}$, and its one proved lower bound is an expectation.
\end{remark}

\begin{remark}[The space axis is pinned at $\Theta(\log n)$]\label{rem:space}
Every rule in this paper is computable in polynomial time and $O(\log n)$ bits: the simulators of §5 run
in $O(\log n)$ bits. Hence $\omega(\log n)$ space is not merely unproved but false
for every rule in the family, and the classification cannot be a separation in space. In the other
direction, a device that \emph{simulates} the run must hold two of the four coordinates $L,R,\ell,r$ in
work tape --- a two-head device can park at most two of them under heads --- so $s(n)=\Omega(\log n)$ for
every rule whose $f$ depends on all four, which includes the canonical rule. The space axis is therefore
closed at both ends at $\Theta(\log n)$, and $T$ is the only free coordinate. This is the formal content
of the manuscript's remark at l.901--903, and it belongs before the conjecture, because it is what makes
the conjecture a statement about repeat reads.
\end{remark}

\begin{remark}[Output and naming conventions]\label{rem:conv}
$\delta\in D_n\cup\{\bot\}$ is written in binary on the write-once output tape. Every claim here has
$s(n)\ge\log_2 n$, so the device can name the answer; the convention matters only at $s(n)=o(\log n)$,
where the class degenerates, and no claim is of that kind. "Pointer" does not mean a pointer machine in
the sense of Tarjan; the word is the manuscript's, and the collision should be disclaimed in a footnote.
\end{remark}

---

## 5. What the manuscript's results are, in this language

Throughout $\rho_{\mathrm{can}}$ is the canonical rule of §2, $\rho_0$ the rule of §6 with
$b_1\equiv0$, and $\rho_1,\dots,\rho_4$ the four rules of §3.2.

\begin{theorem}[Monotone rules: proved]\label{thm:A}
Let $\rho=(\varphi,f)$ satisfy: $f$ never returns a cut or \emph{halt} and never enlarges the window;
$\varphi$ is incrementally computable in the sense of Definition~\ref{def:phistate}; and $\rho$ admits an
initial $\varphi$-state computable from the whole window in $n$ reads with $O(\log n)$ bits. Then
\[
\rho\in\mathsf{Pointer}_1(2n,\,O(\log n)).
\]
\end{theorem}

This is \texttt{thm:monotone}, and it is the only unconditional \emph{upper} membership the paper proves.
The third hypothesis is not implied by Definition~\ref{def:rule}: a finite value set for $\varphi$
supplies no bound whatever on the cost of evaluating it, since $\varphi$ restricted to length $m$ is an
arbitrary table of $3^m$ entries. The theorem therefore establishes membership for a proper subclass of
$\mathcal R$.

\begin{theorem}[Canonical rule: measured, not proved]\label{thm:B}
For $32\le n\le1024$ and every tape tested, $\rho_{\mathrm{can}}$ is computed by a two-head device with
$O(\log n)$ bits in at most $5.4\,n$ reads, and $\mathrm{reads}/n$ is flat over that range. In the
language of Definition~\ref{def:named} the manuscript therefore has
\[
\rho_{\mathrm{can}}\in\mathsf{Pointer}_2(5.4\,n,\,O(\log n))\qquad\text{on the tested range only,}
\]
not for all $n$. The gap is exactly what §3.1 concedes: the device costs $O(1)$ reads per step and the
\emph{number of steps is measured}. Membership is a conjecture supported by a flat ratio over a
thirtyfold range in $n$, and the paper should either prove a step bound or say so inside the claim.
\end{theorem}

\begin{theorem}[The quadratic cost: proved for the trajectory, open for the reads]\label{thm:C}
For $\rho_0$, with $s_\rho(\tau)$ the number of steps of the run of $\rho$ on $\tau$,
\begin{gather*}
\text{(i)}\quad s_{\rho_0}(\tau)\le n(n+1)/2 \quad \text{for every } \tau\in\Gamma^n;\\
\text{(ii)}\quad \mathbb E_{\mu_n}[s_{\rho_0}]\ge 8.1\cdot10^{-6}\,n^{2} \quad \text{for all large } n;
\end{gather*}
and (iii) if $M$ computes $\delta_{\rho_0}$ by following the run and its memory is additive, then $M$
reads $\Omega(n^{2})$ cells in expectation.
\end{theorem}

Statements (i)--(ii) are the corollary and the epoch bound of §6. Statement (iii) is the accounting of
l.845--848: additive memory stores the values of distinct cells in disjoint fields, so on an epoch of
length $\ell_e$ a device with $s$ bits reads at least $\ell_e-O(s)$ cells, and
$\sum_e\max(0,\ell_e-2s)$ at $s=O(\log n)$ is again $\Omega(n^2)$. All three bound the \emph{trajectory}
and the reads of devices that follow it. None bounds the reads of an arbitrary device: a device may
compute $\delta_{\rho_0}$ without simulating anything, and no argument in the paper excludes that. This
is the precise location of the gap between "$\Theta(n^2)$ cost" and "$\rho_0\in\mathsf{Walk}$".

\begin{theorem}[Single head: conditional]\label{thm:D}
Assume
\[
(H_\Sigma)\qquad \sum_{k=1}^{n-1}\mathrm{CC}_k(\delta_{\rho,n})\ \ge\ c\,n^{2}
\quad\text{for all large } n,\ c>0,
\]
where $\mathrm{CC}_k$ is the two-way communication complexity of $\delta_{\rho,n}$ at the cut after cell
$k$. Then every one-head device with $s$ bits that computes $\delta_\rho$ satisfies
\[
T\,\bar J\,s\ \ge\ \sum_k \mathrm{CC}_k(\delta_\rho)-n s\ \ge\ c\,n^{2}-n s ,
\]
with $\bar J$ the mean jump length of its read sequence. Hence, whenever $\bar J=O(1)$ --- in particular
for a device whose head advances one cell per read --- $T=\Omega(n^{2}/s)$, and at $s=O(\log n)$,
$T=\Omega(n^{2}/\log n)$. For random access the bound degrades by the factor $\bar J$ and proves nothing,
which is the manuscript's own caveat at l.897--899. In particular
$\rho\notin\mathsf{Pointer}_1(o(n^{2}/\log n),O(\log n))$ --- \emph{for one head only}. The
classification of §4 is about the two-head machine, so this places no rule in $\mathsf{Walk}_2$.
\end{theorem}

$(H_\Sigma)$ is not proved. What is proved is the \emph{upper} bound
$\sum_k\mathrm{CC}_k\le\log_23\cdot n^{2}/4=.396\,n^{2}$, which is the wrong direction; the lower bound
comes from the measurement that $\sum_k\mathrm{CC}_k$ attains $95$--$98$ per cent of that maximum over all
balanced cuts, giving $c\approx.37$ and $T\gtrsim.38\,n^{2}/s$. The theorem is conditional on a
measurement and should be labelled so. The inequality itself is asserted at l.885 without proof.

\begin{theorem}[$p$-pass: conditional]\label{thm:E}
Call a device \emph{$p$-pass} if the cells it senses, in order, form $p$ monotone sweeps of $[n]$, each
sensing every cell, so that $T=pn$ and $\bar J=1$; a device that sweeps only part of the tape is not
$p$-pass and the bound does not apply to it. If
\[
(H_{\mathrm{rank}})\qquad\text{every balanced cut } k \ (\min(k,n-k)\ge n/4) \text{ has communication
matrix of rank } \ge c' n,
\]
--- the manuscript's rank conjecture, verified exhaustively at $n=8,10$ --- then each cut is crossed at
most $2p$ times, each crossing transfers at most $s$ bits, so $2ps\ge c'n$ and
\[
p=\Omega(n/s),\qquad T=pn=\Omega(n^{2}/s)=\Omega(n^{2}/\log n)\ \text{ at } s=O(\log n).
\]
\end{theorem}

Theorem~\ref{thm:E} is Theorem~\ref{thm:D} at $\bar J=1$ and $T=pn$. The two bounds are one bound at two
values of $\bar J$, and presenting them that way removes the impression that the paper has two independent
lower bounds. Together they are the whole of the linear-time exclusion for one head.

\begin{theorem}[No space lower bound can exist: proved]\label{thm:F}
For every rule $\rho$ in the paper's family, $\delta_\rho\in\mathsf{DTISP}(\mathrm{poly}\,n,O(\log n))$.
Hence no rule is excluded from $\mathsf{Pointer}$ by a space bound, the map
$s\mapsto\mathsf{Pointer}_k(T,s)$ saturates at $s=\Theta(\log n)$ for the rules considered
(Remark~\ref{rem:space}), and a lower bound of the form "$\delta_\rho$ needs $\omega(\log n)$ space" is
false for every rule in the family. The only possible separation is in $T$ at $s=\Theta(\log n)$.
\end{theorem}

\begin{theorem}[General upper bound: restated]\label{thm:G}
Every rule whose $\varphi$ is computable by a device with $O(n^{2})$ bits and $O(n^{3})$ steps is
computable with $O(n^{2})$ bits and $O(n^{3})$ evaluations of $\varphi$; if a window of length $m$ costs
$O(m)$ reads to evaluate, the total is $O(n^{4})$ reads.
\end{theorem}

The manuscript's \texttt{thm:upper} states this with "a step is one evaluation of $\varphi$ at unit cost",
which is the oracle reading (O) and is false under (C) --- under (C) $\varphi$ may be uncomputable, and
the theorem's subject class is empty of content. The restatement is the correct one.

\paragraph{The ceiling, formalised.} A \emph{predictor with budget $B$} is a map from a finite sample of
tapes to a function $\hat g:\Gamma^n\to[0,1]$ whose evaluation on any $\tau$ senses at most $B(n)$ cells;
$R^{2}$ is computed out of sample, conditional on $\delta\ne\bot$. The conjecture of §7 is then
\[
\rho\in\mathsf{Walk}_2\ \Longrightarrow\
\sup\bigl\{\mathbb E[R^{2}(\hat g)]:\ \hat g \text{ a predictor with budget } o(n^{2}) \text{ and }
O(\log n) \text{ bits}\bigr\}\ \le\ \kappa<1 .
\]
Two things must be said. First, the hypothesis is worst case over tapes (Definition~\ref{def:walk}) while
the conclusion is an average over $\mu_n$; nothing in the paper relates the two, and the conjecture is not
well posed until one of them is changed. The natural repair is to state the hypothesis with
$\mathsf{Walk}^{\mathrm{avg}}_2$ (Definition~\ref{def:avg}), or to assert that the two coincide for
memoryless rules --- which the paper neither states nor argues. Second, "predictor" is not defined in the
manuscript at all, so "every predictor" in the current text quantifies over an unlisted set; the
definition above fixes the set.

\paragraph{Summary.} The table collects the statements above. The Status column is the part that matters
for the objection: one row is an unconditional membership, two are measured, two are conditional (on a
measurement and on a conjecture respectively), one is not established at all, and three are theorems
about a different resource than the one the classification is named after.

\begin{center}
\small
\begin{tabular}{@{}lll@{}}
\toprule
Claim & Formal statement & Status\\
\midrule
monotone rules are cheap & $\rho\in\mathsf{Pointer}_1(2n,O(\log n))$ & \textbf{proved}
(Thm~\ref{thm:A})\\
canonical rule is cheap & $\rho_{\mathrm{can}}\in\mathsf{Pointer}_2(5.4n,O(\log n))$, $n\le1024$ &
measured (Thm~\ref{thm:B})\\
interior cut, ignore-window & $\in\mathsf{Pointer}_2(cn,O(\log n))$, $n\le512$ & measured\\
$\rho_0$ costs $\Theta(n^{2})$ & $s_{\rho_0}=\Theta(n^{2})$; reads $\Omega(n^{2})$ if following &
\textbf{proved} on the trajectory (Thm~\ref{thm:C})\\
$\rho_1,\dots,\rho_4$ are expensive & $\rho_i\in\mathsf{Walk}_2$ & \textbf{not established}\\
single head needs $\Omega(n^{2}/\log n)$ & $\rho\notin\mathsf{Pointer}_1(o(n^{2}/\log n),O(\log n))$ &
conditional on $(H_\Sigma)$ (Thm~\ref{thm:D})\\
$p$ passes need $\Omega(n^{2}/\log n)$ & $p=\Omega(n/s)$, $T=\Omega(n^{2}/s)$ & conditional on
$(H_{\mathrm{rank}})$ (Thm~\ref{thm:E})\\
no space lower bound & $\delta_\rho\in\mathsf{DTISP}(\mathrm{poly},O(\log n))$ & \textbf{proved}
(Thm~\ref{thm:F})\\
every rule is simulable & $O(n^{3})$ evaluations of $\varphi$ & \textbf{proved} under (O) only
(Thm~\ref{thm:G})\\
\bottomrule
\end{tabular}
\end{center}

---

## 6. Checklist: where the formal model and the current text disagree

Each item gives the location in \texttt{bthm\_merged.tex}, the deviation, and a one-line fix. Items
marked **[load-bearing]** change what the paper proves, not merely how it says it.

1. **[load-bearing] §3, Definition (l.186--191).** "A strategy is a pointer strategy if its output can be
   computed in $O(n)$ cell reads with $O(\log n)$ bits of state, and a walk strategy if it cannot" has no
   quantifier over $n$, none over tapes, no uniformity condition, no head count, no output convention,
   and no definition of "state". *Fix:* replace with Definitions~\ref{def:pointer}--\ref{def:avg}.

2. **[load-bearing] §2, Definition (l.122--129).** $\varphi$ is "a fixed function of the window with a
   finite value set" and nothing more. A finite value set gives neither computability nor a space
   bound --- the table on length $m$ has $3^m$ entries --- and under the reading that lets $f$ see
   $\varphi$ for free the whole classification collapses (Proposition~\ref{prop:collapse}). *Fix:* insert
   Definition~\ref{def:rule} and one sentence saying that $\varphi$ is a specification and its evaluation
   is charged to the device; cite Definition~\ref{def:phistate} wherever it is assumed incremental.

3. **[load-bearing] §2 "Reads" (l.104--105).** "A read is one access to one cell of the tape" does not say
   whether a re-read is charged again, whether moving a head between cells costs anything, or whether the
   cells a head passes over are charged. *Fix:* replace with Definitions~\ref{def:read}--\ref{def:cost}
   (multiplicity counted; movement free in reads, one step).

4. **[load-bearing] §2, run definition (l.131--145).** No clause says what happens when $f$ returns a move
   that leaves $[L,R]$. Line 147 says the model "does not repair it" but not what it does instead, and
   the two readings --- abort with $\bot$, or leave the state unchanged --- give different $\delta$.
   *Fix:* add clause 2 of Definition~\ref{def:run}, or require $f$ to be total on admissible actions,
   which the paper's rules are by their conventions.

5. **[load-bearing] §3.2 heading and prose (l.285--291), and "The class has four candidates"
   (l.581--583).** The heading says four rules \emph{are in} the walk class; the text says "no $O(n)$
   simulator is known", and §3.1 leaves two further rules unclassified. Under Definition~\ref{def:walk}
   every rule is in exactly one class, so "not classified" can only mean "we do not know which", and
   **no rule in the paper is proved to lie in $\mathsf{Walk}_2$** --- not the four, not any. *Fix:*
   retitle "Candidates for the walk class" and add one sentence: membership in $\mathsf{Walk}$ is
   established for no rule here; what is established is the absence of a known $O(n)$ simulator.

6. **[load-bearing] §2 (l.147--152) against §2, canonical rule (l.157--161).** Line 151 says the cut and
   halt conventions are "choices of those rules, not clauses of the model", yet the canonical rule's
   definition presents them as part of the rule --- and they change $\delta$, since the halt convention
   decides between a value and $\bot$. *Fix:* keep both conventions inside $f$ and say once, in
   Definition~\ref{def:rule}, that every convention is part of $f$.

7. **[load-bearing] §2 (l.107--116) against §3.1 (l.234) and (l.252).** The model has eight \emph{unit}
   moves; §3.1 measures rules "when a move moves the header by two or three cells", and the far-end rule
   uses "a cut at an interior position $p$ \dots\ an action that the machine of §2 does not have". Both
   are outside the stated model, and the paper says so only in passing. *Fix:* define the action set with
   a step size $d$ and an interior cut, and state that the base model is $d=1$, $p\in\{L,R\}$.

8. **[load-bearing] §2 (l.87--92) and §6 (l.679).** The cell law is given but the support is not, and the
   exhaustiveness claim is over "every admissible sequence" --- which is $\mathcal T_n$ of
   Lemma~\ref{lem:support}, a thin set ($141$ of $729$ words at $n=6$). Every worst-case quantifier must
   range over $\mathcal T_n$, not $\Gamma^n$. *Fix:* insert Lemma~\ref{lem:support} and write
   $\tau\in\mathcal T_n$ in Definition~\ref{def:pointer}.

9. **[load-bearing] Theorem 3.1 (\texttt{thm:monotone}, l.397--408).** The hypothesis that the state $s_0$
   determining $\varphi$ is "computable from that window in $O(n)$ reads" is not part of the Definition of
   strategy, so the theorem proves membership for a subclass the Definition does not name. *Fix:* state as
   Theorem~\ref{thm:A} and add that the hypothesis is strictly stronger than being a rule.

10. **[load-bearing] Theorem 6.4 (\texttt{thm:upper}, l.819--821).** "A step is one evaluation of
    $\varphi$ on the current window at unit cost" is the oracle reading; under the computed reading the
    bound is $O(n^{3})$ \emph{evaluations}, hence $O(n^{4})$ reads, and for uncomputable $\varphi$ the
    theorem has no subject. *Fix:* restate as Theorem~\ref{thm:G}.

11. **[load-bearing] §6 (l.849--854, l.880--903).** The conservation law is stated at l.885 without proof,
    and the $p$-pass bound and the single-head bound are both stated in the body as if proved, while the
    rank hypothesis is "verified exhaustively at $n=8,10$" and $\sum_k\mathrm{CC}_k$ is $95$--$98$ per
    cent of its maximum \emph{by measurement}. *Fix:* state both as theorems with $(H_\Sigma)$ and
    $(H_{\mathrm{rank}})$ explicit (Theorems~\ref{thm:D} and \ref{thm:E}), and prove the law or cite it.

12. **[load-bearing] Abstract (l.45--46) against §4 (l.186--191).** The lower bound is for a
    \emph{single-head} machine; the class is defined for the \emph{two-head} machine of §2.
    $\rho\notin\mathsf{Pointer}_1$ does not give $\rho\in\mathsf{Walk}_2$. *Fix:* carry the head count in
    the class symbol and say in the abstract that the bound is for one head.

13. **[load-bearing] §7, Conjecture (l.566--570).** "Every predictor that reads the tape in $o(n^2)$ reads
    with $O(\log n)$ bits of memory has $R^2\le\kappa$": the hypothesis "cannot be computed in $O(n)$
    reads" is worst case, the conclusion is an average over $\mu_n$, and "predictor" is undefined.
    *Fix:* define the predictor class as in §5 and use $\mathsf{Walk}^{\mathrm{avg}}_2$ in the hypothesis.

14. **Table 5.1 caption and all §3.1 tables.** Every read count is "on tapes that produce an output", i.e.
    $\mathsf{Pointer}^{\mathrm{avg}|\mathrm{out}}$, whereas Definition~\ref{def:pointer} charges every tape
    including the $\bot$ ones. *Fix:* report the unconditioned mean alongside, or attach the conditioning
    to the definition and say that the two differ.

15. **§3.1 (l.202--206).** "The number of steps is measured on a finite grid rather than bounded" is exactly
    the gap in Theorem~\ref{thm:B}, and it is stated three paragraphs away from the place where the
    canonical rule is called pointer. *Fix:* move the caveat into the membership claim.

16. **§2 (l.128--129).** "It has no other memory, and in particular it cannot count its own steps" is false
    as written: the geometry records the number of left cuts ($L$), of right cuts ($n-1-R$) and the window
    length ($r-\ell$), so a rule can count cuts and can count moves since the last cut. What it cannot
    count is total moves. *Fix:* "its memory is the geometry $(L,R,\ell,r)$ and nothing else; it may
    therefore depend on the cut counts and the window length, but not on the number of moves."

17. **§2 (l.135) and Definition~\ref{def:run}.** The repeat clause presupposes that $(L,R,\ell,r)$
    determines the future, i.e. that $\varphi$ depends on the current window and not on the history. That
    is true of Definition~\ref{def:rule} but is nowhere stated as a Markov condition, and "a strategy may
    depend on $\ell$, $r$, $L$ and $R$ as well as on the window" (l.127--129) reads as if the list were
    merely permissive. *Fix:* add the Markov condition to Definition~\ref{def:rule}.

18. **§2 (l.122--124).** $f$ takes $n$ as an argument, so a rule is a uniform family $(f_n)$; elsewhere
    (l.107, "Fix $n\ge3$") a rule is treated as an object at a single $n$. *Fix:* say once that a rule is
    a family, and that this is what licenses a single device uniform in $n$ in
    Definition~\ref{def:pointer}.

19. **§5 (l.450--453).** "With cells held scaled by two, both bits become exact integer comparisons" is
    correct for reads, but each such comparison is $O(\log n)$ work-tape operations, so the \emph{step}
    count is $O(n\log n)$ and not $O(n)$. No reported number is wrong, but the cost model should say so.
    *Fix:* cite Definition~\ref{def:cost} and Remark~\ref{rem:arithmetic}.

20. **§6 (l.901--903).** "There is also no room for a space lower bound \dots\ $\omega(\log n)$ space is
    not merely unproved but false." Correct, and it deserves to be a theorem rather than a parenthesis,
    because it is the reason the conjecture must be about repeats. *Fix:* promote to
    Theorem~\ref{thm:F} and pair it with the $\Omega(\log n)$ floor of Remark~\ref{rem:space}, so that the
    space axis is closed at both ends.

21. **§7 (l.575--579).** The paragraph compares three different machines --- the rule (not a machine), the
    exact simulator ($5.7$ reads per step) and the statistical predictor ($\approx200$ reads per step) ---
    and says only that "the two figures measure different machines". *Fix:* name the three
    (Definition~\ref{def:rule}, Definition~\ref{def:pointer}, the predictor of §5) and never compare read
    counts without naming the machine.

22. **Terminology, throughout.** "Pointer" collides with the pointer machine of Tarjan; "memoryless" is
    used for rules whose geometry is four counters; "state" is used for the work tape, for $(L,R,\ell,r)$
    and for the $\varphi$-state. *Fix:* footnote the first, and reserve "state" for the device's work
    tape, "geometry" for $(L,R,\ell,r)$, "$\varphi$-state" for Definition~\ref{def:phistate}.

23. **§3 (l.188--190).** "A pointer strategy sweeps the tape and stops, while a walk strategy's headers
    return to cells they have already read, and the returning is what costs." This is contradicted by the
    paper's own examples: the interior-cut rule re-reads heavily and is pointer (l.238--244), and the
    median rule is linear per step and unclassified (l.273--283). *Fix:* delete, or mark explicitly as a
    heuristic the measurements do not support.

---

## 7. Insertion notes

* Replace the "Reads" paragraph (l.104--105) with Definitions~\ref{def:device}--\ref{def:cost} and
  Remarks~\ref{rem:uncharged}--\ref{rem:heads}.
* Replace the "Definition" paragraph (l.122--129) with Definitions~\ref{def:rule}--\ref{def:output},
  Definition~\ref{def:phistate} and Remark~\ref{rem:repeat}; keep the eight-action display and the run as
  they are, adding only clause 2 of Definition~\ref{def:run}.
* Insert §1 (Definitions~\ref{def:tape}--\ref{def:ensemble}, Lemmas~\ref{lem:cell}--\ref{lem:support},
  Remark~\ref{rem:support}) immediately after "Tape and measure" (l.87--92); the cell law is already there,
  the support is not.
* Replace the Definition of §3 (l.186--191) with Definitions~\ref{def:pointer}--\ref{def:avg} and
  Remarks~\ref{rem:relations}--\ref{rem:conv}. The two paragraphs that follow it (l.193--197) then reduce
  to one sentence: the partition is by Definition~\ref{def:walk}, the evidence is one-sided, and no rule is
  known to lie on the expensive side.
* Insert §5 as a short subsection, "What is proved, in the language of §4", at the head of §6 or as an
  appendix; the summary table is the part a referee will read.
* Nothing in this file requires a change to any measurement, table or number in the manuscript.
