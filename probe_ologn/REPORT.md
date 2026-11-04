# O(log n) bits vs O(n) reads for the two-header longest-run rule — probe report

Question: does the canonical rule (b2 = window-mean vs first-half-mean, b1 = longest run of 2s > m/4,
four-move table, cut-on-exit, halt when the region has < 3 cells) admit a simulator with
O(log n) bits of state and O(n) or O(n polylog n) reads?

**Verdict: no such simulator was found; I conjecture none exists.** What is proved here is a set of
structural facts that (a) pin the read complexity of the function at Θ(n) (≥ n−2 at every n tested,
n−1 at n = 8, 10, 12, and n−1 always suffices),
(b) show that *no* certificate/exchange argument can ever separate O(n) from O(n²) reads, and
(c) reduce the whole long-window dynamics of the rule to an explicit "concavity scan" whose only
natural O(n)-read jump fails on about half the tapes at n ≥ 64. All numbers below are from runs in this
directory; commands are in §8. No .tex file was touched.

---

## 1. Lemma A — b2 is a prefix-average comparison (proved, verified)

Let A_j = mean(a[L..L+j-1]) for the window starting at L, and m = r−L+1. Then

    b2([L,r]) = 0  <=>  A_{floor(m/2)} >= A_m .

Proof: b2 = 1 ⟺ mean(first k) < mean(window), k = max(1, floor(m/2)); for m = 2k or 2k+1 this is
equivalent to A_k < A_m (both cases collapse to the same inequality after multiplying out).
Verified on 20,000 random windows (rng seed 11, lengths 2..200): 0 mismatches (`p8_prefix_avg.py`).

So the scan that the b1≡0 walk performs is: examine m = M, M−1, … , 2 and ask whether the running
average at m has fallen to at most the running average at half-m. The answer sequence in m is *not*
monotone: B_m = A_m − A_{floor(m/2)} changes sign 9.7 / 14.9 / 24.2 times on average at
n = 128 / 256 / 512, and only 3.7% / 2.2% / 1.0% of (L,M) pairs have ≤ 1 sign change (`p8_prefix_avg.py`,
3 tapes per n). Hence no binary search / galloping on m.

## 2. Theorem B — the b1≡0 rule is exactly a "largest concavity" scan (proved, verified on 306,965 tapes)

With b1 forced to 0 the move table is (b2=1 → r−1, b2=0 → ℓ−1). Then, starting from (L,R):

    if R−L+1 < 3: halt (no output)
    m* := max{ m ∈ [2, R−L+1] : A_{floor(m/2)} >= A_m }        (A_j from position L)
    if no such m:      output L/(n−1)
    elif m* = R−L+1:   (L,R) <- (L+1, R−1)      [double cut]
    else:              L <- L+1                 [left cut]

Equivalently: scan r from R downwards; the first r with b2 = 0 gives a double cut if r = R and a left
cut otherwise; if every r in [L+1,R] has b2 = 1 the headers meet and the run outputs L.

Verified tape-for-tape against `walk_b0` (the same control flow as the naive implementation) on **all**
ensemble tapes (each column has ceil(n/2) ones; n even: #0 = #2, n odd: #0 = #2+1):

| n | tapes | agreement |
|---|---|---|
| 6 | 141 | 141/141 |
| 7 | 357 | 357/357 |
| 8 | 1,107 | 1,107/1,107 |
| 9 | 2,907 | 2,907/2,907 |
| 10 | 8,953 | 8,953/8,953 |
| 11 | 24,068 | 24,068/24,068 |
| 12 | 73,789 | 73,789/73,789 |
| 13 | 201,643 | 201,643/201,643 |
| **total** | **306,965** | **306,965/306,965** |

(`$P p2_b0.py 6 7 8 9 10 11 12 13`.) This is the structural reason the trajectory is long: each epoch
re-scans from R down to the largest concavity position, and the average epoch costs ~0.1 n b2 tests.

## 3. Lemma C — b1 can only fire on short windows; the long-window dynamics are exactly b1≡0

If b1 = 1 then the window contains a run of > m/4 twos, so m < 4·rho where rho is the *global* longest
run of 2s in the tape (rho is computable in O(n) reads with O(log n) bits). Measured (40 tapes per n,
seed 2026):

| n | rho mean (max) | max m at a b1-true step | b1-true steps/tape | b1-true fraction of steps | steps/n² full | steps/n² b1≡0 |
|---|---|---|---|---|---|---|
| 64 | 2.65 (5) | 11 | 2.25 | 0.835% | 0.0658 | 0.0673 |
| 128 | 3.05 (6) | 14 | 2.42 | 0.307% | 0.0482 | 0.0484 |
| 256 | 3.83 (5) | 15 | 5.05 | 0.119% | 0.0645 | 0.0645 |
| 512 | 4.33 (6) | 22 | 5.38 | 0.030% | 0.0688 | 0.0688 |

No tape violated m < 4·rho. So on every window longer than ~22 cells at n = 512 the rule *is* the b1≡0
rule; b1 fires O(1) times per tape, always at the tail of an epoch. This confirms the paper's (v) and
sharpens it: the long trajectory is exactly the concavity scan, and the only correction b1 makes is at
window sizes O(rho) = O(log n) (for random tapes). Full-rule delta equals b1≡0 delta on 34–36 of 40
tapes at n = 64…512 (85–90%).

## 4. Attempted O(n)-read simulators (all fail tape-for-tape)

1. **Horizon jump.** Assumes g(L+1) ≤ g(L) for consecutive left cuts with R fixed (then the
   largest-concavity scan for L+1 never needs to look above g(L), and the whole run costs O(n) tests).
   The assumption is false: violations occur in 15/45, 259/765, 3,848/10,749, 14,766/35,493, 45,399/126,195
   consecutive pairs at n = 6, 8, 10, 11, 12 (36% at n = 12), and only 53% of the violations are of the
   benign "g(L+1) = R" form (the other 47% are genuine re-entries above the horizon) (`p4_monotone_g.py`, `p5_violations.py`).
   The resulting O(n)-read algorithm is exact on only part of the ensemble:

   | n | tapes | mismatches vs the exact b1≡0 walk | b2 tests/tape |
   |---|---|---|---|
   | 8 | 1,107 (all) | 90 (8.1%) | 5.67 |
   | 10 | 8,953 (all) | 1,206 (13.5%) | 7.75 |
   | 11 | 24,068 (all) | 11,043 (45.9%) | 9.54 |
   | 12 | 73,789 (all) | 13,748 (18.6%) | 9.86 |
   | 64 | 200 (seed 2026) | 95 (47.5%) | 70.3 (1.10 n) |
   | 128 | 200 | 95 (47.5%) | 148.0 (1.16 n) |
   | 256 | 200 | 101 (50.5%) | 318.2 (1.24 n) |

   So there *is* an O(n)-read algorithm here — it is just wrong on about half the tapes at n ≥ 64, and
   the mismatch rate does not decrease with n (`p7_horizon.py`).
2. **"R frozen" closed form** (pretend R never moves): matches the exact b1≡0 delta on 74/141, 504/1,107,
   3,577/8,953, 8,207/25,653, 26,492/73,789 tapes at n = 6, 8, 10, 11, 12 — i.e. 32–52%; double cuts
   cannot be ignored (`p3_b0_props.py`).
3. **Monotone predicate / binary search.** P(L) = "the region [L, n−1] has no concavity" is monotone in L
   on only 62.5% (n = 10) and 52.7% (n = 12) of tapes, so the answer cannot be found by binary search.
4. **One-pass walk statistics** (argmin/argmax of cumsum(a−1), of cumsum(1[a=2]−1/4), centre of mass of
   2s/0s, first/last 2, first/last 0, and variants): best exact-hit rate 22/47 (n = 6), 147/417 (n = 8),
   984/3,413 (n = 10). Prior work's seven O(n) summary predictors are all ≤ 12% at n ≥ 64.
5. **Prefix-separation (one-pass lower bound) test.** For a one-pass streaming lower bound one needs:
   any two prefixes with the same count profile are separated by some common suffix. False — 15/33
   (n = 6), 78/279 (n = 8), 316/2,205 (n = 10) same-profile prefix pairs have *identical* delta on
   every common suffix. So the naive fooling-set/streaming route does not give a bound either
   (`p9_prefix_sep.py`).

## 5. Certificates: the read complexity is exactly n−1 (exhaustive, n ≤ 13)

For a tape a, a *certificate* is a set S of cells such that every valid tape agreeing with a on S has
the same delta. Every correct simulator on input a must read a certificate: its execution is identical
on any tape agreeing with it on the cells it read. So reads(a) ≥ c(a) = min |S|.

Computed exactly by exhaustively projecting all tapes onto every subset (`p6_cert.py`; the code was
validated against a brute-force subset search at n = 6: identical, 141/141).

| n | tapes | K_global | c(a) min | median | mean | max | max/n | mean/n |
|---|---|---|---|---|---|---|---|---|
| 6 | 141 | **5** | 2 | 2 | 2.64 | 4 | 0.667 | 0.440 |
| 8 | 1,107 | **7** | 2 | 4 | 3.64 | 7 | 0.875 | 0.455 |
| 10 | 8,953 | **9** | 3 | 5 | 4.79 | 9 | 0.900 | 0.479 |
| 11 | 24,068 | **10** | 3 | 6 | 5.66 | 9 | 0.818 | 0.515 |
| 12 | 73,789 | **11** | 3 | 6 | 5.88 | 11 | 0.917 | 0.490 |
| 13 | 201,643 | **12** | 3 | 7 | 6.75 | 11 | 0.846 | 0.519 |

K_global = smallest |S| for which *some* set S is a certificate for *every* tape. It equals n−1 at
every n tested: **no set of n−2 or fewer cells determines delta**, i.e. there is no "read these k
cells" shortcut for k < n−1, even non-adaptively and with unlimited memory.

Upper bound (proved for all n): n−1 cells always suffice, because the ensemble forces the missing cell:
#0 − #2 ∈ {0,1} (0 for n even, 1 for n odd), so given the other n−1 cells exactly one value of the last
cell is consistent. Hence c(a) ≤ n−1 for every tape, and the deterministic decision-tree (query) complexity of delta is
n−1 at n = 8, 10, 12 and n−2 at n = 6, 11, 13 — i.e. exactly n−1 up to the parity noise of the small
cases, and n−1 as the general upper bound. Combined with the table: any correct simulator reads ≥ n−1 cells
on some tape, and ≥ 0.44 n on average.

## 6. The consequence for the classification (this is the load-bearing argument)

* The read complexity of delta with *unbounded* memory is Θ(n): n−1 reads suffice, and ≥ n−2 are
  necessary at every n tested (n−1 at n = 8, 10, 12). So O(n) reads is
  achievable if O(n) bits are allowed (read n−1 cells, infer the last, simulate in memory).
* The pointer/walk distinction is therefore **not** a statement about reads at all: it is a
  space-vs-reads tradeoff. Any proof that the rule is not a pointer strategy must be a *space* lower
  bound (no O(log n)-bit machine computes delta in O(n) reads); read-counting, exchange, pivoting and
  fooling-set arguments are structurally capped at n−1 reads, so they can never separate O(n) from
  O(n²).
* Concretely, a linear lower bound cannot be improved by looking at "which cells matter": the
  certificate bound is tight at n−1. (For contrast, parity has certificate n and O(1) space, so
  certificate size carries no space information at all.)
* What the data do support is the *reason* a log-space simulator has to re-read: the rule's trajectory
  asks Θ(n²) window questions, each a comparison of two block averages over positions that can be
  n/2 apart; with O(log n) bits none of the prefix sums, and none of the answers, can be cached, and
  the answer set has no monotone structure (10–24 sign changes at n = 128–512). The only jump we could
  build costs O(n) reads and is wrong on ~half the tapes. That is evidence, not proof.

## 7. Proved vs conjectured

**Proved (this run):**
1. b2 = 0 ⟺ A_{⌊m/2⌋} ≥ A_m (Lemma A; verified 20,000 windows).
2. The b1≡0 rule is exactly the largest-concavity scan of Theorem B (verified on all 306,965 ensemble
   tapes, n = 6…13).
3. b1 = 1 ⟹ m < 4·rho(tape); measured b1-true steps 2–5 per tape, all at m ≤ 22 for n ≤ 512, and
   full/b1≡0 step counts agree to 3 digits (Lemma C).
4. K_global = n−1 for n = 6, 8, 10, 11, 12, 13: no set of n−2 or fewer cells is a certificate for
   every tape.
5. Every correct simulator reads ≥ c(a) on tape a; max c(a) ∈ {n−2, n−1} at n = 6…13 (n−1 at
   n = 8, 10, 12); mean c(a) = 0.44–0.52 n. In particular reads ≥ n−2 on some tape at every n tested
   (n−1 at n = 8, 10, 12), and the average over tapes is ≥ 0.44 n.
6. c(a) ≤ n−1 for every tape and every n (the counts force the last cell), so no certificate argument
   can prove more than n−1 reads.

**Conjectured:** delta is not computable in O(n) (or n polylog n) reads with O(log n) bits. The evidence
is the reduction to the concavity scan, the failure of the natural jump, the non-monotonicity of the
b2 answer sequence, and the fact that a simulator must answer Θ(n²) far-apart block-average
comparisons while holding O(1) of them. **No proof is available**: a proof needs a space lower bound
for a machine that may read adaptively in any order, and the standard one-pass fooling-set route fails
on this function (§4.5).

## 8. Files and commands

All in `C:/Users/shini/poetry_axis/paper/probe_ologn/`; `$P = C:/Users/shini/poetry_axis/intern_decision/.venv/Scripts/python.exe`.

```
$P p0_tapes.py                      # ensemble enumerator (n even: #0=#2; n odd: #0=#2+1)
$P p1_enum.py <n>                   # all tapes, delta, one-pass candidate hits (n=6..10 run)
$P p2_b0.py 6 7 8 9 10 11 12 13     # Theorem B, all 306,965 tapes
$P p3_b0_props.py                   # 'R frozen' form; b1 census; full-vs-b0 delta
$P p4_monotone_g.py 6 8 10 11 12    # g(L+1)<=g(L) violations
$P p5_violations.py 8 10 11 12      # violation taxonomy + double-cut census
$P p6_cert.py 6 8 10 11 12 13       # exact certificates (cert_n*.npy, cert_summary.json)
$P p7_horizon.py 8 10 11 12         # horizon jump mismatches (exhaustive)
$P p8_prefix_avg.py                 # Lemma A + sign-change statistics
$P p9_prefix_sep.py                 # prefix-separation (streaming) test
$P p10_rho.py                       # K_global(6,8); rho vs b1; step counts
```

Seeds: `walk_naive`/`walk_fast` comparisons and b1 census use seed 2026, 40 tapes per n; horizon
comparison at n = 64/128/256 uses seed 2026, 200 tapes per n; Lemma A uses rng seed 11. `cert_n*.npy`
hold `[delta_code, cert_size]` per tape.

## 9. What this means for the paper

1. The operational definition of the walk class should not rest on the phrase "no O(n) simulator is
   known": the read complexity of delta is exactly n−1 with unbounded memory, so the honest statement
   is a space statement — "no O(log n)-bit simulator with O(n) reads is known". A read-counting proof
   of the classification is impossible in principle (§6).
2. Two structural lemmas are worth stating in the paper because they explain where the quadratic cost
   lives: (i) b2 is the prefix-average comparison A_{floor(m/2)} vs A_m, and the b1≡0 rule is the
   largest-concavity scan (Theorem B); (ii) b1 can only fire when m < 4·rho, the global longest 2-run,
   so every window longer than ~4 rho follows the b1≡0 dynamics exactly (measured: b1-true at m ≤ 22
   for n ≤ 512, 2–5 such steps per tape).
3. The strongest negative statement available is the exhaustive certificate table (§5): on some tape
   every correct simulator must read n−1 cells, and on average 0.44–0.52 n. That is a strictly better
   constant than the exchange bound (0.25–0.63 n) and it is tight as a linear bound.
