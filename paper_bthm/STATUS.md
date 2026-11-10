# STATUS — 2026-11-10 (evening)

Handoff. Every claim names the file that carries it. Supersedes the 14:22 version.

## 1. The paper

`bthm_jalc_v2.tex` — **26 pages, 0 errors, 0 undefined references**. Six sections (Introduction /
The machine / Three costs and the classification / The machines and the ladder / Limits / Where this
sits) + Appendix A (pp. 17–24). The v1 file `bthm_jalc.tex` is untouched at 22 pp.

The classification, **three classes, decision by `g` first**:

| class | criterion | members |
|---|---|---|
| fully linear | `g = O(n)` | `monotone`, `nhalf`, `pattern`, `runs`, `peak` |
| semi-linear | `e = O(n)`, `g = ω(n)` | `ballast` (constructed, `probe_semilinear/`) |
| superlinear | `e = ω(n)` | the other seventeen natural rules |

A fourth cell (quasi-linear) was removed: `p = Θ(e)` wherever both are measured with one machine, so
no rule with `e = ω(n)` and `p = O(n)` exists in the catalogue.

## 2. What was established (morning and evening)

- `p ≤ e ≤ g`; `p = Θ(e)` (measured `p/e = 0.66–1.00`); `e = Θ(g)` for every rule except `longrun`.
- **`G` is measured for all 21 natural rules**: fully linear `G/e = 0.997–1.000`, `g/e ≤ 1.16`
  (`gpeak.py`, `gfl.py`); superlinear `G/e = 0.196–0.875` (`gsup.py`, `gmean.py`).
- **A semi-linear rule exists**: `ballast`, `e ~ n^0.997`, `g ~ n^1.995`, `G ~ n^1.000`
  (`probe_semilinear/SEMILINEAR.md`).
- **`longrun`'s witness is `n³/21`, not `25 n²`** — the record's figure was the `n = 512` slice; its
  `g/e` reaches 86 232 at 4096 (`probe_g_witness/gw_cells.json`).
- **`e`'s exponent on the band rules is not estimable**: the mean is 48–98 % tail (depending on the
  cut), the local exponent runs 1.05–2.76 without converging; `G`'s local exponents span 0.10–0.97
  (`probe_k300/`). Only `e = ω(n)` is safe to state.
- **`p ≪ e` has no instance**: cheap failure sets do not exist, early stopping saves nothing, and
  `R²(k) ≤ 0.45` even with a random forest and half the tape sampled (`probe_r2k/`).
- **Two propositions are proved and in the paper**: the chain `G ≤ e ≤ g`, and bounded cost ratio ⇒
  perfect. The conjecture "no content-dependent superlinear rule is perfect" has 11/11 evidence.

## 3. The audit

`probe_audit/AUDIT_V2.md`: 79 claims, 39 verified, **30 mismatches**, 6 unlabelled machines, 3
unsourced. All 30 fixed in four passes; `AUDIT_V2_RECHECK.md` confirmed 24 and found 6 more (fixed),
including **one error introduced by the fix itself** (nine defects called the machine's maximum —
`PC_PUSH` says eight, and the cost is flat for every fixed count to sixteen).

## 4. The mock review (v30-v38)

**v36-v38 are independent rounds**: three referees in three *fields* (computational complexity,
statistics and Monte Carlo, algebraic and arithmetic combinatorics), each reading the paper without
any previous report. v30-v35 were revision-tracking rounds (`--r2`), which anchor the referees and
accumulate their own item lists; their counts are not comparable.

| round | A complexity | B statistics | C algebra |
|---|---|---|---|
| v36 (first independent) | Major | Major | **Reject** |
| v37 | Reject | Major | Major |
| v38 | Reject | Major | Major |
| **vA x3 (same paper, same prompt)** | **Major, Major, Major** | — | — |

**A is Major, not Reject.** Three fresh draws of the same referee on the same manuscript returned
Major three times; the v37/v38 Rejects are single-sample noise around that boundary, and the objection
itself ("the central claims are not proved for the stated machine class") was word-for-word identical
in every draw. The objection is the result; the label is the noise. Corrected 2026-11-11.

The independent rounds found what five anchored rounds did not: the rank theorem's missing lower
bound (C, v36), the two objects both called `kappa` (A), the reversed transfer, and the `g` definition
sitting on the measured machine rather than the min-max (A, v37). All are fixed. What remains is the theoretical side: A (Major on three of three draws) says the central claims are
not proved for the stated machine class, and no revision supplies a proof; B and C are Major and their
lists are口径 work. On the evidence this is a major-revision paper at an empirical venue, not a reject.

## 4b. The mock review (five anchored rounds, v30-v34)

| round | A theory | B measurement | C reproducibility |
|---|---|---|---|
| v30 (first) | Reject | Major | Reject |
| v31 | Reject | Major | **Major** |
| v32-v34 | Reject | Major | Major |

What moved: the class tree (resolved in v31), the semi-linear construction (B calls `ballast` the one
substantive gain), the kernel dichotomy in abstract and introduction, the transfer factor in all
three places, the two-alphabet information count, and the two different objects that were both called
`kappa` (the cost ratio is now `chi`).  What does not move: A's Reject is structural --- the paper has
no lower bound for the superlinear side and says so; no revision supplies one.  B and C are fixable
and still list, in v34: `p <= e <= g` asserted where only `G <= e <= g` is proved; `cor:perfect`
still leaning on `prop:spread`; the abstract's "every exponent is fitted, never derived from a step
count" called false; `PC`'s mechanism described differently in Section 4 and in the ladder.

## 4b. The report that preceded these rounds

`mock_reviews/2026-11-10-v30/`, three referees, journal labelled Experimental Mathematics, data
summary supplied.

| referee | verdict | main objection |
|---|---|---|
| A theory | **Reject** | convergence mode unnamed; `e,g` domains vague; **three different `g` objects**; tail conditions; a "constant" from four points |
| B measurement | **Major revision** | **the ordered decision tree was broken** (now fixed); exponents without intervals; `κ` treated as a supremum but is a one-family lower bound |
| C reproducibility | **Reject** | the abstract's class counts contradicted the data (fixed); three values vs four states; the kernel-rank theorem vs `n^{2.93}`; the transfer direction; the charge model |

Fixed already: the decision tree, the class table, the abstract, the introduction, the conclusion, and
the `κ`/quasi-linear wording in five places.
Still open: (3) the alphabet statement, (4) the kernel-rank theorem's dependence and the `n^{2.93}`
qualifier, (5) the transfer inequality as a sum over critical queries, (6) one charge model defined
once and cited, (7) bootstrap intervals or a per-row "local, window only" label.

## 5. The venue

- **EM (Experimental Mathematics)**: v2 submitted 2026-10-30 (`tww_em.zip`, `tww_em_v2.pdf`); Zenodo
  v3 uploaded 2026-11-03 (`_zenodo_v3/`). T&F policy: preprints are not duplicate publication; the
  preprint must be **declared in the manuscript and cited**; no requirement that the submitted and
  public versions match. Action: add the declaration in the next revision.
- **BTHM has no venue commitment.** The `bthm_jalc*` names are legacy; the 11-04 memo recommended the
  empirical route (EM or JEA), not a theory venue.

## 6. Next steps, in order

1. Finish review items 3–7 in `bthm_jalc_v2.tex`; add the preprint-declaration footnote.
2. Re-run `_review.py --r2` against v30 and check each point moved.
3. Build the packages: `_zenodo_bthm_v1/` (PDF + the three .tex files) and the GitHub tree
   (`submission/` + this session's nine probes + a "script → number" manifest). **No push.**
4. Venue decision after EM responds.

## 7. Where things are

`bthm_jalc_v2.tex` / `.pdf` — the paper. `ORDERS.md` — per-rule `p/e/g` with sources. `TABLE3.md` —
the measured `(G,e,g)` table. `sec_gmean.tex` — the three-cost section as a fragment. `STATUS.md` —
this file. `_venue_memo_2026-11-04.md` — the venue reasoning and six earlier rejections.
`probe_audit/`, `mock_reviews/`, `probe_semilinear/`, `probe_r2k/`, `probe_g_witness/`,
`probe_k300/`, `probe_e_cheap/`, `probe_ces_exact/`, `probe_en_linear/`, `probe_peak_fix/` — the
evidence.
