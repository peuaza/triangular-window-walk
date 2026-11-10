# CHANGELOG_CLAIMS.md

**The rule.** No edit to the paper is finished until (1) `check_paper.py` passes, (2) every headline
number the edit touched is still in `claims.json`, and (3) the edit is logged below with the claim ids
it touched. An edit that changes a statement but not this file did not happen.

`build.sh` enforces (1) and refuses to compile otherwise. (2) and (3) are on the person editing.

Format: `date | what changed | claim ids touched | checker state`.

---

| date | what changed | claim ids | checker |
|---|---|---|---|
| 2026-11-10 | the action renamed `cut` -> `trim` throughout; the companion-paper paragraph added | all (the word appears in every cost sentence) | NO HITS |
| 2026-11-10 | the class taxonomy cut to three (FL / SL / SU); `WS`/`SS`/`QL` retired in text and in the catalogue table | class_census | NO HITS |
| 2026-11-10 | the catalogue lost `S_n`, `Lambda_2(c=2)`, `R_+` (rows from the earlier paper) and gained `ballast` | class_census | NO HITS |
| 2026-11-10 | `cor:perfect` relabelled: the spread proposition is an identity, its content is the measurement | g_over_e_fl, g_over_e_super | NO HITS |
| 2026-11-10 | the two objects both called `kappa` separated: cost ratio is now `chi` | kappa_range, kappa_star | NO HITS |
| 2026-11-10 | the kernel theorem's conclusion propagated to the abstract and the introduction (`n^3/log n`, not `n^3`) | cx_order | NO HITS |
| 2026-11-10 | `n^{1/log n}` removed from the exponent conjecture: it is a constant | (none numeric) | NO HITS |
| 2026-11-10 | the information ceiling's direction fixed: the smaller bound is the stronger | info_ceiling_sum, info_ceiling_four | NO HITS |
| 2026-11-10 | the status table's `Omega(n^2)` floor relabelled as refuted, with the `n^{3/2}` rung and `ballast` as witnesses | class_census | NO HITS |
| 2026-11-10 | the five fully-linear means stated in Section 3 (they were only in the appendix) | monotone_mean, pattern_ces, runs_ces, peak_exact | NO HITS |
| 2026-11-10 | Figure 1 added: per-tape cost densities, `monotone` + the eleven band rules | g_over_e_super | NO HITS |
| 2026-11-10 | the transfer factor unified in three places (multiplication by `f`, slope `s - 0.99`) | transfer_slope | NO HITS |
| 2026-11-10 | `ballast` construction written into Section 3 | ballast_mean, ballast_exponents, ballast_ratios | NO HITS |

| 2026-11-10 | **notation pass**: `Theta` kept only where both bounds exist. The kernel theorem is now a theorem of upper bounds plus `conj:kernel-floor` for the missing lower one; the abstract, the introduction and the rank paragraph say `O(n^2)` / `O(n^3/log n)`; the endgame buffer is `Omega(n)` forced and `Theta(n)` sufficient; fitted quantities say "measures as", not `Theta` | cx_order | NO HITS |

| 2026-11-10 | `g` split into two objects: the rule's min-max `g^mm` and the reported `g` on the e-best machine `M*`, with `g^mm <= g` stated and the converse denied; the reversed sentence removed; the spread proposition uses `chi` throughout | class_census, g_over_e_fl, g_over_e_super | NO HITS |
| 2026-11-10 | the collapse split into `prop:collapse` (verified on 320,880 queries) and `conj:collapse` (general n); the appendix aligned | monotone_mean | NO HITS |

## Claims deliberately absent from the paper

Recorded here so a later pass does not "restore" them:

- `0.62 n ln n` for `hashmod`: rows, not cell reads (the cell count is `0.093 n^2`).
- `25 n^2` for `longrun`: the `n = 512` slice; the cost is `n^3/21`.
- `0.12 n^2` for `runs`: the pre-CES machine; the paper's figure is `2.77 n`.
- `G/e = 0.196-0.875`: the six-rule supplement, not the paper's table.
- `1497/1500` for `peak`: the corrected count is `1496/1500`.
- `kappa in [0.196, 0.432]`: the pre-INC range.
