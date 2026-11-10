# AUDIT_V2_RECHECK — did the fixes land?

Re-check of `AUDIT_V2.md` after the author's edits to `bthm_jalc_v2.tex` (948 lines) and
`_appendix_evidence.tex` (473 lines). Scope: the **30 MISMATCH** items of AUDIT_V2 §E, the
**6 MACHINE-UNLABELLED** sites (§F), the **3 NO SOURCE** items (§G), plus the four new checks the
brief asks for. Line numbers below are the current ones.

Every "regenerated" value was recomputed in this pass from stored data — the new prose was not
trusted either. Sources used: `probe_en_linear/el_*.json` (per-tape `steps`, all runs, mean over
all drawn tapes), `probe_hashmod/hm_cost_all.json`, `probe_k300/k3_*.json` (`reads` arrays),
`probe_kappa_inc2/ki2_five.json` + `ki2_agg.json`, `probe_kappa_inc/ki_agg.json`,
`probe_g_witness/gw_fit.json`, `probe_ces_exact/ce_patch_*.json`, `probe_peak_fix/gpeak.log`
(a fresh `gpeak.py` run exceeded the 180 s cap; the log is newer than the script),
`HYBRID3.md`, `LONGRUN_EREMEASURE.md`, `PC_PUSH.md`, `probe_cx_down/CX_DOWN.md`. Arithmetic on the
JSONs was redone in Python; `probe_r2k/gfl.py 60` was re-run.

## Verdict summary

| group | fixed | not fixed |
|---|---|---|
| 30 MISMATCH | **24** | 5 STILL WRONG (A1, A5, D8, D9, D14) + 1 NEW ERROR (B12) |
| 6 MACHINE-UNLABELLED | 3 (M1, M3, M4) | 1 partial (M2) + 2 open (M5, M6) |
| 3 NO SOURCE | 1 sourced (C14) | C16 still unsupported; C9 rephrased into a **new false claim** |
| new checks | κ table ✓, abstract/§1/§2.1 consistent ✓ | exponents without a window: 3 places |

The headline: the κ table, the step-count low end, the epoch fraction, the tail cuts, the median,
the witness counts, the longrun law and the level-count table were all genuinely repaired. What the
fix missed is concentrated in **one load-bearing sentence about κ\***, **one headline number
(15.8)**, **two stale scan/ CES values**, and **a newly invented justification for "nine
defects"**.

---

## 1. The 30 MISMATCH items

"old" = value printed before the fix (from AUDIT_V2 §E); "current" = what the file says now.

| # | old | current | regenerated in this pass | verdict | where |
|---|---|---|---|---|---|
| A1 | 0.017–15.8, "factor of nine hundred" | 0.385–15.8, "factor of forty" | 0.385 (NBLK>3, 160-tap el sample, n=4096, all endings) → **15.576** (hashmod seed 2026); **12.541** at seed 777; 15.8 exists in no JSON | **STILL WRONG** (upper end 15.8 vs 15.58; 1.4 % high, seed-dependent) | 38–40, 83–85, 210–212, app 369–372 |
| A2 | 1.50–2.60 | 1.70–2.20 | 1.702 (KCHG) – 2.203 (oddrun), all-runs means, 1024→4096, 15-rule el grid | **FIXED** | 40, 212, app 372 |
| A5 | "at most 1.12 n²" (bound) | 1.33 n² in abstract, §5, §9; **ladder row still "≤1.12 n²"** | HYBRID3 §4.3: 1.33 = 5.3 × 0.2505; max observed e/n² = **0.40**; 1.12 = plateau constant 1.1195 | **STILL WRONG** (ladder row; and "measured 1.12" reads as an observed e/n²) | 50, 633, 867, **681** |
| A13 | 0.093 n² | 0.08–0.096 n² | 0.0801 (seed 2026) / 0.0957 (777), master 0.0956 | **FIXED** | 375–376 |
| A15 | 4.87+0.45 ln n | 7.02+0.15 log₂ n (asympt. 4.8+0.31) | LONGRUN_EREMEASURE: 7.021+0.1494 log₂ n over 512–16384 | **FIXED** | 379–380 |
| A18 | 1.3–274, "eight of nineteen" | 1.3–274, "six of the eighteen" | 6 of 18 lift ≥1.3 (1.346–273.7) | **FIXED** | 362–363 |
| A21 | 0.3–36 % | 0.06–17 % | 0.06 % (oddrun 4096) – 17.10 % (dense2 4096) | **FIXED** | 314 |
| A23 | 820 both | 820 (AMAX), 616 (AMIN) | `_tables.md` M5: 820 / 616 | **FIXED** | 318 |
| A25 | 73–98 under "10× median" | 100 n: 73–98; 10× median: 48–97 | 4096: 100 n = 0.727–0.980; 10× median = 0.482–0.972 | **FIXED** | 322–323, 350–351 |
| A26 | "median 5.3 n" rule unnamed | 5.3 n for BAL>m/2 and KCHG; eleven-rule range 4.1–145 | median/n 5.26 / 5.29; range 4.10–145.11 | **FIXED** | 323–324 |
| A29 | 0.6–1.5, no depth | 0.55–1.9 at 25 % depth | a@25 %, n=4096: 0.55 (AMAX) – 1.92 (oddrun) | **FIXED** | 351 |
| A31 | g/e = 3.4–4.1 | 3.4–6.2 across K=2–16 | 3.38–6.21 | **FIXED** | 347 |
| B3 | 200 meeting tapes at four sizes | 200 at n≤2048 and 100 at n=4096 | nt_meet 200/200/200/**100** | **FIXED** | 306 |
| B7 | pattern g/e ≤ 1.14 | pattern ≤ 1.11; runs ≤ 1.16 | re-ran `gfl.py 60`: pattern 1.092–1.113; runs 1.123–1.153 | **FIXED** | 536 |
| B8 | 2.03→1.07 over 512→2048 | 2.03→1.05 over 512→4096 | `gpeak.log`: 2.032 → 1.048 (G/e 0.984→1.000) | **FIXED** | 538–539 |
| B11 | "fourteen = three + eleven" | eleven rows + the three perfect rules | table has 11 rows; 3 perfect rules live in §3 | **FIXED** | 603–605 |
| B12 | "all eleven have e/G ≥ 4" | "e/G is 3.8 to 12 with one exception at 3.85" | e/G = 3.85 (KCHG) – 12.05 (AMAX); **3.85 is the minimum, not an exception** | **NEW ERROR** (self-contradictory prose) | 604 |
| C1 | 4.48–4.50 | 4.36–4.50 | HYBRID3: 4.36 → 4.48 → 4.50 | **FIXED** | 623 |
| C4 | 4.41–4.86 | 4.41–4.87 | HYBRID3: 4.41 → 4.87 → 4.86 | **FIXED** | 626 |
| C6 | 0.87 attributed to AMAX | 0.655 ramp / 0.87 plateau / 1.12 extrema | HYBRID3 §4.4 exactly | **FIXED** | 633–637 |
| D1 | trip 0.799→0.510 | 0.789→0.472 | `ki2_five` c16 R²: 0.789→0.472 | **FIXED** | 757 |
| D2 | peak 0.782→0.585 | 0.785→0.614 | 0.785→0.614 | **FIXED** | 758 |
| D3 | pattern 0.677→0.502 | 0.693→0.502 | 0.693→0.502 | **FIXED** | 759 |
| D4 | dense2 0.626→0.517 | 0.632→0.509 | 0.632→0.509 | **FIXED** | 760 |
| D7 | "largest 0.927" | 0.999 at n=512, 0.927 at n≥1024 | `ki_agg` inc_r14_c128: max 0.9985 (trip 512); 0.927 = trip 1024 | **FIXED** | 747, 767 |
| D8 | scan 0.119–0.344 | 0.119–0.344 (unchanged) | five rules, all sizes: **0.119–0.361** (pattern 1024); at n=4096: 0.123–0.344 | **STILL WRONG** (range splices n=512 min with n=4096 max) | 777 |
| D9 | caption 0.196–0.432 | caption silent; **appendix status table still prints κ∈[0.196,0.432] "for fifteen of sixteen rules"** | no column gives that pair (50n scan: 0.124–0.458); and there is no 16-rule grid | **STILL WRONG** (survives in the appendix) | app **450** |
| D10 | 140–890 meeting tapes | 373–1 835 | `ki2_five` n=512 meet: 373, 467, 1107, 1197, 1835 | **FIXED** | 784 |
| D14 | "R² = 1 at 3.3–4.5" | caption: "R²≥0.99 at 3.4–4.1 on runs and pattern"; **"3.3 reads per cell" still at 787 and 803** | p-arm (R²≥0.99): pattern 3.58–3.68; runs **58.07 at 512, 3.40 at 1024, 8.45 at 2048, 3.49 at 4096** | **STILL WRONG** (two values for one device; runs not in 3.4–4.1) | 769–770, **787**, **803** |
| D17 | 1.12 repeat in §6 | 1.33 n² (measured 1.12) at §9 | HYBRID3 1.33 | **FIXED** (ladder residue counted under A5) | 867 |

Count: 24 FIXED, 5 STILL WRONG, 1 NEW ERROR = 30.

---

## 2. The 6 MACHINE-UNLABELLED sites

| # | site | current state | verdict | where |
|---|---|---|---|---|
| M1 | `tab:threecosts` (11×6 numbers) | caption now: "Every cost is read from the `hmach3.H3` band device at B=64, charged one read per cell taken or given up plus the opening pass." | **FIXED** | 578–583 |
| M2 | Table 3 (level-count rules) | machine column names `band + frontier, B=64`; charge named only in the following paragraph ("charged one read per cell taken or given up plus the opening pass"). **HYBRID3 §2 records the NBLK>2 and NBLK>3 rows as `band+frontier+certifier`**, so those two rows still name the wrong device. | **PARTIAL** — charge adjacent, device wrong for 2 of 4 rows | 619–631, 633–638 |
| M3 | kappa table | caption now names the ridge regression, the incremental charge, "one to four reads per step, not m"; 41 bits in the prose below | **FIXED** | 764–773, 786 |
| M4 | §5 "1.12 n²" | device + tape family now named: "on the band + frontier device … 0.655 for NBLK=1 on ramp, 0.87 on plateau, 1.12 for the two extrema" | **FIXED** | 633–637 |
| M5 | abstract "3.4–5.3 / 1.12 n²" | abstract now says 1.33 (measured 1.12) but still names **no device** | **NOT FIXED** | 49–51 |
| M6 | ladder table | caption promises "the device that attains each rung named in the third column", but the column holds mechanisms ("Brent cycle detection", "A=s/31 anchors", "a witness tape") and no probe pointer. Worse, the n² row lists AMAX/AMIN under `band + frontier` without the certifier that HYBRID3 says they need. | **NOT FIXED** | 671–695 |

---

## 3. The 3 NO SOURCE items

| # | item | check | verdict | where |
|---|---|---|---|---|
| C16 | "seventeen rules between 0.50 and 0.80, hashmod ≈ 0" | now "**thirteen** rules between 0.47 and 0.80 and hashmod at approximately 0". `ki2_agg` has 13 rules; at n=512, c=16 exactly **9** lie in [0.47, 0.80], **3 above** (NBLK>3 0.819, NBLK>2 0.826, canonical 0.925), 1 below (hashmod 0.260). At n=2048 BAL>m/2 is −0.097 and canonical 0.871. The 0.47 is the **five-rule** min (0.472) and 0.80 the five-rule max (0.789): the sentence splices a five-rule range onto a 13-rule count, and 13 + hashmod = 14 rules in a 13-rule grid. | **STILL WRONG** — and if the grid is the evidence, its maximum is 0.925 (0.826 excluding canonical), so "between 0.47 and 0.80" is not the range that bounds κ\* | 730–733 |
| C14 | A = s/31 anchors, one 31-bit residue | `CX_DOWN.md` l.49/51: `D = ⌊s/31⌋`, "`31 D` bits"; l.278/313 same. | **SOURCED** — minor naming only: the report uses `D`, the paper `A`; "anchors" is not a term in the report | 683, 705–706 |
| C9 | "nine defects" | body 644 is now sourced (`PC_PUSH.md` l.43/224: `1^m + 9` defects, n³/685, 2.90). But the **new appendix justification** says nine is "the largest number of defects that the best known machine can absorb with an O(log n)-bit index, and at ten the machine's cost per test stops being flat." `PC_PUSH.md` §2.2 says the max is **6 plain / 8 gap-coded** and "`k = 9` is outside both encodings"; §3.1 says `C(m)/m` is **flat for every fixed defect count**, k=16 included (0.0012 n³). | **NEW ERROR** — both halves of the new sentence are contradicted by the report it cites | app 350–358 |

Also re-checked and **sourced**: A8's certificate totals. `PC_PUSH.md` l.305 gives n=1024, 8 779
queries, sum 1 285 703, union 1 020, 1260×. The rephrasing in Thm 5 ("On the constant tape at
n=1024 …") introduces no new claim.

---

## 4. The four new checks

**(1) Abstract vs §1 vs §2.1 on the step-count span and exponent — consistent, all three wrong.**
Abstract (38–40), §1 (83–85) and §2.1 (210–212) all print 0.385 → 15.8, "a factor of forty",
1.70–2.20. Consistency passes; correctness does not: the maximum stored value is **15.576**
(hashmod, seed 2026, n=4096, k=8) or **12.541** at seed 777. Nothing stores 15.8. The "forty" is
right for 15.58/0.385 = 40.5. The other two numbers regenerate exactly: 0.385 = mean over all 160
taps (every ending) of NBLK>3; 1.70–2.20 = 1.702–2.203 over the 15-rule el grid, all-runs means, 1024→4096.
Caveat the paper does not state: the ratio is an **all-runs** mean (halts and repeats included),
not a meeting-tape mean; on meeting tapes only, canonical is 3.47 and NBLK>3 is 0.71.

**(2) κ table vs `probe_kappa_inc2/ki2_five.json` at c=16 — matches exactly.** All ten endpoints:
trip 0.789→0.472, peak 0.785→0.614, pattern 0.693→0.502, dense2 0.632→0.509, runs 0.569→0.520.
The "0.999 at n=512 (0.927 at n≥1024), c=128" also regenerates from `ki_agg`.

**(3) Does any table still print an exponent without its window?** `tab:threecosts` does state
"log–log exponents over n=512…4096". **Three prose sites fail:**
- 545–552: "For the band family at **n=4096**, … G=n^{1.09–1.90}, e=n^{1.49–2.15}, g=n^{2.00–2.02}" — an exponent cannot be "at n=4096", and no fitting window is given (it is 512–4096, from `k3_*.json`).
- 363–364: "g=n^{2.00–2.02} for fourteen of nineteen rules, against n^{−0.36–2.17} for the sample maxima" — no window, and two different estimators (witness panel vs sample maxima).
- 327–328: "the measured a is 0.08 to 0.98 depending on how the tail is defined" — no window and no named cut (it is the 10×-median column over 512–4096; the 100 n column gives 0.25–1.57).

**(4) Leftovers the fix created or failed to propagate**
- **"fourteen of nineteen"** at 364 should be **eighteen**: `gw_fit.json` has 18 rules, and 14 of them have g_wit_exp in [1.99, 2.03]. A18's denominator was fixed; A19's was not.
- **"the eleven rules not listed"** (kappa caption, 771): the 13-rule grid has **8** not listed.
- §6 prints the five-rule range twice as **"0.50–0.80"** (778, 802); the table's five rules span **0.472–0.789**.
- The appendix step paragraph still opens "the mean step count at n=4096 divided by n²/128 **is:**" and then has no table — the replaced table left a dangling colon (369–370).
- Appendix status table: "**runs costs 0.12 n², exponent 1.99**" (449) contradicts the body's fully linear runs (3.8–8.5 n). It is the halt-gate number; the gate is not named. Pre-existing, but it is a table claim in a file that was edited.

---

## 5. Still open, priority order

1. **C16 and κ\* ≥ 0.80 (730–733; app 450).** The sentence "thirteen rules between 0.47 and 0.80 and hashmod at approximately 0" is false at every measured size: 9 of 13 in range at n=512, three above 0.80 (canonical 0.925, NBLK>2 0.826, NBLK>3 0.819), one below; at n=2048 BAL>m/2 is negative. Either say "the five re-measured rules span 0.47–0.79" and drop the grid claim, or print the grid and take the bound from its max. The appendix still carries the older, also unsourced κ∈[0.196,0.432] "for fifteen of sixteen rules" and must be deleted, not just the caption.
2. **The 15.8 endpoint (38–40, 83–85, 210–212, app 371).** Replace with 15.6 (seed 2026) or "12.5–15.6 over two seeds"; 15.8 is in no JSON. This is the paper's headline "no clock" number and it is now repeated four times consistently, which makes it harder to dislodge.
3. **D8's scan range (777).** 0.119–0.344 mixes the n=512 minimum with the n=4096 maximum. The five-rule range is 0.119–0.361; at fixed n=4096 it is 0.123–0.344. Pick one window and print it.
4. **D14 (769–770 vs 787, 803).** Two CES values for one device: "3.4–4.1" in the caption, "3.3" twice in the prose. And "R²≥0.99 on runs" is only true at n=1024/4096 (3.40/3.49); at n=512 and 2048 the 0.99 arm costs 58.07 n and 8.45 n. Restrict the claim to pattern, or print the runs cells.
5. **The new nine-defects rationale (app 350–358).** It is the one place where a fix invented a claim: PC_PUSH says the index holds at most 8 defects (9 is outside both encodings) and that C(m)/m stays flat at k=16. Revert to the report's own statement: nine is a convenient fixed count, not a capacity limit.
6. **The ladder's n² row (681).** Still "≤1.12 n²" against the paper's own 1.33 n² bound, and it lists AMAX/AMIN under `band + frontier` without the certifier. Both are one-line fixes.
7. **B12's sentence (604).** "3.8 to 12 with one exception at 3.85" is self-contradictory; write "3.85 to 12, with KCHG below 4".
8. **A19's denominator (364).** "fourteen of nineteen" → "fourteen of eighteen".
9. **M5 / M6.** The abstract still quotes the 1.33/1.12 constant with no device, and the ladder caption still promises devices but lists mechanisms and no probe directory.
10. **Tidy-ups:** "eleven rules not listed" (771) vs 8; "0.50–0.80" (778, 802) vs 0.472–0.789; the dangling colon in `app:steps`; "runs costs 0.12 n²" in the appendix status table.

## 6. Residual risk (unchanged by the fixes)

- The 15.8 endpoint is seed-sensitive: 15.58 (2026) vs 12.54 (777) at n=4096, and hashmod's `k`
  changes across sizes (24 at 512, 16 at 1024–2048, 8 at 4096), so the 1024→4096 exponent mixes two
  protocols. The paper does not state that the "no clock" numbers are all-runs means, nor which seed.
- The five kappa values are a lower bound from one ridge family; the caption now says so, but the
  body still converts "R²≥0.99 at a finite budget" into κ=1 without remarking that κ is a limsup over
  c and n.
- C14's "A" is the report's "D"; if a reader follows the pointer to `CX_DOWN.md` the symbol will not
  match.
