# ORDERS — the order of `p`, `e`, `g` per rule

*2026-11-09, 20:05. Assembled from `probe_table_eg` (five-size grid, meeting-tape domain),
`probe_table_p` (the `p` battery), `probe_ord_band` (K-matched fits, 14 rules),
`probe_ord_prefix` (the `p/e` slopes), `probe_ord_constructed` (constructed rules, audited).
Every number is script output; the evidence class is stated for each entry.*

## 0. Definitions (fixed)

Domain: the meeting tapes of a rule, as `n → ∞`. Each measure takes its own best machine among
those with `O(log n)` bits of state:

| measure | requirement | quantifier over tapes |
|---|---|---|
| `g` | exact on every meeting tape | worst |
| `e` | exact on every meeting tape | mean |
| `p` | only `R² → 1` | the budget at which `R²` reaches 1 |

`p ≤ e ≤ g` is a consequence of the definitions and is used as a per-row consistency check.

## 1. The table

Evidence classes: **construction** (exact by design), **theorem** (proved), **decomposition**
(mechanism), **witness** (an explicit tape or predictor), **fit** (K-matched or five-point OLS with
a t-interval), **bound** (a lower bound on the true order).

**Read the `g` column as a lower bound only.** It is the maximum over a uniform sample of 300–6000
tapes; the structured tapes that maximise the walk are essentially never drawn. The witness panel
(`probe_g_witness`, nine structured tapes built as balanced column pairs) lifts `g` by factors of
1.3–274 at `n = 512` for 8 of 19 rules, and the lifted values have clean exponents
(`g_wit = n^{2.00…2.02}` for 14 of 19 rules, against `n^{-0.36…2.17}` for the sample maxima; the
sample maximum has a *negative* exponent for `BAL>m/2` because at `n = 4096` only two or three
meeting tapes remain per cell). Use the
witness exponent where a witness exists; the sample exponent is noise.

**The `e` column is also systematically high for the superlinear rules.** `probe_en_linear` found an
exact-on-every-tape machine at `n + steps + cuts` reads (dropping `pcr`'s re-span). Measured on the
same tapes (`probe_e_cheap`, 200 meeting tapes per cell, 2026-11-10):

| rule | `e_old/n` @2048 | `e_new/n` @2048 | ratio | legality of `n+steps+cuts` |
|---|---|---|---|---|
| `dense2` | 443.2 | 120.8 | 3.67 | verified (`EN_LINEAR` §8) |
| `pal` | 350.4 | 112.8 | 3.11 | sliding statistic, assumed |
| `LAS` | 256.5 | 66.2 | 3.87 | assumed |
| `oddrun` | 237.4 | 79.0 | 3.01 | verified |
| `canonical` | 236.7 | 65.4 | 3.62 | assumed |
| `trip` | 223.5 | 65.4 | 3.42 | verified |
| `peak` | 212.6 | 47.9 | 4.43 | verified |
| `hashmod` | 193.2 | 193.0 | 1.00 | no re-span to drop |
| `KCHG` | 102.9 | 26.2 | 3.93 | assumed |
| `BAL>m/2` | 60.9 | 15.6 | 3.91 | **invalid**: its anchor is not a carried statistic --- it moves with the cut, so the re-span is necessary (EN_LINEAR §5), and the record's `n²/4` stands |
| `BAL>m/4` | 55.6 | 12.7 | 4.37 | **invalid** (same) |
| `NBLK=1` | 61.6 | 18.4 | 3.34 | **invalid** (same) |
| `NBLK>2` | 13.8 | 4.0 | 3.48 | **invalid** (same) |
| `NBLK>3` | 13.4 | 3.9 | 3.46 | **invalid** (same) |
| `MAXP` | 224.3 | 74.7 | 3.00 | **illegal**: a window extremum needs a scan per step, so the row is a step count |
| `AMAX` | 236.8 | 50.3 | 4.71 | **illegal** (same) |
| `AMIN` | 169.4 | 37.0 | 4.57 | **illegal** (same) |

The class assignments do not change: every ratio is a constant, and every `e_new` exponent over
`n = 512/1024/2048` lands in `1.38`–`2.29`, the same mixture band as before, whose asymptote is
`n²` (the tail's weight rises; see the band-rule note above). What changes is that the paper must not
quote the record's `e` values as minima: on the eight rules where the machine is verified or assumed
legal they are `3.0`–`4.4` times too high, and on the five level rules the record's `n²/4` stands
because their anchor is not a carried statistic.

### Natural rules

| rule | `p` | `e` | `g` | evidence | caveat |
|---|---|---|---|---|---|
| `monotone` | `n` | `n` | `n` | construction (`e = g = 2n−1`) | — |
| `longrun` | `n` | `n log n` | **`n³/21`** | `e`: decomposition (`4.87 + 0.45 ln n` per cell); `g`: **witness** `2^{n/2}0^{n/2}` measured at four sizes: `6 544 735 / 50 408 550 / 395 496 223 / 3 132 928 486` = `n³/20.5, 21.3, 21.7, 21.9` | the "`25 n²`" of the record is the `n = 512` slice of this; the coefficient of `n²` doubles per doubling of `n`, so the order is `n³`; `g/e = 86 232` at `n = 4096`, growing like `n²/log n` |
| `pattern` | `3.6 n` | `3.7–4.2 n` | `Θ(n)` | **CES exact on 100 % of the meeting tapes** at `n = 512/1024/2048/4096` (200 tapes each); `p ≤ e` in every cell | `p` was ~`1 n` too high: `kp_lib._Win.__init__` charges the opening pass twice. A fully exact machine (halts included) does not exist: the terminal `M0` reaches `0.99 n` |
| `runs` | `3.4 n` | `3.8–8.5 n` | `Θ(n)` | CES + the `W2` endgame, `100 %` of the meeting tapes at `n ≥ 1024` | the `n = 512` cell is `58 n`, forced by a terminal at `M0 = 439`; measured, not proved |
| `dense2` | `0.4 n²` | `n²` | `n²` | `p`: witness counter, `p/e` flat 0.44 → 0.41 | a near-machine reaches `2.2 n` but is exact on only `574/700` meeting tapes; the `b2=1` epoch law fails mid-tape (83 343 misses, 2 with epoch length ≤ 4) |
| `peak` | `0.3 n²`? | **`2.0 n` (99.8 % exact) / `≈7 n` (exact)** | `Θ(n)` | `en_linear`: `e/n = 2.0` flat over four sizes; `probe_peak_fix`: `1 496/1 500` at `W2 = 128`, `1 497/1 500` at `W2 = 256` | the "all three failures are first-epoch" of `EN_LINEAR` is **wrong**: two of the four are mid-tape (`M0 = 576`, `1 529`), so making epoch 0 literal is a no-op (exactness unchanged, cost `3.0 n`). Exactness needs `W2 = 0.75 n` and costs `≈7 n` (20-tape sample at `n = 2048`) |
| `trip` | `0.5 n²` | `n²` | `n²` | `p/e` flat 0.57 → 0.51 | near-machine `2.1 n`, exact `690/700`; the `b2=1` law fails at 2 409 epochs |
| `hashmod` | `n²` | `n²` | `n²` | `p = e` on every arm; `g = 0.35 n²` | the record's `0.62 n ln n` is a **row count**, not cell reads |
| `canonical` | `n^~2` | `n^~2` | `n²` | fit | near-machine `3.0 n`, exact `696/700`; `p/e = 0.950` at n=1024 is a grid jump |
| `pal` | `n^~2` | `n²` | `n²` | fit (1.94, 2.08) | — |
| `MAXP` | `Θ(e)` | `n^{1.83}` | `n^{1.89}` | fit, K = 153 | near-machine `47 n` at 4096, exact `696/700`; `b1=1` bucket 0.506 on halting vs 1.000 on meeting tapes |
| `LAS` | `Θ(e)` | `n^{1.94}` | `n^{2.05}` | fit, K = 72 | near-machine `2.1 n`, exact `689/700` |
| `oddrun` | `Θ(e)` | `n^{2.09}` | `n^{1.90}` | fit, K = 108 | near-machine `2.1 n`, exact `690/700` |
| `AMAX` | `Θ(e)` | `n^{2.15}` | `n^{2.21}` | fit, K = 82 | 0 wrong cuts in 700 tapes, but the meeting terminal reaches `M0 = 820` at `n = 2048` (0.4 n) ⇒ `W2 = Θ(n)` ⇒ endgame `Θ(n²)` |
| `AMIN` | `Θ(e)` | `n^{1.83}` | `n^{1.95}` | fit, K = 93 | same terminal-tail failure (`M0 = 820`) |
| `KCHG` | `Θ(e)` | `n^{1.29}` | `n^{1.83}` | fit, K = 160 | near-machine `2.0 n`, exact `694/700`; `g/e` grows 12.6 → 50.7 |
| `BAL>m/2` | `Θ(e)` | `n^{1.57}` | `n^{1.92}` | fit, K = 141 | near-machine `1099 n` at 4096, exact `678/700`; `g/e` 13 → 29 |
| `BAL>m/4` | `Θ(e)` | `n^{1.61}` | `n^{1.96}` | fit, K = 149 | **exact 700/700** but the level query's anchor moves with the cut ⇒ `n²/4` (to within 6 %) |
| `NBLK=1` | `Θ(e)` | `n^{1.79}` | `n^{1.81}` | fit, K = 60 | **exact 700/700**, `n²/4`: the anchor slides (`0` failures above the buffer in 700 tapes, but `Θ(M)` per epoch) |
| `NBLK>2` | `Θ(e)` | `n^{1.75}` | `n^{2.03}` | fit, K = 200 | **exact 700/700**, `n²/4`, same anchor slide |
| `NBLK>3` | `Θ(e)` | `n^{1.69}` | `n^{2.01}` | fit, K = 148 | **exact 700/700**, `n²/4`, same anchor slide |

### The witness panel (`probe_g_witness`)

Nine structured tapes, each built as a pair of balanced columns so that it lies in the ensemble's
support. `g_wit` is the maximum over them; the exponent columns are a four-point OLS over
`n = 512/1024/2048/4096`. Measured 2026-11-10.

| rule | `g_wit/n²` at 512 | `g_wit` exp | `g_K` exp | `g_wit/e` at 4096 |
|---|---|---|---|---|
| `AMAX` | 1.08 | 2.02 | 2.11 | 11.5 |
| `AMIN` | 1.08 | 2.02 | 2.07 | 12.2 |
| `BAL>m/2` | 1.26 | 2.00 | **−0.36** | 54.9 |
| `BAL>m/4` | 1.26 | 2.00 | 0.13 | 42.5 |
| `KCHG` | 0.76 | 2.00 | 0.68 | 20.1 |
| `LAS` | 0.76 | 2.00 | 2.10 | 8.3 |
| `MAXP` | 0.76 | 2.00 | 1.82 | 7.4 |
| `NBLK=1` | 0.85 | 2.01 | 1.77 | 14.4 |
| `canonical` | 0.01 | 1.00 | 2.17 | 0.0 |
| `dense2` | 1.18 | 2.00 | 2.14 | 4.8 |
| `hashmod` | 0.24 | 2.16 | 2.00 | 3.8 |
| `longrun` | 24.97 | 2.97 | 0.78 | **86 232** |
| `monotone` | 0.004 | 1.00 | 1.00 | 1.0 |
| `oddrun` | 0.76 | 2.00 | 1.78 | 5.9 |
| `pattern` | 0.76 | 2.00 | 2.01 | 2.9 |
| `peak` | 0.76 | 2.00 | 1.67 | 6.6 |
| `runs` | 0.63 | 2.00 | 2.02 | 4.8 |
| `trip` | 1.25 | 2.00 | 1.48 | 13.2 |

Three readings:

- **`g = n²` for fourteen of nineteen rules** by the witness, against exponents of `−0.36`–`2.17` from
  the sample maxima. The witness is the estimator; the sample maximum is noise (and at `n = 4096` it
  can even fall with `n`, because only two or three meeting tapes remain per cell).
- **`longrun` is the one genuine separation**: `g/e = 86 232` at `n = 4096`, growing like
  `n² / log n`, because its witness is `2^{n/2}0^{n/2}` at `n³/21` while its mean is `n log n`. The
  record's `25 n²` is the `n = 512` slice of that sequence.
- **`canonical` is the one rule the panel misses** (`0.01` of its sample maximum): its worst tape is
  not among the nine --- it is generic (`311` runs of length `1`–`8`, counts `132/248/132` at
  `n = 512`, `gw_canonical.py`). So the panel supplements sampling for some rules, not for all.

### Constructed rules

| rule | `p` | `e` | `g` | evidence | caveat |
|---|---|---|---|---|---|
| `nhalf` | `1` read | `n^{0.48}` | `n^{0.48}` | measurement (0.484) + `e = g` | the claim is 0.5 |
| `PC` | `0.999 e` | `n^{2.44}` (Brent) / `n^{2.90}` (scan) | `n^{2.5}` (Brent) | fit ±0.04 | the rung is machine-dependent |
| `DC` | `1.000 e` | `n^{2.86}` ±0.40 | `n^{2.86}` ±0.40 | fit | `r=1` regime; asymptote `n³/√log n` |
| `PRF` | `1.000 e` | `n^{2.94}` | `n^{2.94}` | fit, rising; log excluded at 10.8σ | conditional on the PRF heuristic |
| `CX` | `1.000 e` | `n^{2.63}` ±0.14 | `n^{2.79}` ±0.17 | fit | `g` within 1.1σ of `n³/log n`; `e` 5.5σ below |
| `KRAND` | `1.000 e` | `n^{2.64}` ±0.13 | `n^{2.79}` ±0.17 | as `CX` | — |
| `GRID` | `0.999 e` | `n^{2.90}` ±0.02 | `n^{2.90}` ±0.02 | fit | `n³/log n` |
| `GRID_D2` | `0.999 e` | `n^{3.06}` ±0.04 | `n^{3.06}` ±0.04 | fit | `n³` |
| `GRID_sqrt` | `1.000 e` | `n^{2.55}` ±0.03 | `n^{2.55}` ±0.03 | fit | `n^{2.5}` |
| `GRID_lin` | `0.999 e` | `n^{1.83}` ±0.02 | `n^{1.83}` ±0.02 | fit | `n²/log n` |

## 2. What is solid, and what is not

**Solid.** `monotone` (`e = g = 2n−1`); `longrun`'s `e = n log n` and its `n²` witness; `pattern` and
`runs` with `p = O(n)` by witness; `nhalf` at `n^{1/2}`; `PRF` at `n³`; `CX`'s `g` at `n³/log n` to
1.1σ; the four `GRID` variants as the `D`-dial ladder.

**Not solid.** Every entry marked *fit* is an effective exponent over `256…4096` with a t-interval
that is often ±0.2–0.7; the K-matched maxima are `K/(K+1)` quantiles, not suprema; and `p` is a
minimum over a *declared battery*, so every `p` is an upper bound.

## 3. What the table already says about the paper

1. **The GRID family is not `n²`.** Its four variants sit on four different rungs (2.899, 3.059,
   2.547, 1.831): the `D` dial is the ladder, not a constant.
2. **`PC`'s rung depends on the machine**: the scan is `n^{2.9}`, Brent is `n^{2.5}`. A table that
   quotes one number must name the machine.
3. **`CX`'s claim is about the worst case.** Its own `e` exponent (2.63) is 0.22 below `n³/log n`,
   outside the interval; the `g` exponent is inside.
4. **`longrun`'s `g` cannot be read off a sample.** The K-matched exponent is 0.82; the witness tape
   gives `n²`. Any rule with a named witness must have its order taken from the witness.
5. **Four rules have `e` exponents near 1.3–1.6** (`KCHG` 1.29, `BAL>m/2` 1.57, `BAL>m/4` 1.61,
   `NBLK>3` 1.69). If those survive a larger K they are the paper's first sub-quadratic means.

## 4. Gaps to close next

- `p` for the constructed rules: measured as `Θ(e)` but without an exponent of its own.
- `e` for `DC`, `PRF`, `GRID*`: not fitted in this batch.
- Larger `K` (300+) for the band family, so the quantile exponent moves toward the supremum.
- `peak`: CES arm, to decide whether its `p/e` fall is a power or a slow constant.
