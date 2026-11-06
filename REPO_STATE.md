# Repository state, and the decision about it

*2026-11-06. Written so that a later session does not have to rediscover this.*

## Where things are

| Item | Where |
|---|---|
| Public repository | `https://github.com/peuaza/triangular-window-walk`, tip **`9ed31e1`** (2026-11-03), 35 files: `bthm0/` (8) + paper-1 root scripts + `README.md`/`RUNTIME.md` |
| Local `main` | 95 commits ahead of `origin/main`; **none pushed**, and none should be until the BTHM paper is written |
| Paper 2's code (`family/`, 20 files) | Last present at local commit **`dca6555`** ("Add the six table scripts of paper 2 and map every table to its script"); **absent from the public tip**, because the rollback to `9ed31e1` removed it |
| Paper 2's preprint | Published; `family/README.md` is titled "Code for paper 2: *Memory and randomness control the shape*" and maps every table to its script |

## The decision (user, 2026-11-06)

Paper 2 will be **set to private** for now. The repository will be dealt with **in one pass
after the current paper is written**, not piecemeal.

## What was already checked, so it need not be rechecked

* `9ed31e1` is an **ancestor** of `dca6555`; the rollback therefore dropped `family/`.
* The public `README.md` at `9ed31e1` describes **paper 1 only** ("a conjectured arcsine limit").
* Pushing local `main` would publish all 95 commits, including unpublished work (Conjecture 2,
  `NBLK`, the kappa measurements, the four-valued model). **Do not push `main`.**
* The safe shape for a later fix is a branch built on `9ed31e1` that adds only `family/` (taken
  from `dca6555`) plus a README update, pushed as `restore-family:main`, so that exactly one
  commit appears on the public tip.
* `family/README.md` refers to `_dist/<n>.jsonl`; `_dist/` was never public, so a restore must
  decide whether to include the data.

## The rule that stays in force

Unpublished work does not go to the public repository. "Commit" means `git commit`, which is
local; it never implies a push.
