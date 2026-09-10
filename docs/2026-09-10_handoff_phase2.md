---
date: 2026-09-10
status: SPENT 2026-09-11. E24 is written, sourced and reported; the result note
        docs/2026-09-10_e24_mitopps_rank_comparison.md is the entry point for
        Phase 2 from here. This file is kept only as the record of the state
        E24 started from, and is deletable.
state: main clean at 5dcd0b6, three commits unpushed. Phase 1 CLOSED. Phase 2
       PLANNED, NOT STARTED - GATE 0 passed, no number below it computed, no
       script written.
plan: /Users/gs/.claude/plans/witty-snuggling-dewdrop.md - the E24 design, held
      outside the repo because the session that wrote it never left Plan Mode.
supersedes: nothing. docs/2026-09-08_handoff.md remains the entry point for
       everything OUTSIDE Phase 2, and section 4 of docs/2026-09-04_handoff.md is
       still the live consolidation list.
---

# Handoff - 2026-09-10: Phase 2 is planned, not started

> **SPENT 2026-09-11.** `E24` was written, dry-run, sourced by the author and
> written up. **G1 passed, G2 FAILED (the rank is size-structural, so R1 was
> scored on the pre-specified size-matched band), G3 falls but its interval
> covers zero, and R2's cohorts disagree.** Everything below describes the state
> before that happened. Read
> `docs/2026-09-10_e24_mitopps_rank_comparison.md` instead.

**Nothing was analysed, nothing was computed, nothing was run.** The session that
produced this ended in Plan Mode without writing a script, a result or a figure.
The only thing that exists is GATE 0's answer and a design.

> **The plan itself is NOT in this file.** It lives at
> `/Users/gs/.claude/plans/witty-snuggling-dewdrop.md` and is the plan of record
> for `E24`. This handoff is the state; that file is the design. Read this first,
> then the plan.

## State

| | |
|---|---|
| repo | `myc_human_exploratory`, `main`, **`5dcd0b6`** |
| tree | **clean** |
| unpushed | **3 commits** — `9b8d349`, `852731b`, `5dcd0b6`. Nothing has been pushed since Phase 1 opened |
| Phase 1 | **closed.** Corrections and input audit done; do not revisit |
| mouse artefact | snapshotted at `data/from_myc_mouse/`, md5 `b8a125af4bc0d5909ad03f6e126e2890`, provenance README tracked, artefact gitignored |
| Phase 2 | **planned, not started.** No `E24`, no results, no figure, no note |

## What happened this session

**GATE 0 was run and it PASSED.** The full table is in the plan file. The
headline: the two panels are the same object — same catalogue, the **identical
eight pathways** drop under the filter in both species, the `mt-` strip is
identical and was verified **against the saved object** rather than any script,
and `OXPHOS subunits` is **89 in both**. Intersection is **142 pathways**.

**One question was asked and answered.** The mouse panel carries two extra
pathways, `Apoptosis-PRO` (25) and `Apoptosis-ANTI` (9), a disjoint partition of
mouse `Apoptosis` (34). You asked for the same split on the human side.

**That answer was researched and it changes the build.** Three findings, all in
the plan:

1. **It costs a human mitoPPS recompute.** `functions/mitopps.R` puts the
   universe size `P` in the denominator of both `.mitopps_universe` (`P - 1`) and
   `.mitopps_query` (`P`), so adding two pathways changes **every** human score.
   The saved 142-row `$mitopps_universe` cannot be extended.
2. **It is buildable human-natively** from
   `data/genesets_celldeath_human/cell_death_genes_consolidated.csv`, which has
   first-class human `effect` calls. **No ortholog function, no mouse-derived
   tree.**
3. **It is not the mouse's set.** PRO agrees exactly at **25**; **ANTI does not —
   9 mouse against 6 human**, and the human set sits just above the `min_genes = 3`
   floor. Coverage is 34/34 mouse against 31/35 human. An independently curated
   analogue, not a translation, so any Apoptosis rank read across species is
   **descriptive only** and no rule is pre-specified for it.

The plan puts this in as **PART A2, a clearly bounded third leg**, with its own
control: rebuild the human 142 universe and assert it reproduces the saved object
bit-equal **before** any 144 number is read.

## The one thing pending

**The plan was presented for approval and you declined to proceed, to exit.** It
is unmodified below. Next session: read it, approve or redirect, then `E24` gets
written. **No number below GATE 0 has been computed.**

## Do this first, next session

1. **Decide the plan** — approve, or redirect PART A2, which is the only part
   that grew beyond the original brief. It is at
   `/Users/gs/.claude/plans/witty-snuggling-dewdrop.md`.
2. **Then `E24`.** Option A: I write it, you source it in Positron. The note is
   written only after your run reproduces the dry run object by object.

## Two things noticed and deliberately not fixed

- **`docs/2026-09-10_mouse_corrections.md` §7.4 is now stale.** It says "No mouse
  artefact is snapshotted here", which was true when the audit ran and is false
  since `5dcd0b6`. Phase 1 is closed and I was told not to revisit it, so it
  stands. **One line to fix whenever the file is next open.**
- **The artefact is gitignored, following `data/from_validation/`.** That sibling
  is ignored because it is ~560 MB; this file is **511 KB** and could be tracked,
  which would make it drift-proof rather than drift-detectable. Convention was
  followed rather than re-derived. Flagged in
  `data/from_myc_mouse/README.md`; a one-line `.gitignore` change reverses it.

## Constraints that bind next session

Rank never value, no cross-species value in any output. No pooling of mouse and
human into one scoring run. TCGA and SCAN-B separate, never averaged. The mouse
artefact is read-only input — **do not add the mouse repo to the session and do
not recompute mouse mitoPPS**. `myc_human_validation` untouched. No FDR across
pathways. N3. **G3 is genuinely open**: if OXPHOS's percentile does not fall from
6W WT to 12W WT, the analysis stops there.

---
