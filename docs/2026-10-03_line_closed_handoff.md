---
date: 2026-10-03
status: >
  HANDOFF at a clean stopping point. The E39/E40/E41 line is COMPLETE, RUN,
  VERIFIED, WRITTEN UP and PUSHED. NO PHASE IS OPEN and nothing is half-done.
  `main` is at 652b842 and level with origin; the tree is clean. This note is
  about what to do NEXT, not about what was done - the results live in
  docs/2026-10-03_e39_e40_e41_result.md and are not repeated here.
posture: >
  EXPLORATORY and POST-HOC throughout. Everything the line produced is
  hypothesis-generating. The ONE clause it licenses carries three
  qualifications that are not optional, and the list of things it does NOT
  license is longer than the list of things it does.
branch: main at 652b842, pushed, level with origin/main, tree clean.
supersedes:
  - docs/2026-10-03_e41_handoff.md (sections 3-6 of it are still the E41
    record; only its section 7 is out of date, and it is bannered)
relates-to:
  - docs/2026-10-03_e39_e40_e41_result.md (THE RESULT. Read it before anything)
  - docs/2026-10-02_E39_respiratory_axis_decomposition_declaration.md (governing)
  - docs/2026-10-02_e39_e40_handoff.md (the E39/E40 leg; its open integrity item
    is closed in section 8 of the result note)
  - docs/2026-10-02_handoff.md (the E36/E37/E38 line; its "no phase is open" was
    briefly wrong and is TRUE AGAIN)
next-action: >
  Nothing in this repo is blocked. The one live obligation is on the
  MANUSCRIPT, not on the code: section 2 below.
---

# Handoff - the respiratory-axis line is closed, and the live item is on the manuscript

> **SPENT. SUPERSEDED by `docs/2026-10-04_handoff.md`.** Its section 2
> obligation - replacing the draft sentence *"respiration did not itself predict
> response"* - **STANDS**, but now carries a fourth qualification from `E42`:
> the subtype scope was checked and the clause stands unqualified, with
> `HRneg_HER2pos` never fitted. **Sections 3 to 7 here are still good as
> written.** Two things postdate this note entirely: `E42` and `X01`.

**Option A throughout: Claude Code wrote the scripts, the author sourced them.**
Nothing was written to `myc_human_validation`, frozen at `d3ac60e`, or to
`myc_mouse`.

---

## 1. State, in four lines

- **No phase is open.** The respiratory-axis line finished and was written up.
- **`main` is at `652b842`, pushed and level with `origin/main`. Tree clean.**
- **The next free script number is `E42`.** `E31` is still the only analysis off
  `main`, on `bh3-mimetic-oxphos`, deliberately.
- **`CLAUDE.md` is current.** It carries the line, its three travelling
  qualifications and the two new traps.

## 2. THE ONE LIVE OBLIGATION, and it is not code

**The draft contains the sentence *"respiration did not itself predict
response"*. That sentence must be replaced or dropped.** It was the reason the
whole line was run, and the finding is that **its supporting quantity had never
been fitted** - only three numbers printed once in an `if (FALSE)` sandbox at
`scripts/13_outcome_models.R:745`.

Now fitted, that quantity is `m0`, and it **fails the declared reading rule in
both directions**: signs `+ - +` across the three cohorts, pooled CI including
zero, I2 = 48.8, and the only cohort-level `m0` excluding zero is the smallest
cohort disagreeing in sign with the largest - the exact Block F1 shape section 9
cites as a warning. **So the sentence cannot be re-sourced and cannot simply be
kept.**

**What may replace it is one clause**, and the three qualifications travel with
it (result note 7.1):

> Once proliferation is adjusted for, higher OXPHOS subunit expression is
> associated with lower pathological complete response across three neoadjuvant
> cohorts (pooled -0.16 log odds per within-cohort SD, 95% CI -0.28 to -0.04,
> sign consistent in all three, I2 = 0).

1. **SINGLE-INSTRUMENT.** mitoPPS exists in none of the four cohorts.
2. **The unadjusted coefficient does NOT support it** and is reported beside it.
3. **No coefficient compared numerically with a TCGA or SCAN-B one** - GSE25066
   covers 57 of 89 OXPHOS subunits on GPL96.

**Two further licensed statements** are in result note 7.1: the forkscale
INTERMEDIATE statement, and that the identifiability limit was measured rather
than assumed.

## 3. What a next session MUST NOT write

These are the ones most likely to be reached for, and each is forbidden by a
declared rule rather than by taste. **Result note 7.2 is the full list.**

- **That METABRIC replicated, confirmed or extended the pCR finding** - 10.5,
  whatever the sign. Different estimand, and the two spines differ so the two
  `m0` are not the same model.
- **Any distant-outcome claim.** Nothing meets the rule in any pre-specified
  set, in either cohort family.
- **The `s508` DRFS cells** (HR 0.776 / 0.774 / 0.802, CIs excluding zero).
  **That set has NO SUBTYPE TERM.** The pre-specified 465 and the 490 have
  nothing. The coefficient crosses the line exactly when subtype adjustment is
  removed.
- **METABRIC's ER-positive `m1`** (HR 0.8962, p 0.0492). It excludes zero by
  3.7e-04, **its own declared sensitivity does not reproduce it** (p 0.186), and
  `m2` and `m3` include zero. One boundary cell of thirty.
- **`m1 -> m2` as "what remains after proliferation is MYC."** Declared out, and
  E39's `rho(MYC, PROLIF_DISJOINT) = 0.8236` in GSE164458 - inside the
  declaration's own NOT IDENTIFIABLE band - makes it mandatory.
- **`m3` as "buffering capacity."** It is adjustment for two anti-apoptotic
  transcripts whose limbs correlate at 0.019, 0.120 and -0.016.
- **Two-instrument agreement**, a causal claim, an `ER x OX` interaction (never
  fitted), a per-gene FDR, "primed" language, or any p value called a near miss.

## 4. If anything from this line is to be pushed further

**Four falsifier lists are already written**, in result note section 10, as
`CLAUDE.md` requires - before the next analysis, not after. Nothing needs
inventing:

- **10.1** four falsifiers for the licensed pCR clause. The cheapest is the
  second: swap `PROLIF_DISJOINT` for another of the arm's proliferation scores
  and see whether the pooled `m1` CI crosses zero. The clause asserts a
  proliferation-adjusted association, not one adjusted for a particular
  318-gene list.
- **10.2** the three conditions under which METABRIC's ER-positive cell would
  become a finding. All three are currently unmet.
- **10.3** two falsifiers for the forkscale statement.
- **10.4** records deliberately that the two unanticipated forkscale readings -
  `rho(PROLIF, MB1) = 0.619` exceeding `rho(OX, MB1) = 0.416`, and
  `rho(MB1, MB2) = 0.865` - have **no** falsifier list **because they are not
  claims**. They invite a reading about what the published fork axis actually
  measures. **That reading is not licensed and must be declared first.**

## 5. Two traps this line added, and they generalise beyond it

Both are in `CLAUDE.md` and in declaration 10.3a, and both will bite the next
person who harmonises symbols or touches the configuration in METABRIC.

1. **An alias is not automatically safe, and the guard belongs on the ROUTE, not
   the symbol.** `COX7A2L -> SCAF1` and `POLR1B -> RPA2` each point at the
   HGNC-approved symbol of a **different** gene (58506, 6118). Taking the first
   would have put a splicing factor inside the OXPHOS exposure; both are
   rejected, which is why OXPHOS subunits recover 69 of 89 and not 70. **But a
   blanket ban on those symbols is also wrong**: `RPA2` is a legitimate *direct*
   member of `PROLIF_DISJOINT`, and banning the symbol would have silently
   altered the proliferation score in every model. **The E41 script found this by
   stopping on it in a dry run.**
2. **`BBC3` is absent from METABRIC.** Costless here, because the configuration
   is not fitted - but `BBC3` is the trigger limb the guardian-switch finding is
   built around, so **any future METABRIC analysis touching the 6-gene signed
   configuration is blocked on it**, and a five-gene stand-in is forbidden.

## 6. Housekeeping, and one thing deliberately left alone

1. **`results/` and `outputs/` are gitignored and do not survive a clone.**
   Re-source `E39`, `E40`, `E41` in that order - all three deterministic, about
   five minutes in all, and `E41` reproduces the 26 digests in its handoff
   section 3. **The result note is the durable record.**
2. **The declaration's `next-action` frontmatter is stale** - it still reads
   "write the script. Nothing is fitted until this note is committed."
   **I have deliberately NOT corrected it.** The declaration is the
   pre-registration for this line; editing its frontmatter after the fact, even
   to fix something plainly false, is the one kind of edit that a declare-then-run
   repo should not make casually. **If it is to be fixed, fix it as a dated
   amendment that says what changed and when**, as `ccadfa5` and `82659ff` did
   for its body.
3. **`E41` has no figure**, as `E40` had none. The briefs asked for tables only,
   so `CLAUDE.md`'s "from `E31` on, write the key figure twice" has nothing to
   apply to. **If a figure is wanted, the obvious one is the four-rung ladder
   with its CI, two endpoints, both cohort families** - and it would need a
   decision about whether the non-licensed cells appear on it.
4. **`A7`'s Figure 7 and Figure 3C claims remain AUTHOR-SOURCED and
   unverified.** The Menegollo PDF is in no repository, so no script here can
   check them. The author undertook to re-read the panel labels **before that
   section is cited**.
5. **MB3 is unused and its identity still unverified** (10.3): whether
   `MB3.forkscale == MB3.pc1.rev / MB3.index`. Verify before any use.
6. **METABRIC's provenance rests on the author's authority** - Synapse
   `syn1757063` / `syn1688369`; nothing inside the file names Synapse, the readr
   `col_spec` is consistent-with and not evidence, and the md5 is the only
   check. Recorded in `data/menegollo_biclusters/README.md` and declaration 4.5.
7. **Two items still open on purpose** from `docs/2026-09-11_handoff.md`:
   whether standalone mouse-only figure reproductions are wanted (a `myc_mouse`
   task, not buildable from here), and whether `data/from_myc_mouse/` should be
   tracked rather than gitignored.

## 7. The commit trail for this leg

| commit | what |
|---|---|
| `ccadfa5` | declaration amended for METABRIC (4.5, 10.1, 10.1a, 10.2, 10.3a, 10.4, 10.5, 10.6) and the Menegollo README provenance corrected, in one commit |
| `82659ff` | 10.3a amended: the alias guard is on the route, not the symbol - forced by the script stopping on it |
| `750af99` | `scripts/E41_metabric_ladder.R`, written and not yet run, carrying its 26 pre-run digests |
| `a364e36` | display-only: CI to 4 dp and `excl_0` computed, after the run exposed a rounding hazard on the one boundary cell |
| `954664e` | the E41 handoff, now bannered as superseded |
| `09331cc` | **the result note for all three scripts** |
| `652b842` | `CLAUDE.md` brought up to date; next free number `E42` |

**Note the order**: declaration amended and committed **before** the script was
written; script committed **before** it was run, with its digests in the
message; result note **after** the run was verified against them. That sequence
is the audit trail and it survives in the history.

## 8. The shortest possible summary

**One quantity met the declared reading rule** - the pooled pCR coefficient at
`m1`, `-0.1616 [-0.281, -0.042]`, sign consistent in three cohorts, I2 = 0.
**The unconditional coefficient did not**, which is why the draft sentence has
to go. **Distant outcome produced nothing that met the rule**, in 1,981
METABRIC patients with 623 cause-specific events or in 465 GSE25066 patients
with 103 - and both of its near-misses are traps with named causes rather than
results. **The forkscale question got a third reading and INTERMEDIATE turned
out to be a property rather than an instrument artefact.**

Everything else this line produced is a limitation, a trap or a falsifier, and
all of them are written down.
