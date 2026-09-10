# E24 - where OXPHOS sits in the mitochondrial panel, in mouse gland and in human tumours

**Run 2026-09-10 23:42; note written 2026-09-11. EXPLORATORY, POST-HOC,
DESCRIPTIVE. Nothing here is pre-registered and nothing here is a hypothesis
test.** Both R1 outcomes were
written into the script header before any number was computed; neither was
predicted.

Script: `scripts/E24_mitopps_rank_comparison.R`. **Sourced by the author
2026-09-10**; `results/mitopps_rank_comparison.rds`,
`outputs/tables/E24_ranks_by_group.csv` and
`outputs/mitopps_rank/E24_oxphos_percentile.pdf` are on disk. The numbers below
were first written from a dry run in the scratchpad. The author's run reproduces
**all 21 saved elements identically** - `identical() == TRUE` element by
element, the seeded 5,000-draw bootstrap intervals included - so every number
here is read from the saved object and **the dry-run caveat is discharged**.

**RANK, NEVER VALUE.** No mitoPPS value is compared across species anywhere in
this analysis. The saved object carries percentile rank, rank correlations and
rank differences only; it has no column named for a score or a value, and the
figure carries none. Mouse and human were never scored together, and TCGA and
SCAN-B are never averaged.

---

## 0. The answer in six lines

| | |
|---|---|
| **Do the two panels support a rank comparison at all?** | **YES.** G1 passed on the strongest available evidence: the **identical eight** pathways drop under `min_genes = 3` in both species, `OXPHOS subunits` is 89 in both, and `mt-`/`MT-` genes form one 13-gene synthetic pathway in both with none inside the primary set or the umbrella. |
| **Is the primary's rank readable against the panel?** | **NO. G2 FAILED.** `rho(rank, log size)` is **0.54** in TCGA and **0.57** in SCAN-B across the intersected panel, and the primary sits at the **94th percentile of the panel by size**. R1 is therefore scored on a size-matched band, which is the remedy the rule prescribed in advance. |
| **Does the mouse reference exist?** | **Weakly.** The OXPHOS percentile falls from 6W WT to 12W WT in every panel definition, but **every interval covers zero** - `-37.3 [-78.2, +3.9]` on the intersected panel, `-47.7 [-86.4, +11.4]` size-matched. n = 6 per arm. The two reference positions are not well separated. |
| **Where do human tumours sit?** | **Not at the 12W position, in either cohort.** On the size-matched panel: mouse 6W WT **54.5**, mouse 12W WT **6.8**, TCGA **38.6**, SCAN-B **47.7**. |
| **Do the cohorts agree?** | **NO, and R2 pre-declared that as the result.** TCGA reads as **between the two mouse positions and near neither**; SCAN-B reads as **near the 6W position**. |
| **What does the apoptosis split add?** | **Almost nothing, and that is worth recording.** The mouse `Myc+` arms collapse on **both** halves at once (PRO 10.4 / 12.2, ANTI 5.2 / 9.0), so the split does not separate a pro-death from a pro-survival move. It reports a genotype effect on the whole apoptosis pathway's rank. |

**The one-sentence reading: neither human cohort sits at the mouse's lowered
12W position, but the mouse reference is too imprecise and the primary set too
large for that to be more than one descriptive observation.**

---

## 1. G1 - the panels are the same object

GATE 0 was answered against the **saved mouse artefact**, never against the
producing script, per `docs/2026-09-10_mouse_corrections.md` section 7.2. PART 0
re-asserts it in code and saves it as `$panel_reconciliation`, so the comparison
carries its own provenance.

| | human (both cohorts) | mouse |
|---|---|---|
| pathways before filtering | 150 | 152 |
| pathways after filtering | **142** | **144** |
| the filter | `min_genes = 3L` | **the identical 8 drop** |
| `mt-` / `MT-` handling | one synthetic pathway, n = 13 | one synthetic pathway, n = 13 |
| `mt-` inside `OXPHOS subunits` | none | none |
| `OXPHOS subunits` | **89** | **89** |
| `OXPHOS` umbrella | 156 | 155 |
| `OXPHOS assembly factors` | 68 | 67 |
| input scale | linear DESeq2-normalised | linear DESeq2-normalised |

The eight dropped in both: `Cholesterol-associated`, `Cytochrome C`,
`Glycerol phosphate shuttle`, `mtDNA modifications`, `OXA`,
`Vitamin B1 metabolism`, `Vitamin B6 metabolism`, `Vitamin C metabolism`.
**That the filter lands on an identical set is stronger evidence of a shared
rule than the threshold value is.**

The two residual ontology differences sit **outside** the primary set. The rule
fixed in advance says they do not downgrade the placements to descriptive, and
they did not.

**Intersection: 142 shared. Human-only: none. Mouse-only: `Apoptosis-ANTI` and
`Apoptosis-PRO`** - the two halves of a disjoint partition of the mouse's
`Apoptosis` pathway, redundant panel members that the intersection drops
automatically.

---

## 2. PART A2 - the human apoptosis split, and what it cost

The split was added at the author's request. It required **recomputing human
mitoPPS from scratch at 144 pathways in both cohorts**, because
`functions/mitopps.R` carries the universe size `P` in the denominator of both
`.mitopps_universe` (`P - 1`) and `.mitopps_query` (`P`). Adding two pathways
changes every human score; the saved 142-row universe cannot be extended.

**The control came first and it passed.** Rebuilding the human **142** universe
from the linear matrices reproduces the saved `$mitopps_universe` **bit-equal,
`max |delta| = 0`, in both cohorts** - including TCGA's, which was originally
built in the validation repo. That is what licenses reading anything from the
144 run.

What the two extra pathways cost the primary's percentile:

| cohort | median at 142 | median at 144 | shift |
|---|---|---|---|
| TCGA | 46.13 | 45.49 | -0.64 |
| SCAN-B | 50.35 | 50.35 | -0.005 |

Negligible, as expected for a 2-in-142 change to the universe.

**The human split is an ANALOGUE, not a translation.** It is human-native, from
`data/genesets_celldeath_human/cell_death_genes_consolidated.csv`, whose
`effect` column is a first-class human call. **No ortholog function is called;
the tripwire is clean.** But it is not the mouse's set:

| | mouse | human |
|---|---|---|
| PRO | 25 | **25** |
| ANTI | **9** | **6** |
| unclassified or unannotated | 0 | 4 |
| coverage of `Apoptosis` | 34 of 34 | 31 of 35 |

PRO agrees exactly. **ANTI does not, and the human ANTI set is 6 genes, just
above the `min_genes = 3` floor.** So the readout is **descriptive only and no
reading rule is pre-specified for it**, as fixed in the plan.

---

## 3. G2 FAILED - the size gate, and the remedy

This is the finding that changes how everything below it reads.

| panel | `rho_panel` | 95% CI | within-sample mean |
|---|---|---|---|
| TCGA, intersected | **0.536** | 0.400 to 0.651 | 0.071 |
| SCAN-B, intersected | **0.571** | 0.444 to 0.679 | 0.071 |

`rho_panel` is `rho(rowMeans(percentile), log(size))` across the 142 pathways -
"across the whole intersected panel", which is what the rule asked for.
`rho_within_sample_mean` is the mean of the same correlation computed inside
each sample; it is reported beside it and is **not** the gate statistic. The two
answer different questions and differ by a factor of eight, which is why the
gate names which one it is.

**`OXPHOS subunits` carries 88 genes as scored** - 89 in the catalogue, 88
present in the TCGA matrix - **against a panel median of 11. It is the 94th
percentile of the panel by size.** A large set's mitoPPS rank is partly a
statement about its size, so the primary's position against the panel is
substantially structural.

**The remedy was pre-specified in the script header, not chosen after seeing
the number.** Two size-matched bands around 88 genes:

| band | n | gene range |
|---|---|---|
| f = 1.5 | 12 | 59 to 121 |
| f = 2.0 | **22** | 44 to 155 |

R1 is scored on the larger band, `size-matched (f=2.0, n=22)`.

### 3.1 A limitation of the remedy, recorded rather than argued away

MitoCarta is hierarchical, and size-matching concentrates the primary's own
family. Of the 142-pathway panel, **11 members (7.7%) are OXPHOS-family**. Of
the f=2.0 band, **5 of 22 (23%) are** - `Complex I`, `Complex IV`, `OXPHOS`,
`OXPHOS assembly factors` and `OXPHOS subunits` itself. The band also carries
`Protein import, sorting and homeostasis` alongside both of its children. So the
primary competes against its own parent, its own children and its own siblings,
at three times the density of the full panel.

**What this does and does not damage.** It means the band is not a set of
independent comparators and the absolute percentile on it should not be read as
"OXPHOS is at the 39th percentile of mitochondrial biology". It does **not**
differentially bias the mouse-human comparison: the band is defined once on
human gene-set sizes and applied to the same 22 named pathways in all six arms,
so whatever distortion it introduces is common to every arm and the *relative*
placement survives.

---

## 4. G3 - the mouse reference falls, but weakly

The gate was: if the OXPHOS percentile does not fall from 6W WT to 12W WT, the
analysis stops. **It falls, in every panel definition, and every interval covers
zero.**

| panel | pathway | 6W WT | 12W WT | diff | 95% CI |
|---|---|---|---|---|---|
| intersected (142) | `OXPHOS subunits` | 56.0 | 18.7 | **-37.3** | -78.2 to +3.9 |
| intersected (142) | `OXPHOS` | 53.9 | 26.1 | -27.8 | -69.0 to +5.6 |
| own full panel | `OXPHOS subunits` | 55.9 | 18.4 | -37.5 | -78.5 to +3.8 |
| apoptosis-split (144) | `OXPHOS subunits` | 55.9 | 18.4 | -37.5 | -78.5 to +3.8 |
| size-matched (f=1.5) | `OXPHOS subunits` | 58.3 | 8.3 | -50.0 | -91.7 to +8.3 |
| size-matched (f=2.0) | `OXPHOS subunits` | 54.5 | 6.8 | **-47.7** | -86.4 to +11.4 |

**The direction is consistent across all six definitions and the magnitude is
large. The precision is not there.** n = 6 per arm. The script's `READABLE` flag
evaluates `TRUE` because the fall is in the pre-specified direction and the band
is large enough to rank in, **but the two reference positions are not well
separated, and every sentence below inherits that.** The reading rules were not
relaxed after this was seen.

**This is the first time the mouse reprioritisation has been asked to survive
translation into a within-panel rank.** It was established on composites scored
against expression-matched nulls. Direction survives; precision does not.

---

## 5. R1 and R2 - the placements

**Scored on `size-matched (f=2.0, n=22)`, because G2 fired.** The intersected
panel is shown beside it as the pre-specified sensitivity, never as the primary.

| arm | n | size-matched f=2.0 | IQR | intersected (142) | IQR |
|---|---|---|---|---|---|
| mouse 6W WT | 6 | **54.5** | 13.6 to 88.6 | 56.0 | 21.7 to 88.7 |
| mouse 12W WT | 6 | **6.8** | 2.3 to 21.6 | 18.7 | 11.8 to 22.9 |
| mouse 6W Myc+ | 6 | 86.4 | 63.6 to 92.0 | 81.7 | 71.3 to 84.7 |
| mouse 12W Myc+ | 6 | 31.8 | 14.8 to 69.3 | 36.3 | 21.8 to 57.6 |
| **TCGA** | 1,095 | **38.6** | 11.4 to 75.0 | 46.1 | 23.6 to 73.6 |
| **SCAN-B** | 3,207 | **47.7** | 15.9 to 75.0 | 50.4 | 28.5 to 71.5 |

- **R1, TCGA: BETWEEN the two mouse positions, and near neither.** Reported as
  such and **not assigned to the nearer end**, which the rule forbade in
  advance.
- **R1, SCAN-B: near the 6W position.**
- **R2: THE COHORTS DISAGREE - and R2 named disagreement as the result.** They
  are never averaged. On the intersected sensitivity the gap narrows (46.1
  against 50.4) and both drift toward the middle, but neither cohort approaches
  the 12W position on any panel.

**Neither human cohort sits at the lowered 12W position.** That is the whole of
the positive content, and it is one descriptive observation.

Two things visible in the table that were not part of any gate, recorded because
they are coherent and free:

1. **The `Myc+` arms sit above their WT counterparts at both ages** - 86.4
   against 54.5 at 6W, 31.8 against 6.8 at 12W. The transgene raises the
   respiratory chain's rank at both timepoints.
2. **The developmental fall happens in both genotypes**, 6W to 12W, WT and
   `Myc+` alike. G3 was WT-only by design; the `Myc+` arms were not used to
   define the reference and are not being promoted to evidence here.

The `OXPHOS` umbrella tracks the primary throughout as the pre-specified
secondary (TCGA 47.5, SCAN-B 51.1, mouse 6W WT 53.9, 12W WT 26.1 on the
intersected panel).

---

## 6. The apoptosis readout - reported, not scored

**No reading rule was pre-specified and none is applied.** From the
apoptosis-split 144 panel:

| pathway | 6W WT | 12W WT | 6W Myc+ | 12W Myc+ | TCGA | SCAN-B |
|---|---|---|---|---|---|---|
| `Apoptosis-ANTI` | 78.8 | 90.6 | 5.2 | 9.0 | 48.3 | 46.9 |
| `Apoptosis-PRO` | 74.3 | 84.4 | 10.4 | 12.2 | 55.9 | 57.3 |
| PRO minus ANTI | -4.5 | -6.2 | +5.2 | +3.1 | +7.6 | +10.4 |

**The split does not separate.** In the mouse the two halves move together and
massively: both WT arms sit high on both halves, both `Myc+` arms collapse to
single digits on both. What the split reports is a genotype effect on the rank
of the whole apoptosis pathway, not a shift in the balance between pro-death and
pro-survival members.

**The PRO-minus-ANTI gap is a within-species quantity only.** The mouse gap uses
a 9-gene ANTI set and the human gap a 6-gene one, so the last row must not be
read across the species boundary. Within each species the direction is
consistent (WT slightly ANTI-leaning, `Myc+` and both human cohorts slightly
PRO-leaning) and the magnitudes are small against n = 6 mouse arms.

---

## 7. What would falsify this

Written down now, before the next analysis, as this repo requires.

1. **The mouse fall.** A rank comparison with more animals, or a paired
   within-litter design over the same 24, that places the 12W WT median at or
   above the 6W WT median. Equally: showing the fall disappears when the
   comparator band is both size-matched **and** family-free, which would make it
   a statement about the OXPHOS family's internal ordering rather than about the
   respiratory chain's priority in the compartment.
2. **The human placement.** A human breast cohort scored on this panel whose
   median lands near the 12W position. Or - the more likely and more
   interesting failure - showing that the human percentile is carried by
   non-tumour content. **CLAUDE.md trap 2 applies directly**: in TCGA
   `rho(OXPHOS subunits, purity) = 0.214` and adipose is OXPHOS-high, so a
   tumour-cell-intrinsic percentile could sit well below the whole-tumour one.
   **TCGA has a purity estimate and SCAN-B does not, which is also the leading
   candidate explanation for R2's disagreement.**

**The named next test, and it is cheap:** re-read the TCGA percentile in
high-purity tumours only, against the same size-matched band. If TCGA moves
toward SCAN-B, the cohort disagreement is a purity artefact. If it moves away
from both mouse positions, the human placement was never about tumour cells.

---

## 8. Constraints observed

- **Rank only.** No mitoPPS value crosses the species boundary in the object,
  the table or the figure. Verified after the run: every numeric column in every
  saved data frame is a percentile, a rank correlation, a rank difference or a
  count.
- **No pooling.** Mouse and human were never scored together. TCGA and SCAN-B
  are separate everywhere and never averaged. Species is a cohort (trap 6).
- **The mouse artefact is read-only input.** md5 `b8a125af4bc0d5909ad03f6e126e2890`
  re-checked against `data/from_myc_mouse/README.md` at run time and still
  `CURRENT` against source after it. **Mouse mitoPPS was not recomputed here**
  (`$settings$mouse_artefact$recomputed = FALSE`) and the mouse repo was not
  read at run time.
- **No ortholog call anywhere.** Tripwire clean. The human apoptosis split is
  human-native (trap 11).
- **SCAN-B scored through `scanb_pheno.rds$symbol_map`** (trap 7).
- **`mt-` / `MT-` genes held separately in both species** and never pooled with
  nuclear-encoded subunits (trap 8, standing convention).
- **N3.** Nothing here is called primed. Every number is a transcript rank.
- **Nothing was written to `myc_human_validation` or `myc_mouse`.** The mouse
  repo's working tree and HEAD (`129b0c0`) are exactly as found.

## 9. Provenance and reproduction

- Script: `scripts/E24_mitopps_rank_comparison.R`, seed 1, 5,000 bootstrap
  draws.
- Mouse input: `data/from_myc_mouse/mitopps_scores.rds`, produced by
  `myc_mouse` `scripts/08_mitoPPS_analysis.R` at commit `0b9b28e`. See that
  directory's README for why the repo HEAD does not pin the file.
- Human inputs: the saved 142-pathway universes, TCGA's from
  `data/from_validation/tcga_brca_mito_scores.rds` and SCAN-B's from
  `results/scanb_scores.rds`; the linear DESeq2-normalised matrices from
  `data/from_validation/`, used to rebuild the 142 universe as the A2 control
  and then to build the 144 one. SCAN-B symbols are remapped through
  `scanb_pheno.rds$symbol_map` before scoring.
- Saved: `results/mitopps_rank_comparison.rds` (21 elements),
  `outputs/tables/E24_ranks_by_group.csv`,
  `outputs/mitopps_rank/E24_oxphos_percentile.pdf`.
