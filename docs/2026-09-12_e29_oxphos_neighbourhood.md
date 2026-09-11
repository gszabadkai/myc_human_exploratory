# E29 - the map: position matches are coincidental, and the panel is not organised around OXPHOS

**A MAP, NOT A READOUT.** No pairing produced here may be used as a comparator
without a separately declared validation across the full estimator ladder, both
cohorts and the mouse. That sentence is in the saved object, the verdict and
every figure caption.

**2026-09-12. EXPLORATORY, POST-HOC, DESCRIPTIVE.** Nothing pre-registered,
nothing a hypothesis test. R1 to R5 were fixed before any number so the outcome
could not be chosen afterwards, and none was a prediction.

Script: `scripts/E29_oxphos_neighbourhood.R`. **Sourced by the author
2026-09-12**; the author's run reproduces **all 27 saved elements identically**,
so every number here is read from the object. `E25`, `E27` and `E28` are not
modified and nothing is re-scored.

---

## 0. The answer in seven lines

| | |
|---|---|
| **Is mitoPPS a closed composition?** | **No, and this is the useful fact.** The per-sample sum ranges 135.8 to 184.6, and pairwise rho between the 142 scores sits on zero. The design's premise is tested away. |
| **Is anything constrained?** | **Yes, one layer up.** Percentiles sum to exactly 7100 per group, so the change vector sums to zero and the neutral point is -0.42 to -0.18, not zero. |
| **Do similar positions mean similar configurations?** | **No. R1 is negative.** **Not one of 28 pairs** exceeds its rho-matched null. Every overlap is what the whole-ordering correlation already predicts. |
| **What about the named pair?** | 12W wild-type and TCGA MYC-high share **7 of 20** neighbours inside a null interval of **[1, 7]**. It does not exceed prediction either. |
| **Are the off-diagonal pathways a coherent set?** | **No. Scattered.** 81 distinct pathways fill the extreme-residual slots across 12 pairs; only 8 recur in half or more. |
| **Is there an opposing programme?** | **A count, not a set, and the claim is declined.** 52 move with the reference in both species, 51 oppose in both, 38 are species-specific. |
| **Does the same thing track the reference in both species?** | **Among pathways that actually move, yes.** 84 of 101 keep their sign, against 19 of 40 among those that barely move. |

**The one-sentence reading: the panel is not organised around the respiratory
chain in any way this map can resolve, and two groups sitting at the same
position does not mean they sit in the same place.**

---

## 1. Both controls pass

1,136 positions reproduce `E25` and 88 reproduce `E27`, **max |delta| = 0** in
both.

---

## 2. Two premises about the score: one tested away, one found

**The design was written expecting mitoPPS to be a closed composition** -
ratio-based, centred at 1 per pathway, therefore inducing spurious negative
correlation between pathways. **That concern is reasonable a priori and it is
not borne out.**

| cohort | per-sample sum | pairwise rho, median | percent negative |
|---|---|---|---|
| TCGA | 142.0, range 135.8 to 184.6 | -0.005 | 51.0 |
| SCAN-B | 142.0, range 137.2 to 169.7 | -0.009 | 52.1 |
| mouse | 142.0, range 140.2 to 150.3 | +0.004 | 49.8 |

**The sum is not fixed, so the panel is not closed and the standard
compositional objection does not apply.** Across 10,011 pathway pairs the
correlation distribution sits on zero in every cohort. **No Phase 2 note records
that the sum is unconstrained, and that is the durable fact here** - more useful
than the correction that prompted it.

R5 says the diagnostic goes into `CLAUDE.md`'s trap list **if strongly
negative**. It is not, so it is recorded here and **`CLAUDE.md` is not edited.**

**What is constrained is the percentile transform.** The 142 percentiles sum to
exactly 7100 in every group, so the change vector sums to exactly zero and the
neutral point for aligned change is `-|delta reference| / 141`, running **-0.42
to -0.18** across the six contrasts. It is drawn on the figures rather than
assumed. **It is small here only because the forcing scales as 1/(P-1): on a
20-pathway panel the same reference change would force -2.7 percentile points
and would have to be handled rather than noted.**

**PART B's contrast-profile design is kept, for a different reason than it was
written for:** it is the estimand every other Phase 2 contrast uses, which keeps
`E29` comparable with `E25`, `E27` and `E28`.

---

## 3. Where the reference sits, and what its neighbourhood is made of

| group | percentile | rank of 142 | neighbours retained under resampling |
|---|---|---|---|
| mouse 6W WT | 60.9 | 87 | **5.5** of 20 |
| mouse 12W WT | 7.4 | 11 | 12.9 |
| mouse 6W Myc+ | 89.1 | 127 | 12.2 |
| mouse 12W Myc+ | 32.7 | 47 | **7.6** |
| TCGA MYC-low | 16.5 | 24 | 12.5 |
| TCGA MYC-high | 68.0 | 97 | 11.7 |
| SCAN-B MYC-low | 15.1 | 22 | 15.7 |
| SCAN-B MYC-high | 74.3 | 106 | 14.7 |

**The n = 6 caveat is a number rather than a sentence.** Mouse groups retain 5
to 13 of 20 neighbours under resampling against 12 to 16 in the human, and mouse
6W wild-type at 5.5 is the least stable neighbourhood in the study. Weight the
mouse panels accordingly.

**The neighbourhood is largely the reference's own MitoCarta family.** The
OXPHOS family is 20 of 142 panel members, 14%, but fills **10 to 40%** of the
neighbourhood depending on the group. **Two groups therefore share neighbours
partly because the OXPHOS branch sits near OXPHOS in both, which is a property
of the ontology and not of biology.** That is measured and reported before any
overlap is read.

---

## 4. R1 and R2 - the null is what changed the answer

**Not one of the 28 pairs exceeds its rho-matched, position-matched null.**

The plain random null would have told a different story, and the difference is
the whole point of building the better one:

| pair | shared | plain null | rho-matched null |
|---|---|---|---|
| mouse 12W WT and 6W Myc+ | 13 | 2.84 | **13.6 [11, 16]** |
| TCGA MYC-low and SCAN-B MYC-low | 10 | 2.84 | **9.6 [6, 13]** |
| mouse 12W WT and SCAN-B MYC-low | 8 | 2.84 | **6.6 [3, 10]** |
| mouse 12W WT and TCGA MYC-high | 7 | 2.84 | **3.6 [1, 7]** |

Against the plain null every one of these looks like strong conservation.
Against a null that knows the two orderings' correlation and the reference's
position in each, **all four are exactly what is predicted.**

**R1 is therefore negative: position matches are coincidental and must not be
read as biological similarity.**

**R2's named pair does not escape it.** Mouse 12W wild-type and TCGA MYC-high
share 7 of 20 neighbours, at the top edge of a [1, 7] interval but not beyond
it. Their whole orderings sit at rho -0.51 and their residuals run -96 to +98
percentile points. **The identical central-dogma gap that `E28` found for these
two groups came from very different positions, and this is what sits behind it:
two orderings that disagree, sharing no more than their disagreement predicts.**

---

## 5. R2b - the off-diagonal pathways are scattered

Across the twelve ordering scatters, **81 distinct pathways** fill the
extreme-residual slots and **only 8 recur in half or more**. Those eight are
`Biotin utilizing proteins`, `Branched-chain amino acid dehydrogenase complex`,
`Calcium uniporter`, `Chaperones`, `Cytochromes`, `Molybdenum cofactor synthesis
and proteins`, `Transcription` and `mtDNA replication`.

**The orderings differ diffusely and no subset carries the difference.** Eight
recurring names out of 81 is not a programme, and they are not named as one.

---

## 6. R3 and R4 - co-movement

**R3, read three ways as the rule requires:** 52 pathways move with the
reference in both species, 51 oppose it in both, 38 are species- or
state-specific.

**The coherence claim is declined, and that is a departure from what the rule
would have allowed.** The rule offered "a coherent opposing set is a map entry
worth naming" on a count threshold. A count of 51 is not a coherent set:
coherence would require the opposing pathways to be related to each other, which
this script does not test, and **R2b had just found the off-diagonal structure
scattered.** So the sign distribution is reported and no opposing programme is
named.

**R4 is stronger than its headline, and the headline is where it is weak.** 103
of 141 pathways keep their sign between mouse and human, 73% against 70
expected under independence. But that figure mixes two populations:

| pathways | sign agreement |
|---|---|
| **101 with mean magnitude 10 or more in both species** | **84 of 101, 83%** |
| 40 with mean magnitude under 10 in either | 19 of 40, **48%** |

**The near-zero pathways agree at exactly chance, as noise should.** The
conservation is carried entirely by the pathways that actually move. That split
is descriptive, computed from the saved object, and is not part of R4's rule.

---

## 7. Two consistency checks, one of which is half tautology

**`E27`'s subunits-against-assembly-factors split reappears at finer grain, and
the gene overlap decides which half of it is informative.**

| OXPHOS-branch member | genes shared with the reference | behaviour |
|---|---|---|
| `CI`, `CII`, `CIII`, `CIV`, `CV subunits` | **100% of each** | all move **with** the reference |
| `Complex I` to `Complex V` | 38 to 79% | all move **with** it |
| `CI`, `CII`, `CIII`, `CIV`, `CV assembly factors` | **0%** | 2 oppose in both, 2 species-specific, 1 moves with |
| `OXPHOS assembly factors` | 1 of 68 | **opposes in both** |
| `Respirasome assembly` | 1 of 4 | **opposes in both** |
| `mtDNA-encoded OXPHOS subunits` | **0 of 13** | **opposes in both, most strongly** |

**The first two rows are largely tautological** - the subunit sets are wholly
contained in the reference, so "moving with it" is arithmetic. **The informative
half is the gene-disjoint sets**, and of the seven of them four oppose the
reference in both species, two are species-specific and one does not. **The
respiratory chain and the machinery that assembles it move apart, and the
mtDNA-encoded transcripts move furthest apart of all**, which is `E27`'s R4b and
trap 8 arriving independently.

**`Protein import, sorting and homeostasis` behaves as `E28` said it would.** It
ranks 106 to 121 of 141 on aligned change across all six contrasts, that is
firmly among the pathways moving with the reference. `E28` found it the closest
comparator to OXPHOS in the mouse wild-type window, and closeness and
co-movement are the same statement. **This is a consistency check on a known
observation, not a new finding.**

---

## 8. What this is for, and what it is not

**It is a map.** It says where the respiratory chain sits in each group's
ordering, what is around it, how stable that is, and which pathways move with or
against it. Every one of those is descriptive.

**It is not a source of comparators.** Nothing here may be promoted to a readout
without a separately declared validation across the full estimator ladder, both
cohorts and the mouse. `E28` is why: two pre-declared comparators failed, and
scanning 141 candidates and keeping the ones that behave would be selecting on
the outcome.

**Three specific things it must not be used to say.**

1. **That any two groups share a mitochondrial configuration.** R1 is negative
   across all 28 pairs.
2. **That a set of pathways opposes the respiratory chain.** The counts are
   reported; the set is not established.
3. **Anything about redox.** No redox annotation was applied and none was
   pre-declared. Names were read after the ranking was produced. Any redox
   reading is a separate, separately declared analysis.

**And one thing about the figures.** A point on a scatter's diagonal shares a
**rank** in both orderings. It says nothing about its mitoPPS value, which does
not cross the species boundary. Percentiles are uniform within a group by
construction, so C1's pale band is a reference range and not a density. The two
panels whose whole-ordering rho is near zero are labelled in the figure as
calibration panels, so an absence of structure reads as what no relationship
looks like rather than as a failed panel.

---

## 9. Constraints observed

- **A map, not a readout**, in the object, the verdict and every caption.
- **Nothing re-scored.** Re-reading and re-ranking existing objects.
- **No mitoPPS value crosses the species boundary.** Percentiles and percentile
  changes only; cross-species comparison is by pathway name.
- **TCGA and SCAN-B separate**, never averaged.
- **No ortholog call.** Tripwire clean.
- **No FDR**: a ranking exercise, not a screen.
- The mouse timeline is **batch-confounded**; reported in full, nothing
  withdrawn.
- **N3.** Transcript associations. "Primed" appears nowhere.

## 10. Provenance

- Script: `scripts/E29_oxphos_neighbourhood.R`, seed 1, 5,000 bootstrap draws
  and 5,000 null simulations.
- Reads `results/mitopps_group_rankings.rds`, `results/category_readouts.rds`,
  `results/gap_readout.rds`, `data/from_myc_mouse/mitopps_scores.rds` (md5
  `b8a125af4bc0d5909ad03f6e126e2890`), the validation snapshots and
  `data/mitocarta_human/Human.MitoCarta3.0.xls`.
- Saved: `results/oxphos_neighbourhood.rds` (27 elements),
  `outputs/tables/E29_aligned_change.csv`,
  `outputs/tables/E29_diagonal_residuals.csv`, and seven figures under
  `outputs/mitopps_rank/`.
