# E16 - G1: is `ox_rel` a different ruler, and does the configuration survive it?

**2026-09-04. EXPLORATORY. Nothing here is pre-registered.**

Script: `scripts/E16_respiratory_rulers.R`. Written this session, **not yet run
by the author** - the numbers below are from a dry run with the three output
paths redirected to the scratchpad, and the repo was left untouched (`git status`
clean apart from the new script). Sourcing it in Positron takes a few seconds:
it reads six saved objects and builds two composites. That source is what puts
`results/respiratory_rulers.rds` and the two CSVs on disk.

**This is a MEASUREMENT CHECK, not an analysis.** No interaction was fitted, no
MYC stratification made, no subtype split taken, no new endpoint or ratio
introduced. mitoPPS was not recomputed - the scores exist for both cohorts and
were read.

**N3 throughout.** Every number here is a transcript correlation. The word
"primed" appears nowhere as a description of a tumour.

---

## 0. The answer in five lines

| | |
|---|---|
| **R1** are `ox_rel` and `ox_ppd` interchangeable? | **NO.** They correlate at 0.913 / 0.885 - above the 0.85 bar - but their MYC entanglements differ by up to 0.089, above the 0.05 bar. **Distinct instruments.** |
| **R2** does human reproduce the mouse ordering `ox_rel < ox_ppd < ox_lvl`? | **NO. 2 of 6 cells.** It holds on the MSigDB signature and fails on MYC mRNA and on the regulon estimator, in both cohorts. **The mouse's case for `ox_rel` does not transfer.** |
| **R3** does the three-gene reversal hold on `ox_rel`? | **YES, both cohorts.** `BBC3` +0.312 / +0.145, `BCL2L1` +0.360 / +0.306, `MCL1` -0.159 / -0.221. |
| **R4** does `BCL2L1` minus `BBC3` stay positive? | **YES, and on every ruler tried.** +0.049 / +0.161 on `ox_rel`, and positive in 8 of 8 ruler-by-cohort cells. **The load-bearing quantity does not depend on the ruler.** |
| **R5** is the configuration ruler-stable? | **Mostly, and `ox_rel` is the least stable of the three.** Rank agreement with the incumbent 0.867 / 0.783, against 0.930 / 0.902 for mitoPPS. |

**Primary ruler for the human arm going forward: the incumbent
`OXPHOS subunits` GSVA stays primary, and `ox_rel` is added beside it as a named
third instrument, never substituted for it.** Section 8 says why, and says
plainly which part of that is a rule and which part is a judgement.

---

## 1. What `ox_rel` is here, and how the recipe was matched

Read from `git -C /Users/gs/G/data/MK_myc_2022/myc_mouse show
HEAD:scripts/48_gate_model_mouse_verification.R`, HEAD `e348dd8`, branch
`paper-final`, PART A. The mouse repo was read only; nothing was written to it,
and **no ortholog function was called in either direction** - the tripwire is
clean.

```
mito_all  = union of every MITOCARTA_* set in the mouse GMT
ox_sub    = MITOCARTA_OXPHOS_SUBUNITS
rest_mito = setdiff(mito_all, ox_sub)
comp_e(e) = colMeans(t(scale(t(L[e, ]))))      L = log2(DESeq2-normalised + 1)
ox_rel    = comp_e(ox_sub) - comp_e(rest_mito)
```

Three questions the audit had to settle before a line was written, each answered
against the mouse GMT rather than assumed:

**(a) Where do the 13 mtDNA-encoded genes go?** *Excluded from the numerator
only.* The mouse `MITOCARTA_OXPHOS_SUBUNITS` set carries 89 symbols and not one
`mt-` gene. But `mito_all` is the union of every `MITOCARTA_*` set, and that
union includes `MITOCARTA_ALL`, `MITOCARTA_MTDNA_ENCODED` and
`MITOCARTA_OXPHOS_MT` - so all 13 sit in `rest_mito`, the **denominator**. That
is consistent with this repo's standing convention, which forbids *pooling* the
13 with the nuclear subunits and says nothing about a compartment-wide
background. The script carries `ox_rel_nomt` as a sensitivity rather than
arguing the point: **Spearman 0.9997 / 0.9995** against `ox_rel`. The decision
does not matter, and now that is measured.

**(b) Is `mito_all` the inventory or the union of the pathways?** *The
inventory.* In the mouse the two coincide - the GMT's `MITOCARTA_ALL` is 1,140
and so is the union of every `MITOCARTA_*` set, because the GMT also carries
`MITOCARTA_UNASSIGNED`. In human they do **not**: the sheet-2 inventory is 1,136
and the union of the 149 MitoPathways is 1,035, leaving 101 genes in no pathway.
The inventory is the matching object. It is also the same pinned set the
`__MITOSTRIP` estimators were stripped with, which is why `ox_rel` shares
**exactly 0 genes** with `FELSHER__MITOSTRIP` - on both halves, by construction.

**(c) What scale is the composite on?** `log2(linear DESeq2-normalised + 1)`,
which is the mouse's `L` and this repo's gene-level scale. It is **not**
`zmean_arms`, which E02 builds on the VST. Those two are not the same object:
Spearman between them is 0.9998 / 0.9995, high but not 1, and `ox_lvl` is built
here rather than read so the recipe matches.

### What is actually being measured

| | nominal | TCGA used | SCAN-B used | coverage |
|---|---|---|---|---|
| numerator `OXPHOS subunits` | 89, **0 mtDNA-encoded** | 88 | 88 | **0.989 / 0.989** |
| denominator, MitoCarta minus those | 1,047, **13 mtDNA-encoded** | 1,034 | 1,033 | 0.988 / 0.987 |

`COX8C` is the single missing gene, and it is missing from **both** matrices.
Unharmonised, SCAN-B's numerator would cover **0.775** instead of 0.989 and
Complex V would lose its F1 head and c-ring - `scanb_scores.rds$symbol_map` is
used throughout and the coverage is asserted.

For the shape of it and **never as a value**: the mouse built `ox_rel` from 87
subunits against 999 background genes.

### What overlaps what

| set | n | with the numerator | with the denominator |
|---|---|---|---|
| `FELSHER__MITOSTRIP` | 61 | **0** | **0** |
| `MYC_UP.V1_UP__MITOSTRIP` | 178 | **0** | **0** |
| `HALLMARK_MYC_TARGETS_V1__MITOSTRIP` | 177 | **0** | **0** |
| `M_b__MITOSTRIP` (CollecTRI) | 811 | **0** | **0** |
| `PROLIF_DISJOINT` | 318 | 0 | **10** |

The ten are `AK2 DUT MTHFD2 PAICS PRDX4 TBRG4 UNG DTYMK LIG3 POLQ`. Dropping
them from the denominator moves `ox_rel` by nothing: Spearman **0.9999 /
0.9999**.

---

## 2. Check 1 - do the rulers agree?

Spearman, **within each cohort and never across** - mitoPPS's baseline is
composition-dependent and GSVA is cohort-relative, so a value comparison between
cohorts is meaningless here and a value comparison between species is forbidden.

| pair | TCGA | SCAN-B |
|---|---|---|
| `ox_rel` vs `ox_ppd` | 0.913 | 0.885 |
| `ox_rel` vs `ox_lvl` | 0.941 | 0.924 |
| `ox_rel` vs `ox_gsva` | 0.922 | 0.902 |
| `ox_ppd` vs `ox_lvl` | 0.889 | 0.838 |
| `ox_ppd` vs `ox_gsva` | 0.880 | 0.825 |
| **`ox_lvl` vs `ox_gsva`** | **0.986** | **0.975** |

Two readings. **The four rulers are one family**, none below 0.825. And **the
absolute level and the GSVA score are all but the same ruler** (0.986 / 0.975),
which retrospectively justifies the G1 prompt's own framing: the exploratory
arm's GSVA and the validation arm's absolute level really were both `ox_lvl`
analogues, and now that is a number rather than an assertion.

---

## 3. Check 2 - separability. The check that decides G3.

**RAW marginal correlations, not proliferation-adjusted** - entanglement is the
quantity being measured, so partialling it out would answer a different
question. `prolif` is a column of the table instead. This is the only section of
the script where the standing adjustment does not apply.

The three set-based estimators are `__MITOSTRIP`, matching what the mouse did to
its two signature scores. **N7**: MYC mRNA is not MYC activity, both are
reported, neither collapses into the other.

Spearman (the standing measure here):

| cohort | ruler | myc_mRNA | myc_msigdb | myc_regulon | myc_felsher | myc_lowent | prolif |
|---|---|---|---|---|---|---|---|
| TCGA | `ox_lvl` | -0.016 | 0.525 | 0.279 | 0.367 | 0.471 | 0.208 |
| TCGA | `ox_ppd` | -0.014 | 0.497 | 0.292 | 0.306 | 0.396 | 0.240 |
| TCGA | **`ox_rel`** | -0.103 | **0.442** | 0.333 | **0.259** | **0.329** | **0.159** |
| TCGA | `ox_gsva` | -0.032 | 0.558 | 0.240 | 0.388 | 0.479 | 0.237 |
| SCAN-B | `ox_lvl` | -0.003 | 0.552 | 0.190 | 0.422 | 0.452 | 0.383 |
| SCAN-B | `ox_ppd` | -0.025 | 0.484 | 0.175 | 0.305 | 0.338 | 0.416 |
| SCAN-B | **`ox_rel`** | -0.048 | **0.470** | 0.254 | **0.328** | 0.338 | **0.334** |
| SCAN-B | `ox_gsva` | -0.011 | 0.577 | 0.117 | 0.427 | 0.430 | 0.394 |

The same table in **Pearson**, which is what the mouse table is, is in the saved
object and in `outputs/tables/E16_ruler_separability.csv`. It moves nothing: the
orderings are identical on every row that matters. The two are reported side by
side as **orderings only** - the mouse values come from 24 animals where
everything correlates with everything, these from 1,095 and 3,207 tumours, and a
species is a cohort.

### 3.1 The result, and it is not the one the mouse predicts

**`ox_rel` does not order the way the mouse says it does.** The mouse has
`ox_rel < ox_ppd < ox_lvl` on all three of its estimators. In human that holds
on **2 of 6** ruler-cohort cells - the MSigDB signature, in both cohorts - and
fails on the other four:

| cohort | against | `ox_rel` | `ox_ppd` | `ox_lvl` | as mouse? |
|---|---|---|---|---|---|
| TCGA | myc_mRNA | 0.103 | 0.014 | 0.016 | no |
| TCGA | **myc_msigdb** | **0.442** | 0.497 | 0.525 | **yes** |
| TCGA | myc_regulon | 0.333 | 0.292 | 0.279 | no |
| SCAN-B | myc_mRNA | 0.048 | 0.025 | 0.003 | no |
| SCAN-B | **myc_msigdb** | **0.470** | 0.484 | 0.552 | **yes** |
| SCAN-B | myc_regulon | 0.254 | 0.175 | 0.190 | no |

Values are `|r|`, because entanglement is a distance from zero and a human ruler
may lean the other way.

**The pattern is coherent and it replicates: `ox_rel` buys separability from
signature-based MYC estimators and loses it against the regulon estimator.**
Against the incumbent ruler, gene set by gene set:

| against | `ox_rel` minus `ox_gsva`, TCGA | SCAN-B |
|---|---|---|
| myc_msigdb | **-0.116** | **-0.107** |
| myc_felsher | **-0.129** | **-0.099** |
| myc_lowent | **-0.150** | **-0.092** |
| prolif | **-0.078** | **-0.060** |
| myc_regulon | **+0.093** | **+0.137** |
| myc_mRNA | +0.071 | +0.037 |

All six signature cells go down, and so do both proliferation cells; both
regulon cells go up.
`myc_mRNA` moves up too, but both values there are near zero (0.103 against
0.032; 0.048 against 0.011) and the direction of a near-zero pair is not worth
leaning on.

**So "separability from MYC" is not a property of a ruler in human. It is a
property of the ruler-estimator pair.** That is CLAUDE.md trap 3 restated on a
new axis: the MYC estimators disagree about what MYC is by a factor of thirty in
proliferation entanglement, and they disagree about which respiratory ruler is
furthest from MYC as well. In the mouse the ordering was consistent across all
three estimators; in human it is not, and any G3 that names a ruler must name an
estimator in the same breath.

One property of `ox_rel` is unambiguous and worth keeping: **it is the least
proliferation-entangled of the four rulers, in both cohorts** (0.159 / 0.334
against 0.208-0.240 and 0.383-0.416).

---

## 4. Check 3 - the twelve on every ruler

Partial Spearman on `PROLIF_DISJOINT`, Fisher-z intervals with the
Bonett-Wright variance - E10's `.cor_block`, copied verbatim. **All twelve are
reported.** Selecting the three of interest after seeing them is the
grid-of-cells trap.

**Reproduction check, and it is exact.** `ox_gsva` and `ox_ppd` are the two axes
E10 already carries as `OXPHOS` and `OXPHOS_mitopps`. All **48 cells reproduce
E10 with max |delta| = 0**. That is the load block's positive control: a failure
there would be a broken input, not a new result.

**TCGA**, ordered by `ox_rel`:

| gene | side | `ox_rel` | `ox_ppd` | `ox_lvl` | `ox_gsva` |
|---|---|---|---|---|---|
| `BAD` | pro | +0.479 | +0.451 | +0.531 | +0.525 |
| `BID` | pro | +0.372 | +0.301 | +0.359 | +0.329 |
| **`BCL2L1`** | **anti** | **+0.360** | +0.327 | +0.404 | +0.426 |
| **`BBC3`** | **pro** | **+0.312** | +0.293 | +0.366 | +0.365 |
| `BIK` | pro | +0.289 | +0.255 | +0.327 | +0.338 |
| `BCL2A1` | anti | +0.143 | +0.081 | +0.041 | -0.004 |
| `BMF` | pro | -0.113 | -0.089 | -0.177 | -0.200 |
| `PMAIP1` | pro | -0.150 | -0.074 | -0.161 | -0.147 |
| **`MCL1`** | **anti** | **-0.159** | -0.216 | -0.227 | -0.241 |
| `BCL2L2` | anti | -0.230 | -0.149 | -0.164 | -0.145 |
| `BCL2` | anti | -0.269 | -0.100 | -0.203 | -0.160 |
| `BCL2L11` | pro | -0.324 | -0.299 | -0.331 | -0.314 |

**SCAN-B**, ordered by `ox_rel`:

| gene | side | `ox_rel` | `ox_ppd` | `ox_lvl` | `ox_gsva` |
|---|---|---|---|---|---|
| `BAD` | pro | +0.437 | +0.423 | +0.478 | +0.482 |
| **`BCL2L1`** | **anti** | **+0.306** | +0.305 | +0.347 | +0.350 |
| `BID` | pro | +0.200 | +0.067 | +0.137 | +0.091 |
| **`BBC3`** | **pro** | **+0.145** | +0.158 | +0.173 | +0.170 |
| `BIK` | pro | +0.120 | +0.066 | +0.210 | +0.213 |
| `BCL2A1` | anti | +0.070 | +0.004 | -0.062 | -0.107 |
| `BMF` | pro | -0.047 | -0.060 | -0.183 | -0.223 |
| `BCL2L2` | anti | -0.048 | -0.007 | -0.013 | -0.008 |
| `BCL2L11` | pro | -0.048 | -0.026 | -0.053 | -0.025 |
| `PMAIP1` | pro | -0.092 | +0.016 | -0.101 | -0.064 |
| **`MCL1`** | **anti** | **-0.221** | -0.257 | -0.299 | -0.292 |
| `BCL2` | anti | -0.232 | -0.030 | -0.130 | -0.074 |

The 3.1 reading of the synthesis document survives the ruler change intact:
**the family does not move as a block, functional class does not predict
position, and pro- and anti-apoptotic members interleave across the whole
range.** `BAD` is the top gene on every ruler in both cohorts; the bottom is
`BCL2L11` on all four rulers in TCGA, and `MCL1` on three of four in SCAN-B.

### 4.1 Where `ox_rel` disagrees with the incumbent, gene by gene

**11 of 12 genes keep their sign in both cohorts.** The single exception is the
same gene in both: **`BCL2A1`**, which is +0.143 / +0.070 on `ox_rel` and
-0.004 / -0.107 on `ox_gsva`. It is also the only cell in the whole grid that
changes sign under the infiltrate covariate (section 6). `BCL2A1` / A1 is
highly expressed in leukocytes, so a gene that moves with the composition of the
non-tumour compartment behaving differently under a compartment-relative ruler
is unsurprising - but that is an observation, not a claim, and nothing in this
note rests on it.

The largest per-gene shifts, `ox_rel` minus `ox_gsva`:

| TCGA | shift | SCAN-B | shift |
|---|---|---|---|
| `BCL2A1` | +0.148 | `BCL2A1` | +0.177 |
| `BCL2` | -0.109 | `BMF` | +0.176 |
| `BMF` | +0.088 | `BCL2` | -0.158 |
| `BCL2L2` | -0.085 | `BID` | +0.109 |
| `MCL1` | +0.082 | `BIK` | -0.093 |

**The three genes the mouse model names are among the more stable.** `BBC3`,
`BCL2L1` and `MCL1` shift by at most **0.082** in TCGA and **0.070** in SCAN-B,
against 0.148 and 0.177 for the largest shift in each cohort. `MCL1` is fifth on
the TCGA list above; the other two are outside both top-fives.

---

## 5. Check 3, continued - R5 and the load-bearing quantity

**Cross-ruler rank agreement over the twelve**, Spearman of the 12 rho values:

| pair | TCGA | SCAN-B |
|---|---|---|
| `ox_rel` vs `ox_ppd` | 0.958 | 0.860 |
| `ox_rel` vs `ox_lvl` | 0.930 | 0.853 |
| **`ox_rel` vs `ox_gsva`** | **0.867** | **0.783** |
| `ox_ppd` vs `ox_gsva` (reference) | 0.930 | 0.902 |

**No numeric threshold was pre-specified for R5 and none is invented here.** The
verdict is stated as a comparison against a reference this repo already carries:
GSVA-vs-mitoPPS agreement over these same twelve genes, 0.930 / 0.902, and
CLAUDE.md trap 5, which puts cross-instrument agreement over the 18
mitochondrial arms between 0.24 and 0.94.

Read that way: **the configuration is ruler-stable enough that a three-gene
claim is safe, and `ox_rel` is the least stable of the three against the
incumbent.** 0.783 in SCAN-B is the lowest number in the table and is the one to
quote if anyone claims the rulers are interchangeable in practice.

### `BCL2L1` minus `BBC3` - the quantity that carries the whole argument

| ruler | TCGA | SCAN-B |
|---|---|---|
| `ox_rel` | **+0.049** | **+0.161** |
| `ox_ppd` | +0.034 | +0.147 |
| `ox_lvl` | +0.038 | +0.175 |
| `ox_gsva` (incumbent, the published pair) | +0.061 | +0.180 |

**Positive in 8 of 8 ruler-by-cohort cells.** This is the strongest single
statement in the check. The gap is why `PRIME = BBC3 - BCL2L1` comes out flat on
OXPHOS, which is what explains the validation arm's pre-registered null - and
that explanation does **not** depend on which respiratory ruler is used. Swapping
in the mouse's bridge ruler does not rescue `PRIME`; it makes the gap slightly
smaller than the incumbent in both cohorts and leaves it positive throughout.

**N2 is untouched.** The sibling study is pre-registered, closed and frozen at
`d3ac60e`. Nothing here reopens it. This work *explains* its null and is written
that way.

---

## 6. Check 4 - infiltrate control, TCGA only

n = **1,007** with both purity and leukocyte fraction. **SCAN-B has no purity
estimate** (trap 2) and it is never imputed, so nothing here replicates and
nothing is written as if it did.

**Both rows are computed on the same 1,007 samples.** Comparing an adjusted
value at n = 1,007 against the section-4 value at n = 1,095 would confound the
covariate with the sample set, and the difference wanted here is the covariate
alone. Every ruler is carried, not just `ox_rel`: a shift all four share is a
property of the samples, not of the ruler.

| ruler | mean \|shift\| | max \|shift\| | sign changes |
|---|---|---|---|
| **`ox_rel`** | **0.010** | **0.045** | 0 |
| `ox_ppd` | 0.022 | 0.097 | 0 |
| `ox_lvl` | 0.027 | 0.103 | 0 |
| `ox_gsva` | 0.038 | 0.136 | 1 (`BCL2A1`, -0.009 to +0.128) |

**`ox_rel` is the most infiltrate-robust of the four**, by a factor of two to
four on the mean and two to three on the worst cell. The load-bearing
gap barely moves: `BCL2L1` minus `BBC3` goes +0.053 to +0.057 on `ox_rel`, and
+0.057 to +0.058 on the incumbent.

That robustness is the clearest practical argument `ox_rel` has in human, and it
is worth stating precisely: it is a **within-compartment share**, so anything
that scales the whole mitochondrial compartment - including a shift in how much
of the sample is tumour - largely cancels. It is one cohort and it cannot
replicate.

---

## 7. The verdicts, on rules fixed before any number was looked at

| rule | pass | evidence |
|---|---|---|
| **R1** are the rulers interchangeable? | **NO** | rho 0.913 / 0.885 (bar 0.85, passes); largest MYC entanglement gap **+0.089** on `myc_mRNA` in TCGA (bar 0.05, fails). Excluding the near-zero mRNA row the largest is still +0.079 (`myc_regulon`, SCAN-B), so the verdict does not hang on that row |
| **R2** does human reproduce the mouse ordering? | **NO** | full mouse ordering on **2 of 6** MYC cells; `ox_rel` is the least entangled in 2 of 6 |
| **R3** does the reversal hold on `ox_rel`? | **YES, both cohorts** | `BBC3` +0.312 [+0.256, +0.365] / +0.145 [+0.111, +0.179]; `BCL2L1` +0.360 [+0.306, +0.412] / +0.306 [+0.274, +0.338]; `MCL1` -0.159 [-0.216, -0.100] / -0.221 [-0.255, -0.188] |
| **R4** does the gap hold? | **YES** | +0.049 / +0.161 on `ox_rel`; positive on all four rulers in both cohorts |
| **R5** is the configuration ruler-stable? | **no pass/fail** | `ox_rel` vs `ox_gsva` rank rho 0.867 / 0.783 against the in-repo reference 0.930 / 0.902. No threshold was pre-specified, so this is reported as a comparison, not as a test |

**R1 deserves one sentence of care.** It fails on its *entanglement* clause, not
its correlation clause: the two rulers correlate above the bar and are still not
interchangeable, because what they are entangled with differs. That is the
useful shape of the result - "highly correlated" and "the same instrument" are
different claims, and the rule was written to separate them.

---

## 8. Which ruler is primary for the human arm

**The decision rules specify one branch and not the other, and that is said here
rather than papered over.** R2 reads: *"If yes, `ox_rel` is primary for G3
regardless of what R1 returned."* It fixes what happens when the human
reproduces the mouse ordering. It fixes nothing for the case where it does not.
Inventing a rule for that branch after seeing the numbers is exactly the move
this repo's posture exists to prevent, so the script reports the state and names
the decision instead of manufacturing it.

### What the rules DO fix

1. `ox_rel` and `ox_ppd` are **separate instruments** (R1) - neither may absorb
   the other, and a result on one is not a result on the other.
2. `ox_rel` is **not** the least MYC-entangled ruler in human on two of the
   mouse's own three estimators (R2), so the mouse's argument for promoting it
   does not transfer.
3. The configuration **survives every ruler tried** (R3, R4), so nothing in the
   synthesis document needs restating on ruler grounds.

### The recommendation, labelled as a recommendation

**Keep `OXPHOS subunits` GSVA (`ox_gsva`) as the primary human ruler, and add
`ox_rel` to the reported panel as a named third instrument.** Three reasons:

- **Nothing licenses a promotion.** R2 failed. Promoting `ox_rel` on the
  strength of the one estimator where it does order as the mouse predicts would
  be picking the cell that agrees.
- **Nothing requires a change.** Every number in
  `docs/2026-09-03_human_arm_for_mouse_reconciliation.md` is on `ox_gsva`, and
  R4 - the number that carries the argument - holds on it as well as on
  `ox_rel`. A ruler swap would cost a full re-derivation and buy no claim.
- **`ox_rel` earns a place in the panel anyway.** It is genuinely distinct (R1),
  it is the least proliferation-entangled ruler in both cohorts, it is the most
  infiltrate-robust in TCGA, and it is less entangled with every *signature*
  MYC estimator than the incumbent. CLAUDE.md trap 5 says report all four
  instruments or justify the one; `ox_rel` makes that five, and its
  disagreements are informative rather than noise.

### What this means for G3, which is the question G1 was asked to serve

**The ruler is not the thing to choose. The ruler-estimator pair is.** In the
mouse, `ox_rel` was less MYC-entangled than the alternatives on all three
estimators, so "pick `ox_rel`" was a complete instruction. In human it is less
entangled on the MSigDB signature, on `FELSHER` and on the low-entanglement
signature, and **more** entangled on the CollecTRI regulon - in both cohorts.
Any G3 that fixes a ruler without fixing an estimator has fixed half a
measurement, and the half it left open is the one CLAUDE.md trap 3 already says
is worth a thirty-fold spread.

If G3 wants the most separable available pair on these numbers, it is
**`ox_rel` with `FELSHER__MITOSTRIP`** (0.259 / 0.328) or **`ox_rel` with
`MYC_UP.V1_UP__MITOSTRIP`** (0.329 / 0.338) - and that choice should be declared
before the model is fitted, not after. Stating it here does not pre-register it;
nothing in this repo is pre-registered.

---

## 9. What would change these answers

1. **A third cohort.** R2's failure rests on two cohorts that agree with each
   other; the mouse ordering held on three estimators there and on one here.
2. **A different denominator.** `rest_mito` is the whole MitoCarta inventory
   minus the numerator. A denominator restricted to, say, the mitoribosome would
   be a different instrument with a different entanglement, and the mouse's
   `relify()` helper builds exactly those. Not asked here.
3. **`COLLECTRI_MYC_STIM`**, 739 genes, in the snapshot and never scored
   (synthesis section 10.5). The regulon estimator is the one where `ox_rel`
   loses, and there is a second regulon estimator sitting unused. **This is the
   cheapest thing that could overturn section 3.1's pattern**, and it needs a
   pipeline re-run.
4. **SCAN-B purity.** Section 6 is TCGA only and structurally cannot replicate.

## 10. What was NOT done, listed so it stays a decision

No interaction model. No MYC stratification. No Johnson-Neyman. No subtype
split. No new endpoints or ratios. No mitoPPS recomputation. Nothing written to
`myc_human_validation` or to `myc_mouse`. No ortholog function called, in either
direction - `grep -rnE "(mouse_to_human|human_to_mouse|ortholog[s]?)[[:space:]]*\("
scripts/` returns nothing.

---

## 11. Where the numbers live

| | |
|---|---|
| script | `scripts/E16_respiratory_rulers.R` |
| object | `results/respiratory_rulers.rds` - `$check1`, `$separability`, `$twelve`, `$repro`, `$ruler_stability`, `$the_gap`, `$infiltrate`, `$infiltrate_shift`, `$r1_delta`, `$r2`, `$r3`, `$r4`, `$r5`, `$verdicts`, `$primary_ruler`, `$set_sizes`, `$overlap_audit`, `$ruler_sensitivity`, `$sign_flips` |
| tables | `outputs/tables/E16_ruler_separability.csv`, `outputs/tables/E16_twelve_on_rulers.csv` |
| figures | none - every number here is a table |
| mouse source | `myc_mouse` `e348dd8` (`paper-final`), `scripts/48_gate_model_mouse_verification.R` PART A, and `docs/2026-09-02_myc_oxphos_priming_gate_model.md` section 2, both read read-only |
