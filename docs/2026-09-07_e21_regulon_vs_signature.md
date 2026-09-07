# E21 - why the regulon estimator dissents, and it is not about MYC

**2026-09-07. EXPLORATORY. Nothing here is pre-registered.**

Script: `scripts/E21_regulon_vs_signature.R`. **NOT YET SOURCED BY THE AUTHOR.**
Every number below is from a dry run (script copied to the scratchpad with
`PATH_E21` and the table path redirected; repo `results/` and `outputs/`
untouched, `git status` clean). **The dry-run caveat is discharged only when the
author's run is compared object by object.**

**This is a NEW analysis, not a re-read.** It builds two new instruments over
six gene sets and a size-matched random null. Three reproduction controls tie it
to what already exists, and all three hold exactly.

**N3 throughout.** Every number here is a transcript correlation. The mouse is
not read, compared or cited by this script.

---

## 0. The answer in six lines

| | |
|---|---|
| **Does scoring `COLLECTRI_MYC_STIM` resolve E20's dissent?** | **NO.** The activation-only half gives `d = -0.086 / -0.110` - still negative, both cohorts, material. It is 736 genes, the same size as the whole regulon, and it behaves the same way. |
| **Is it the instrument (ULM vs GSVA)?** | **No.** All three signatures stay positive under ULM (+0.064 to +0.152), 6 of 6 cells. ULM is not what makes `M_b` dissent. |
| **Is it the sign?** | **No.** Splitting the regulon by its own sign and scoring each half unsigned leaves both halves negative. |
| **Is it the SIZE?** | **YES, for the two large sets, and completely.** The size-matched random null's median `d` is **negative at every size** and grows with size (-0.045 at n=61 to -0.089 at n=811). `regulon_all` and `regulon_stim` sit **inside** their own null's 95% interval. **`M_b`'s dissent is what a random 811-gene composite does.** |
| **Does that damage E20's signature result?** | **No - it strengthens it.** All three signatures are **outside** their size-matched nulls, positive, in both cohorts: 6 of 6. |
| **Is anything in the regulon real?** | **Yes, the repressed half.** 72 genes - signature-sized - at `d = -0.171 / -0.177`, far outside its null `[-0.130, +0.048]`. MYC's **repressed** targets genuinely track OXPHOS share rather than content. |

**So the estimator that "disagreed about MYC" was not disagreeing about MYC. It
was too big.** E20's verdict of ESTIMATOR-DEPENDENT stands exactly as scored -
no rule is rewritten - but the dissent that produced it is now explained, and it
is an artefact of composite size measured against a within-sample contrast
ruler.

---

## 1. What `M_b` actually is, and the four ways it differs

This is the part the earlier notes never separated.

| | the regulon estimator `M_b` | the three signature estimators |
|---|---|---|
| **gene set** | MYC's **CollecTRI** regulon: curated TF-target edges from OmniPath (Muller-Dott et al. 2023, `data/collectri_human/`, retrieved 2026-08-28). 886 targets, **811** after the MitoCarta strip | experimentally derived co-expression sets: `HALLMARK_MYC_TARGETS_V1` 177, `MYC_UP.V1_UP` 178, `FELSHER` 61, all `__MITOSTRIP` |
| **instrument** | **decoupleR `run_ulm`** - a per-sample linear regression of the whole expression vector on a target indicator | **GSVA**, rank-based enrichment on the VST, `kcdf = "Gaussian"` |
| **sign** | **signed**: `mor = +1` on a stimulatory edge, `-1` on an inhibitory one, and an edge flagged **both** takes `+1` (E02's rule, reproduced not improved) | unsigned. No GSVA set carries a direction |
| **size** | **811 genes, 4.5% of the matrix** | **61 to 178 genes** |

Four differences at once. Nobody had separated them, and until this script the
phrase "the regulon disagrees" quietly meant all four.

**One thing was already known** (`E06`, `docs/2026-09-01_phase2_estimator_findings.md`
E3): `M_b` is empirically its own **activated** half. It correlates 0.955 / 0.934
with the 736 activated targets scored unsigned, and 0.068 / 0.086 with the 72
repressed ones. The activated targets outnumber the repressed **10.2 to 1**, so
the signed fit is dominated by them and the sign barely operates. That is why
"score the activation-only set" looked like a promising test - and why, in
hindsight, it could not have been one.

### What each estimator shares with what it is measured against

| set | n | x numerator | x rest-of-MitoCarta | x mitoribosome | x `PROLIF_DISJOINT` |
|---|---|---|---|---|---|
| `regulon_all` | 811 | 0 | 0 | 0 | 75 |
| `regulon_stim` | 736 | 0 | 0 | 0 | 73 |
| `regulon_repr` | 72 | 0 | 0 | 0 | 2 |
| `sig_msigdb` | 177 | 0 | 0 | 0 | 41 |
| `sig_felsher` | 61 | 0 | 0 | 0 | 0 |
| `sig_lowent` | 178 | 0 | 0 | 0 | 2 |

Every mitochondrial column is **0 by construction** - `MITOCARTA_ALL` is the
strip set and both halves of every ruler live inside it. The proliferation
column is not zero, which is why the adjusted panel is the one read.

---

## 2. The design: vary one thing at a time

**The statistic is E20's**, one number per estimator per cohort:

```
d = rho(est, ox_lvl) - rho(est, ox_rel)
    positive = the signature behaviour ("content")
    negative = the M_b dissent
```

**The design.** Score every gene set on **one common instrument** so the gene
set varies and nothing else; then apply the **other** instrument to the
signature sets so the instrument varies and nothing else.

| gene set | `zmean` | `ulm` | `gsva` |
|---|---|---|---|
| `regulon_all` (811) | built | built (signed) = **`M_b`** | - |
| `regulon_stim` (736) | built | built | - |
| `regulon_repr` (72) | built | built | - |
| `regulon_diff` (stim - repr) | built | - | - |
| `sig_msigdb` (177) | built | built | read |
| `sig_felsher` (61) | built | built | read |
| `sig_lowent` (178) | built | built | read |

`zmean` is the common instrument - mean of per-gene z on `log2(linear + 1)` -
because it is exactly what the four rulers are built from, and what `E06`
already used to split this regulon.

**GSVA is not run here**, and section 5 tests whether that was allowed rather
than assuming it.

### The three controls, and all three hold exactly

| control | result |
|---|---|
| the rebuilt **signed** ULM must reproduce the saved `M_b__MITOSTRIP` | **max abs diff = 0**, both cohorts. The ULM leg is `M_b`, not a lookalike |
| the two regulon halves must reproduce `E06`'s `$halves` | **max abs diff = 0**, both cohorts |
| the three GSVA signatures and the saved `M_b`, on the four rulers, must reproduce `E20` | **64 of 64 cells, max \|delta\| = 0** |

---

## 3. The panel - `d` for every set on every instrument

Proliferation-adjusted, with the paired sample-level bootstrap (B = 2,000, seed
20260831). **Every cell below is material** (`|d| >= 0.05` and an interval
excluding zero).

| set | instrument | TCGA `d` | SCAN-B `d` | says |
|---|---|---|---|---|
| `regulon_all` | **ulm (= `M_b`)** | **-0.089** [-0.109, -0.069] | **-0.098** [-0.112, -0.084] | share |
| `regulon_all` | zmean | -0.099 [-0.121, -0.078] | -0.121 [-0.135, -0.107] | share |
| **`regulon_stim`** | **zmean** | **-0.086** [-0.108, -0.065] | **-0.110** [-0.124, -0.096] | **share** |
| `regulon_stim` | ulm | -0.105 [-0.125, -0.085] | -0.117 [-0.131, -0.103] | share |
| **`regulon_repr`** | **zmean** | **-0.171** [-0.193, -0.150] | **-0.177** [-0.191, -0.163] | **share** |
| `regulon_repr` | ulm | -0.164 [-0.185, -0.143] | -0.164 [-0.178, -0.151] | share |
| `regulon_diff` | zmean | +0.119 [+0.097, +0.139] | +0.107 [+0.094, +0.120] | content |
| `sig_msigdb` | gsva | +0.083 [+0.065, +0.102] | +0.089 [+0.075, +0.102] | content |
| `sig_msigdb` | **ulm** | **+0.064** [+0.046, +0.084] | **+0.107** [+0.092, +0.121] | content |
| `sig_msigdb` | zmean | +0.066 [+0.048, +0.084] | +0.078 [+0.065, +0.091] | content |
| `sig_felsher` | gsva | +0.109 [+0.088, +0.132] | +0.091 [+0.077, +0.104] | content |
| `sig_felsher` | **ulm** | **+0.088** [+0.067, +0.111] | **+0.098** [+0.084, +0.111] | content |
| `sig_felsher` | zmean | +0.092 [+0.071, +0.115] | +0.084 [+0.070, +0.097] | content |
| `sig_lowent` | gsva | +0.142 [+0.122, +0.165] | +0.116 [+0.102, +0.130] | content |
| `sig_lowent` | **ulm** | **+0.147** [+0.127, +0.170] | **+0.152** [+0.138, +0.166] | content |
| `sig_lowent` | zmean | +0.108 [+0.087, +0.130] | +0.089 [+0.075, +0.102] | content |

Three readings, in order of how much they matter.

**The instrument is not the cause.** Put the signature genes through ULM - the
regulon's own instrument - and all three stay positive, in both cohorts, 6 of 6.
`sig_lowent` is actually *more* positive under ULM (+0.147 / +0.152) than under
GSVA. Whatever makes `M_b` dissent, it is not `run_ulm`.

**The sign is not the cause, and this is the direct answer to "would
`COLLECTRI_MYC_STIM` resolve it".** The activation-only half gives `-0.086 /
-0.110`. It does not flip. Nor does the repressed half, which goes the same way
harder. There is no sign split that rescues the regulon into agreeing with the
signatures.

**Every regulon set is negative on the common instrument and every signature set
is positive on it.** Holding the instrument fixed, the split survives. So on the
rule fixed before the numbers existed, the cause is attributed to the **gene
set** - and section 4 is what says which property of the gene set.

---

## 4. The size null, and it is the answer

`ox_lvl` is a mean-z composite. `ox_rel` is that same composite **minus a
1,047-gene mean-z**, so it cancels whatever two large composites share - and
what a large *random* gene set shares with another large composite is a global
expression component with nothing to do with MYC.

**So a random gene set may produce a non-zero `d` on its own.** 200 random sets
per observed size, drawn from matrix genes that are not MitoCarta members (so a
random set cannot touch either half of any ruler), scored by `zmean`:

**Null distribution of `d`, adjusted:**

| size | TCGA median | TCGA 95% | SCAN-B median | SCAN-B 95% | frac > 0 (T / S) |
|---|---|---|---|---|---|
| 61 | -0.055 | [-0.127, +0.035] | -0.045 | [-0.129, +0.071] | 0.105 / 0.195 |
| 72 | -0.056 | [-0.130, +0.048] | -0.047 | [-0.144, +0.057] | 0.100 / 0.150 |
| 177 | -0.074 | [-0.129, -0.000] | -0.064 | [-0.133, +0.025] | 0.025 / 0.085 |
| 178 | -0.067 | [-0.126, +0.013] | -0.067 | [-0.144, +0.014] | 0.040 / 0.045 |
| 736 | -0.089 | [-0.123, -0.047] | -0.085 | [-0.129, -0.041] | 0.000 / 0.000 |
| 811 | -0.089 | [-0.119, -0.043] | -0.085 | [-0.124, -0.041] | 0.000 / 0.000 |

**The null is negative at every size, and it gets more negative as sets get
bigger.** At 736 and 811 genes, **not one random set in 200 produced a positive
`d`, in either cohort.** The raw null is the same shape.

**Every observed estimator against its own size-matched null, adjusted:**

| cohort | set | size | observed `d` | null median | null 95% | outside null? |
|---|---|---|---|---|---|---|
| TCGA | `sig_felsher` | 61 | **+0.092** | -0.055 | [-0.127, +0.035] | **yes** |
| TCGA | `regulon_repr` | 72 | **-0.171** | -0.056 | [-0.130, +0.048] | **yes** |
| TCGA | `sig_msigdb` | 177 | **+0.066** | -0.074 | [-0.129, -0.000] | **yes** |
| TCGA | `sig_lowent` | 178 | **+0.108** | -0.067 | [-0.126, +0.013] | **yes** |
| TCGA | `regulon_stim` | 736 | -0.086 | -0.089 | [-0.123, -0.047] | **NO** |
| TCGA | `regulon_all` | 811 | -0.099 | -0.089 | [-0.119, -0.043] | **NO** |
| SCAN-B | `sig_felsher` | 61 | **+0.084** | -0.045 | [-0.129, +0.071] | **yes** |
| SCAN-B | `regulon_repr` | 72 | **-0.177** | -0.047 | [-0.144, +0.057] | **yes** |
| SCAN-B | `sig_msigdb` | 177 | **+0.078** | -0.064 | [-0.133, +0.025] | **yes** |
| SCAN-B | `sig_lowent` | 178 | **+0.089** | -0.067 | [-0.144, +0.014] | **yes** |
| SCAN-B | `regulon_stim` | 736 | -0.110 | -0.085 | [-0.129, -0.041] | **NO** |
| SCAN-B | `regulon_all` | 811 | -0.121 | -0.085 | [-0.124, -0.041] | **NO** |

**`regulon_all` and `regulon_stim` sit inside their own size-matched null in
both cohorts.** `M_b`'s dissent is not a statement about MYC's regulon. It is
what a random gene set of that size does against these two rulers.

**And the three signatures are outside their nulls, positive, 6 of 6.** Their
`d` is not a size artefact. E20's content reading, for the estimators that
carry it, survives its own hardest control.

### The one regulon object that is real

**`regulon_repr` - MYC's 72 repressed targets - is signature-sized and sits far
outside its null**, at `-0.171 / -0.177` against a null of `-0.056 / -0.047` and
an interval that comfortably covers zero. It is the only regulon set that says
something the null does not.

That is `E3`'s finding arriving on a new axis. `E3` established that the
MYC-OXPHOS signal in this regulon is carried by the repressed targets going
**down** (`rho` -0.518 / -0.578 against OXPHOS) far more than by the activated
targets going up (+0.246 / +0.100). E21 adds the share-versus-content reading:
**those repressed targets track OXPHOS share, not OXPHOS content**, and they do
it materially and outside any size explanation.

### The balanced contrast, reported and not leaned on

`regulon_diff` = `zmean(stim) - zmean(repr)`, `E3`'s preferred estimator, gives
`d = +0.119 / +0.107` - it sides with the signatures. **But it is itself a
within-sample contrast of two composites, structurally the same kind of object
as `ox_rel`, and it has no size-matched null** because a difference of two sets
has no single size. It is reported because it is informative and it is not part
of any verdict here.

---

## 5. Was `zmean` allowed to stand in for GSVA?

The gene-set leg rests on `zmean` and no GSVA was run. That is licensed only if
the two agree where both exist - the three signature sets. Measured:

| cohort | set | zmean `d` | gsva `d` | gap |
|---|---|---|---|---|
| TCGA | `sig_msigdb` | +0.066 | +0.083 | -0.017 |
| TCGA | `sig_felsher` | +0.092 | +0.109 | -0.017 |
| TCGA | `sig_lowent` | +0.108 | +0.142 | -0.034 |
| SCAN-B | `sig_msigdb` | +0.078 | +0.089 | -0.011 |
| SCAN-B | `sig_felsher` | +0.084 | +0.091 | -0.006 |
| SCAN-B | `sig_lowent` | +0.089 | +0.116 | -0.028 |

**All six the same sign, largest gap 0.034 against a bar of 0.05. `zmean` is a
fair stand-in.** GSVA is systematically a little more positive than `zmean` on
every one of the six, which is worth knowing and does not change a sign.

---

## 6. Attribution, on the rule fixed before any number existed

| cause | fires |
|---|---|
| **SIGN** - `COLLECTRI_MYC_STIM` resolves it | **NO** |
| **GENE SET** - the split survives holding the instrument fixed | **YES** |
| **INSTRUMENT** - the signatures go negative under ULM | **NO** |
| **SIZE** - the null already gives the *signatures'* own `d` | **NO** |

**One honest note about the rules.** The `SIZE` row was written to ask whether
the **signatures'** `d` is explained by size; it is not, so it does not fire.
The finding that the **regulon's** `d` *is* explained by size comes from the
per-estimator null table in section 4, not from that rule. The pre-fixed rule
did not anticipate the shape the answer took, and saying so is cheaper than
rewriting the rule afterwards.

Read together: **`GENE SET` fires, and the null says which property of the gene
set it is - the size.** The regulon is not a different view of MYC. It is a
composite four to thirteen times larger than any signature, and at that size the
statistic has a built-in negative bias that swamps whatever MYC signal the set
carries.

---

## 7. What this does to E20, and to what I said about it

**E20's verdict is unchanged.** ESTIMATOR-DEPENDENT was scored on a rule fixed
in advance, the regulon did dissent, and no rule is being rewritten now.
`docs/2026-09-07_e20_share_vs_content.md` stands.

**Three statements around it do change.**

1. **"`COLLECTRI_MYC_STIM` is the cheapest thing that could settle this" was
   wrong**, and E21 is the test that shows it. That set is 736 genes - the same
   size class as the regulon it came from - so it inherits the same size bias
   and returns the same null-consistent answer. E20 section 0 and section 6
   item 1 name it as what would settle the question; it does not. The reason it
   looked promising was `E3`'s sign finding, and the sign turns out not to be
   the operative difference.
2. **E20 section 6 item 4 asked for a technical control on `ox_nuc_mtrib`** -
   whether a difference of two correlated composites carries a real signal.
   E21 answers the same question for `ox_rel` and the answer is that **the
   contrast rulers have a size-dependent null that must be checked before any
   `d` is read.** That is now measured rather than flagged.
3. **The regulon's standing as an independent view of MYC is weaker than the
   panel implies.** `M_b` remains a legitimate estimator - E16, E20 and the
   synthesis all report it correctly - but on *this* statistic it carries no
   information beyond its size, and that should travel with it.

**What does not change:** `M_b`'s behaviour elsewhere in the repo. Every other
result that uses it - E16's separability table, E10's twelve genes, the
synthesis's estimator panel - is a correlation with `M_b` itself, not a
difference between two contrast rulers, and nothing here bears on those.

---

## 8. What would change these answers

1. **A regulon of signature size.** Every conclusion about the regulon here is
   confounded with its size, and the snapshot already carries the column that
   would break the confound: `curation_effort`, never used anywhere in this
   repo. Of MYC's 886 targets, **313 clear `curation_effort >= 5` and 101 clear
   `> 10`** - so a high-confidence regulon in the signatures' own size range is
   directly available, needs no new data, and would separate "regulon" from
   "large" for the first time. **This, not `COLLECTRI_MYC_STIM`, is the cheap
   test that was wanted.** Its falsifier is already implied: if a 100-300 gene
   MYC regulon still gives a negative `d` outside its size-matched null, the
   dissent is about regulons after all and section 4's reading is wrong.
2. **A GSVA run over the three regulon sets.** Section 5 licenses `zmean` as a
   stand-in on the signatures and cannot license it on sets four times larger.
   One batch with the `.PIN_A` / `.PIN_B` pins would close it.
3. **The null's mechanism.** This script establishes the size-dependent negative
   bias empirically and does not explain it. The obvious candidate is that a
   large mean-z composite is dominated by a global expression axis that `ox_rel`
   cancels by construction, but that is an interpretation and it is untested
   here.
4. **`regulon_repr` at n = 72 is small**, and its own null interval is
   correspondingly wide. It clears it comfortably in both cohorts, but 72 genes
   is 72 genes.

---

## 9. What was NOT done, listed so it stays a decision

No GSVA run. No interaction model, no MYC stratification, no subtype split, no
new endpoint or ratio. The sign rule was **reproduced, not improved** - the 86
both-flagged edges still count as activating, because changing that would break
the `M_b` reproduction control. No score was pooled or compared across cohorts.
The mouse was not read, compared or cited. Nothing written to
`myc_human_validation` or to `myc_mouse`. No ortholog function called, in either
direction.

---

## 10. Where the numbers live

| | |
|---|---|
| script | `scripts/E21_regulon_vs_signature.R` |
| object | `results/regulon_vs_signature.rds` - `$dstat`, `$panel`, `$vs_null`, `$null_summary`, `$null_d`, `$instr_check`, `$attribution`, `$verdicts`, `$mb_repro`, `$halves_repro`, `$repro20`, `$overlap_audit`, `$settings`, `$rules` |
| table | `outputs/tables/E21_regulon_vs_signature.csv` |
| figures | none |
| the controls | `$mb_repro` and `$halves_repro` at 0; `$repro20` at 64 cells, max \|delta\| = 0 |
| what it explains | `docs/2026-09-07_e20_share_vs_content.md` sections 0, 3 and 6 |
| what it builds on | `docs/2026-09-01_phase2_estimator_findings.md` E3; `E06`'s `$regulon` and `$halves` |
