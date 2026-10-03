---
date: 2026-10-03
status: >
  EXPLORATORY LEADS. Hypothesis-generating for future work. NOT part of the
  declared E39/E40/E41/E42 line and NOT reported in the manuscript.
posture: >
  No result in this note or in scripts/X01_lead_subtype_oxphos.R may enter any
  manuscript document, figure, legend or supplementary file without its own
  declaration written and committed first. It carries no reading rule, no
  pass, no fail and no falsification criterion. It is a survey of where a
  predictive-marker hypothesis might be worth building, and nothing more.
why-not-an-E-number: >
  The E-numbered line is manuscript-bound: every E-script's result is reported
  whatever it shows, under a declaration committed before it ran. X01 is not.
  Keeping it off the numbering makes the difference visible in the file listing
  so nobody downstream mistakes it for part of the declared work.
relates-to:
  - 2026-10-02_E39_respiratory_axis_decomposition_declaration.md (the declared
    line this is NOT part of; sections 13.2 and 13.4 in particular)
  - results/e41_metabric_ladder.rds (the ER-positive estimate this surveys)
  - results/e42_pcr_subtype_scope.rds (the subtype pattern that points here)
---

# X01 - where might a predictive marker around OXPHOS be?

## 1. What prompted this, and what it is not

Two marginal results point at the same stratum.

**E41, METABRIC.** In ER-positive disease, proliferation-adjusted, higher
respiration tracks lower breast-cancer mortality: HR 0.896 [0.803, 1.000],
1,505 patients, 431 events. The upper bound rounds to the null and the declared
ER.Expr sensitivity gives 0.929 [0.833, 1.036].

**E42, the pCR scope check.** The HR-positive/HER2-negative stratum carries the
largest pooled estimate of the three fitted: -0.276 [-0.530, -0.022] across two
cohorts, against E40's pooled -0.162 and TNBC's random-effects -0.110.

So the same stratum shows both arms - worse chemotherapy response, better
long-term outcome - while the receptor-negative strata are null or incoherent on
both.

**That observation is not a finding.** E42 section 13.2 forbids any subtype
result becoming a finding in its own right, and that prohibition is not relaxed
by moving the work into this note. What the convergence does is tell this survey
where to look.

## 2. The three cautions that govern reading anything here

**Caution 1 - winner's curse.** With a pooled METABRIC estimate of HR 0.896,
individual PAM50 subtypes will scatter by sampling variation alone. **The
largest stratum estimate will look impressive and will be partly the maximum of
several noisy draws.** No subtype is the lead merely for being the largest.

**Caution 2 - the luminal pCR estimate is probably inflated.** E42 declared the
detectable effect for HRpos_HER2neg at 0.312 log-odds (0.441 at twofold variance
inflation). The observed pooled estimate, 0.276, is **smaller than the effect
the stratum was declared able to detect**. An underpowered analysis that clears
significance typically overstates magnitude. The true luminal effect is likely
nearer E40's pooled -0.162 than -0.276.

**Caution 3 - the s508 explanation already failed once.** E42 section 13.7
measured rho(OX, subtype2) in GSE25066 at -0.054, refuting the proposed
explanation for that cohort's DRFS divergence. The plausible story was wrong.
Treat every mechanism offered in this note the same way.

## 3. What would make a lead worth building a proposal on

Not the point estimate. Three things, declared here so they are not chosen after
the numbers:

1. **Mechanistic coherence.** Does the subtype where it is strongest make
   biological sense? LumA is the most differentiated and most oxidative;
   Her2-enriched has a distinct metabolic profile. A result that fits the
   biology is worth more than a larger one that does not.
2. **Agreement across the two ER calls.** METABRIC has two partitions,
   ER_IHC_status and ER.Expr, disagreeing on 127 of 1,940 patients. If the same
   subtype leads under both, that is a cheap internal check. It is **not
   independent replication** and must never be described as such.
3. **A proliferation contrast.** The point is that this is proliferation-
   adjusted. A subtype where m0 is null and m1 is protective is a better lead
   than one where both are protective, because the first shows the adjustment
   doing work.

## 4. The treatment question, which may matter more than the subtype split

In METABRIC, 1,111 of 1,505 ER-positive patients (73.8%) received endocrine
therapy and 147 (9.8%) chemotherapy; of 435 ER-negative, 270 (62.1%) received
chemotherapy.

E40 measured **response to neoadjuvant chemotherapy**. E41 measures **prognosis
under mixed, largely historical care**. The two arms point opposite ways, and in
METABRIC they are partly different patients.

**If respiration predicts outcome in opposite directions depending on treatment,
that is a more fundable hypothesis than a single prognostic association**,
because it is actionable: a marker that says which treatment, not only how bad
the tumour is.

The cells will be small and the treatment variable has eight levels. This is a
survey of whether the pattern is there at all, not a test of it.

## 5. What this note cannot license, ever

- No manuscript sentence, figure, legend or supplementary table.
- No claim that any subtype or treatment stratum "shows" or "demonstrates"
  anything. The vocabulary here is *suggests*, *is consistent with*, *would be
  worth testing*.
- No comparison of any coefficient across cohorts as a raw value (E39 section
  4.2 still holds).
- No reinstatement of anything on the retired list.
- No interaction test, no formal effect-modification test. Not because of
  posture but because the event counts do not support one.
