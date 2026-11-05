---
date: 2026-11-05
status: >
  The author's ruling that the MYC/BCL2L1/MCL1 co-amplification result
  (validation study G2) becomes a main panel of Fig. 4 (Fig. 4x), sized like
  the 26-mm correlation panel, and the four choices that specify it. It
  AMENDS decision 5 of 2026-11-02, which admitted the result as an Extended
  Data panel. Recorded verbatim and committed alone, before anything is built.
  It is E54's specification.
relates-to:
  - docs/2026-11-02_separation_line_brief.md (section 2, decision 5; section 9, the scope and the checks before any panel)
  - docs/2026-11-05_E52_prism_size_ruling.md (9c2ca44; the 26-mm panel this one is sized like)
  - myc_human_validation d3ac60e, docs/2026-08-28_G2_result.md and docs/2026-08-28_D4_scoping_and_G2_design.md (read only)
  - docs/handoff-2026-10-21_validation_audit.md (G2's pre-registration record)
---

# The co-amplification panel becomes Fig. 4x (the author, 2026-11-05)

## 1. What this is

**Sections 2 and 3 are the author's, recorded unchanged.**
- **Section 2** is the author's message.
- **Section 3** gives the author's answers to four questions Claude Code asked about it. Each question and the option the author chose are quoted as they were shown, with the preview where there was one.

**What it changes:**
- **Decision 5 is amended.** The co-amplification result enters the paper as one main-text sentence and **one main panel, Fig. 4x**, in place of one Extended Data panel (question 1). Nothing else from the validation study is admitted.
- **The panel** is drawn by a new display script, `E54`, with a 26-mm axis, like E52's panel. It shows the adjusted odds ratios for BCL2L1 gain and MCL1 amplification with MYC amplification:
  - first in all tumours, under both pre-specified covariates;
  - then in Basal, Her2, LumA and LumB (question 2).
- **Its numbers come from G2's saved results,** snapshotted into this repo with a provenance record. The script checks every number it draws against the committed G2 result note (question 3).
- **The sentence waits on a later text ruling** (question 4). The caveats of the 2026-11-02 brief, section 9, travel with it.

**Nothing in sections 2 and 3 is edited here.**

## 2. The author's message, verbatim

> OK. back to the main. I have inserted this sentence from a previous version, according to decision 5. `In our analysis, MYC-amplified tumours showed
>   enrichment of BCL2L1 co-amplification across breast cancer subtypes, whereas MCL1 co-amplification was restricted to the basal subtype (Fig. 4x).` Can
>   you find the previous panel on this and adjust it to fit in a similar size as the previous correlation plot in the main figure (26mm axis length).

## 3. The author's answers, verbatim

**Question 1:**

> Decision 5 (your ruling of 2026-11-02) admitted the MYC/BCL2L1/MCL1 co-amplification result as one main-text sentence and one Extended Data panel, 'not as a main panel'. Your sentence cites 'Fig. 4x' and you asked for the main figure's 26-mm size. Where should the co-amplification panel go?

**The author's choice:**

> Main panel, Fig. 4x: Amends decision 5: the panel joins Fig. 4, sized like the 26-mm correlation panel (E52). Recorded as a new ruling before anything is built.

**Question 2:**

> The only earlier co-amplification figure is the validation study's G2 gate diagnostic: four facets, about 35 odds ratios (primary tests, threshold grid, PAM50 strata, source and covariate checks). A 26-mm panel holds about ten rows. Which rows should the new panel show?

**The author's choice:**

> Pooled, then subtypes (Recommended): For BCL2L1 gain and MCL1 amplification: all tumours (the pre-registered primary test, with both pre-specified covariates), then Basal, Her2, LumA and LumB (declared secondary in G2). Matches the sentence's subject. Note that it shows BCL2L1 near 1 in Her2 and LumB.

**The preview shown with it:**

```
Adjusted odds ratio, MYC amp vs not (log)
                 0.3     1      3     10
BCL2L1 gain       .      |      .      .
  All tumours     .      |  *-o-*      .   1.77; FGA 1.36
  Basal           .      |  *----o-----.   2.69
  Her2           *------o|------*      .   0.99
  LumA            .     *|---o----*    .   1.70
  LumB            .  *---o---*     .   .   0.90
MCL1 amp          .      |      .      .
  All tumours     .      | *-o-*       .   1.69; FGA 1.60
  Basal           .      |   *----o----*   3.12
  Her2          *--------|o--------*   .   1.22
  LumA          *-----o--|*        .   .   0.68
  LumB            .  *---|o----*   .   .   1.12

(o = estimate, * = 95% CI ends; Normal stratum,
 n = 35, left out as artefact-prone)
```

**Question 3:**

> G2's saved results (its odds ratios, strata and tests) are gitignored in the validation repo, so they are not in its frozen commit d3ac60e. They sit in its working tree, written by its pre-registered script on 2026-08-28, and their numbers match the committed G2 result note. Where should the new panel take its numbers from?

**The author's choice:**

> Snapshot G2's results (Recommended): Copy G2's saved results into this repo's data/ with a provenance README and md5; the display script checks every number it draws against the committed G2 note. The panel shows the registered result exactly (n = 1,043).

**Question 4:**

> Your sentence says MYC-amplified tumours 'showed enrichment of BCL2L1 co-amplification across breast cancer subtypes, whereas MCL1 co-amplification was restricted to the basal subtype'. In G2: BCL2L1's pooled OR is 1.77 (1.26-2.47) adjusted for aneuploidy but 1.36 (0.98-1.88) for fraction genome altered, both pre-specified; by subtype only Basal's interval excludes 1, and 'across subtypes' rests on a homogeneity test that did not reject (p = 0.33). MCL1 is enriched only at high-level amplification and inverts at any gain (OR 0.50). The subtype tests were declared secondary and exploratory. Your brief of 2026-11-02 (section 9) said to state MCL1's inversion or drop that clause. How should the sentence be handled?

**The author's choice:**

> Rule on it later (Recommended): Build the panel now; the sentence is settled with the text ruling, which can then match what the panel shows.
