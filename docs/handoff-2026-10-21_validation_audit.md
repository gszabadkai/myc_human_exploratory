---
date: 2026-10-21
status: >
  READ-ONLY AUDIT of the frozen sibling repository `~/code/myc_human_validation`.
  Nothing in that repository was modified. It was read only through
  `git -C ~/code/myc_human_validation ...` (log, show, grep, ls-tree, diff), with
  no checkout and no R. This file is uncommitted, and it proposes no banners
  and no edits: any correction to the frozen repository is a separate decision
  for the author.
scope: >
  Built only from the validation repository. No other repository's account of
  it was used.
---

# Audit of `myc_human_validation` at `d3ac60e`

**Citation convention.**
- `<file>:L<n>` is a line of that file at HEAD `d3ac60e`, read with
  `git -C ~/code/myc_human_validation show HEAD:<path>`.
- `<file>@<hash>:L<n>` is a line of an earlier version.
- "msg `<hash>`" is that commit's message.
- `plan` is `docs/2026-08-27_human_validation_plan.md` (the build
  specification). Other dated notes are cited by filename without `docs/`.

**Authorship.** Every commit's git author is `gszabadkai`. The messages carry
`Co-Authored-By: Claude Opus 5`, and the notes are written in the assistant's
first person. Where a note credits the author with a decision it is quoted.
Otherwise git cannot say who decided.

---

## 1. Freeze state

| check | result | source |
|---|---|---|
| HEAD | `d3ac60ec06c8e07df0626cbc063a52186ee167c5` on `main` | `git rev-parse HEAD` |
| refs | `refs/heads/main`, `refs/remotes/origin/main` and `origin/HEAD` all at `d3ac60e` | `git for-each-ref` |
| anything above HEAD | **nothing**: `git log --all --not HEAD` is empty; 76 commits on all refs, 76 reachable from HEAD | `git rev-list` |
| tags | **none** | `git tag -l` |
| working tree | clean (read with `--no-optional-locks`, so no index write); no stash | `git status --porcelain` |
| history | linear, no merges; author date equals committer date on the key commits | readers' checks |

**The last commits** (2026-08-31, +0100):

| commit | time | what it did |
|---|---|---|
| `b8ef43f` | 00:47:46 | "Record Block F1: H4 is not supported, and the BIM replication is blocked" |
| `6dc6a47` | 00:48:49 | handoff at Block F1, "with all four hypotheses tested" |
| `a0c5519` | 11:47:46 | "Take option 1: declare the BIM replication for SCAN-B, pre-data" |
| `a5d892c` | 11:59:44 | script 16, fetch SCAN-B (GSE202203) |
| `23c6748`, `5a0bd8b` | 12:09:36, 12:17:58 | symbol harmonisation addendum, then its correction |
| `37816ea` | 12:28:43 | script 17, the BIM replication in SCAN-B |
| `7cc3737` | 13:05:55 | "Record the SCAN-B result: BIM does not replicate, and the sign is reversed" |
| `7707aeb` | 13:13:50 | "Draft the human arm write-up as a well-formed negative" |
| `d3ac60e` | 13:24:15 | "Rewrite the handoff for the write-up phase" |

**What the repo's final handoff says the arm concluded**
(`2026-08-31_handoff.md` at `d3ac60e`):

- "**The arm is finished. There is no remaining analysis.**" (L12)
- "**No hypothesis is supported, and nothing survives replication.**" (L72)
- CLAUDE.md agrees: "Scripts 00–17 are all written and run. All four hypotheses
  are tested and none is supported. The BIM replication ran in SCAN-B and
  FAILED ... The arm is analytically complete" (CLAUDE.md:L32-35), and "Three
  of four falsification criteria are met" (CLAUDE.md:L51).

**The repo does not describe itself as frozen.**
- The handoff says only "Everything is committed and pushed" (L10) and "Do NOT
  re-derive, re-audit or re-run anything" (L25).
- It still expects work: the write-up, and a display-item decision that "gates
  script 18" (L140-147). Script 18 does not exist at HEAD.
- The handoff's own status is "disposable - delete once used" (L3).

---

## 2. The pre-registered hypotheses

All four are in `plan@b542e70`, the plan's first commit (2026-08-27 23:17:23),
which precedes every script (the first script commit is `6382e77`, 2026-08-28
01:17:34). The H1-H3 text and all four falsification bullets are byte-identical
at HEAD. The only change inside the hypothesis section is an H4 paragraph added
in `67f144a` (see H4).

| | prediction (verbatim) | tested in | outcome |
|---|---|---|---|
| **H1** | "`MYC:OXPHOS` interaction on `PRIME` is positive and significant, **and** the MYC-high / OXPHOS-high quadrant is enriched for `BUFFER`" (plan:L79-80) | clause 1: Block C, TCGA n = 938, both instruments. Clause 2: Block B, the escape reading | **NOT SUPPORTED, on both clauses** (see below) |
| **H2** | interaction "null or attenuated overall, but is **restored in the PIK3CA-wild-type / PTEN-intact stratum**, and FOXO3 regulon activity is low" in MYC-high/OXPHOS-high (plan:L92-94) | Block C strata and FOXO3 | **NOT SUPPORTED.** PI3K-intact is "the flattest stratum of all" (`2026-08-29_block_c_result_H1_not_supported.md:L171-172`) |
| **H3** | "The `MYC:OXPHOS` interaction on `PRIME` persists in the **TP53-mutant** stratum." (plan:L104) | Block C strata | **NOT SUPPORTED.** "No CI excludes zero anywhere in that table" (block_c:L169) |
| **H4** | "**OXPHOS-high is protective unless buffered.**" (plan:L124); the unbuffered MYC-high/OXPHOS-high state is "**sensitive** (high pCR, chemo benefit)" (plan:L117-122) | Block F1, three neoadjuvant cohorts | **NOT SUPPORTED, as a null rather than an inversion** (section 4) |

### H1 in numbers

**Clause 1, Block C (`2026-08-29_block_c_result_H1_not_supported.md`):**
- Verdict: "H1 is NOT supported in TCGA-BRCA. H2 and H3 are not supported
  either." (L10)
- `M_a:OXPHOS` on PRIME: GSVA **0.0204** (95% CI -0.0332 to 0.0740, p = 0.455);
  mitoPPS **0.0099** (-0.0452 to 0.0649, p = 0.726) (L99-100).
- All four decision clauses "**FAIL**" (L113-116).

**Clause 2, Block B (`2026-08-29_block_b_result_escape_not_supported.md`):**
- "H1 has now failed on BOTH clauses" and "The escape reading is not
  supported, and fails in the REVERSED direction" (L10-11).
- E1, E2a, E3 and E4 FAIL (L30-34).

**Supporting functional test, Block G
(`2026-08-30_block_g_result_and_bcl2l11_declaration.md`):**
- "G-a and G-b FAIL" (L9).
- In breast cell lines (n = 51) and pan-cancer (n = 1,130, L37), the MCL1 and
  BCL2L1 dependency coefficients include zero (L31-45).
- The plan had said a negative "weakens H1 considerably" (plan:L913-915).

**Forkscale stratum, Block D (`2026-08-30_F3pre_and_block_D_result.md`):**
- "**D-a is FALSE.**" (L63)
- MB2 three-way: GSVA 0.0181, p 0.822; mitoPPS 0.0145, p 0.830 (L59-60).

### Falsification criteria (plan:L147-161, unchanged since `b542e70`)

| criterion | status at HEAD | source |
|---|---|---|
| 1. interaction on PRIME null in the full cohort and in every pre-specified stratum | **MET**, "on its complete stratum list" once Block D added forkscale | plan:L947-951; F3pre_and_block_D:L119 |
| 2. no BUFFER enrichment in the MYC-high/OXPHOS-high quadrant | **MET** | block_b:L144-148 |
| 3. FOXO3 regulon shows no relationship to OXPHOS | **MET** | block_c:L365-370 |
| 4. "**H4 fails in the informative direction**: OXPHOS-high predicts chemo*resistance* regardless of buffering status" (plan:L156) | **NOT MET**: H4 failed as a null | `2026-08-31_block_F1_result_H4_not_supported.md:L12-13` |

The plan's text "If that is the outcome, report it, or drop the human arm. Do
not reach for a fifth post-hoc hypothesis." (plan:L160-161) stands.

---

## 3. The cheap gates

The plan's section 6 is titled "Gates (run first, cheap, they decide what gets
built)" (plan:L488).

| gate | what it checked | rule (first committed) | result | what it decided |
|---|---|---|---|---|
| **G1**, MYC signature overlap audit | Felsher and CollecTRI MYC sets against MitoCarta and proliferation; stripped gene counts | "if the stripped Felsher signature falls below ~50 genes or loses its correlation structure ... the CollecTRI regulon plus the 8q24 CNV instrument carry the MYC axis" (`plan@b542e70:L432-434`) | **count criterion PASSED** (Felsher 67 → **61**; `2026-08-28_G1_result_and_decisions.md:L35`). "**G1 is NOT fully discharged.**" (L26) | Felsher (M-a) stays primary; the fallback was not triggered. **The correlation criterion is recorded as passed only in msg `ccc2ef1`** ("stripped mean rho 0.249 vs null 95th pct 0.013") and in a since-deleted handoff. **plan:L12 still lists it "open".** |
| **G2**, CNV co-occurrence | MYC amplification with MCL1 / BCL2L1 amplification, aneuploidy-adjusted | question and branch: `plan@b542e70:L440-448`, "if yes, H1 becomes Panel b and H4 gains its stratifying variable. If no, H2 moves to primary". Numeric pass rule (adjusted OR > 1, CI excluding 1, in ISAR) first committed in `77968ac` (`2026-08-28_D4_scoping_and_G2_design.md:L324-327`) | "**G2 PASSES.**" (`2026-08-28_G2_result.md:L30`). MCL1 adj OR 1.69 (1.17-2.44); BCL2L1 1.77 (1.26-2.47); BBC3 FAIL (L42-45) | H1 becomes Panel b; "`BCL2L1` carries Panel b" (L35). **Disclosed:** crude ORs were computed before the threshold was fixed (D4_scoping:L301-317) |
| **G3**, forkscale availability | whether the Menegollo forkscale exists for TCGA | `plan@b542e70:L452-458`: if absent, "a failed G3 is a reason to drop ED2" | **PASSES** (`2026-08-29_G3_result_forkscale_availability.md:L10`); TCGA complete at n = 1,037 (plan:L529-537) | no MCbiclust re-derivation; F3 split into F3-pre (TCGA) and F3 proper (METABRIC) |
| **F3-pre** (an exposure-side gate added later) | Spearman of MB1_forkscale against OXPHOS subunits | REDUNDANT if abs(rho) >= 0.70 on both instruments, INDEPENDENT if <= 0.30 on both, else INTERMEDIATE (G3 note §4b) | **INTERMEDIATE**, GSVA 0.529 / mitoPPS 0.418 (F3pre_and_block_D:L39-41). "**The stop gate does not close.**" (L45) | the F3 stop gate stays in force |
| **F3 proper**, the stop gate | incremental value of `STATE` over `MB1_forkscale` (nested Cox) | "If `STATE` adds nothing over `MB1_forkscale` ... **Run this before building Panel c.**" (plan:L876-878); "Do not build 13 before F3 returns." (plan:L1050) | **NEVER RUN**: blocked on a METABRIC identifier file (`2026-08-31_handoff.md:L153-154`) | **script 13 (H4) was built and run anyway.** Msg `259efa9` puts the gate "to the author rather than deciding it ... Reading the gate as blocking F1 would not defer F1, it would end it. That is the author's call". **No document records the author's decision.** |

**One further threshold.** The MYC estimator agreement threshold, rho >= 0.30,
is called "pre-registered before the number was seen"
(`2026-08-28_myc_estimator_validation.md:L56-57`). It is present in script 06 at
`f9d3861` before the result `87d9ad1`, and it passed at 0.613. But a
since-deleted handoff calls it "not formally signed off".

---

## 4. The conditional-chemosensitivity test (H4)

**The label is the repo's own.** The plan heading is "H4 - conditional
chemosensitivity (NEW in v2)" (plan:L110).

**Commit `b8ef43f`.** Confirmed: `b8ef43fc35a26fef584ea3b9b665f0be7ad7bd35`,
2026-08-31 00:47:46 +0100.
- Subject: "Record Block F1: H4 is not supported, and the BIM replication is
  blocked".
- Files: `CLAUDE.md` and the new `2026-08-31_block_F1_result_H4_not_supported.md`.
- **It is the RESULT commit, not the declaration or the script.**
- No file in the repo cites the hash (`git grep b8ef43f HEAD` is empty), so the
  external "project record" that names it cannot be checked from here.

**Declaration chain** (all +0100):

| commit | date and time | what |
|---|---|---|
| `b542e70` | 08-27 23:17:23 | H4 and criterion 4 declared. Primary = "level 3 vs level 4", categorical; continuous "as sensitivity, not as an alternative to be swapped in" (`plan@b542e70:L525-529`) |
| `67f144a` | 08-28 12:25:18 | D5: primary becomes the continuous three-way `MYC:OXPHOS:BUFFER`, with the categorical `STATE` contrast as "pre-specified secondary" (plan:L128-130). **This reverses the original guard, before any neoadjuvant data existed** |
| `f1f0c0a` | 08-30 23:23:06 | `STATE` frozen. `BUFFER_c = mean( z(log2 MCL1), z(log2 BCL2L1) )` declared (`2026-08-30_STATE_frozen_and_H4_buffer_declaration.md:L230`), with `MYC:OXPHOS:BUFFER_c` predicted **NEGATIVE**. The note's own lock: "Nothing in this note may be revised after any of those outcome columns has been read" (L17-18) |
| `4f761ff` … `f2386f5` | 08-30 23:36 to 08-31 00:03 | script 12 fetches the three cohorts, including their pCR columns |
| `f299c59` | 08-31 00:10:27 | amendment A1 (H4 "GSVA alone", `2026-08-31_H4_amendment_instruments_and_coverage.md:L88-90`), A2, C3. **Edits the f1f0c0a note after the pCR columns were read**; it justifies this as before any score-versus-pCR association (amendment:L34-36) |
| `259efa9` | 00:20:21 | script 13, with A4/A5 "**Added 2026-08-31 while writing script 13, before it was run.**" (amendment:L222). New lock: "After script 13 computes its first score-versus-pCR association, nothing in this note or in the 2026-08-30 declaration may be amended." (amendment:L348) |
| `cc34f74` | 00:33:27 | **script 13 changed after its first run**, to add the per-gene companion that declaration §6.3 required, "found only when the first run's limb correlations were inspected" (msg `cc34f74`). No declared direction changed |
| `b8ef43f` | 00:47:46 | result |

**The primary test** (amendment A4, L249): `pCR ~ MYC * OXPHOS * BUFFER_c +
subtype + treatment`.
- Run per cohort: GSE194040 (I-SPY2, primary), GSE164458 and GSE25066, with
  "n = 988 / 482 / 470" (F1 result:L27). The estimates are meta-analysed; the
  scores are never pooled.
- GSVA only.
- The prediction was a NEGATIVE three-way coefficient.

**Primary result** (F1 result:L54-60):

| cohort | three-way coefficient |
|---|---|
| GSE194040 | **-0.0082** |
| GSE164458 | +0.4066 (p 0.0085) |
| GSE25066 | +0.2399 |
| fixed effect | **+0.1629**, p 0.052 |
| random effects | **+0.1952** (-0.069 to 0.459), p 0.148 |

I2 is **58.1%**. Verdict: **FAIL**. The pooled estimate is positive where a
negative was predicted, and the primary cohort sits at zero.

**The pre-declared secondary limbs** (F1 result note):

| limb | declared in | returned | verdict |
|---|---|---|---|
| F1-b, same sign in the primary cohort | script 13 | -0.0082 | "pass", then corrected in §6: "A coefficient indistinguishable from zero cannot corroborate a direction." (L187) |
| F1-c, two-way `MYC:OXPHOS` | declaration §6.5 | meta FE +0.079 (-0.036 to +0.193), p 0.176 | **FAIL** |
| F1-d, specificity (not the assembly factors) | amendment A2 | the control arm sits "at 63% of the primary arm's magnitude" | pass, **vacuously** |
| F1-e, intersection-scoring sensitivity | amendment A2 | meta FE +0.154 | pass (the only genuine one) |
| F1-f, `STATE` level 4 vs 3 (the frozen categorical secondary) | plan §7.5 | "Two of three run **opposite**" (L142) | **FAIL** |
| F1-g, the informative failure (criterion 4) | declaration §6.6 | **-0.132 (-0.362 to +0.097), p 0.259** (L146-147) | **FAIL**: "**This is a null, not an inversion.**" (L148-149) |
| F1-h, per-gene companion | declaration §6.3; added in `cc34f74` | GSE164458's result is MCL1 (+0.4383, p 0.00057) and not BCL2L1 (+0.0099, p 0.926) (L80-82); MCL1 heterogeneous across cohorts (I2 82.2%) | pass, vacuously. MCL1 "NOT PURSUED" |

**Declared but not reported.** The amendment promised cohort-specific covariate
sensitivities (amendment:L253-255), an `mp` sensitivity (L260-263) and
imputation counts. None appears in the result note. **Unverified** whether they
were run.

**Was the frozen inversion criterion met? No.**

The criterion (STATE declaration:L314-316, frozen at `f1f0c0a`):

> "The OXPHOS effect on pCR **within MYC-high** - the model's OXPHOS slope
> evaluated at `MYC = +1 SD` - is **negative with a CI excluding zero**, on both
> instruments, **while `MYC:OXPHOS:BUFFER_c` is null**."

The result: "It fails as a NULL, not in the informative direction. Falsification
criterion 4 is therefore NOT met" (F1 result:L12-13).

**One unaddressed tension.** The criterion requires the slope "on both
instruments". After amendment A1, H4 is GSVA-only. So, as written, the
criterion could not have been met even by a significant negative slope. It was
evaluated on GSVA in the primary cohort alone. No document addresses this.

**"The BIM replication is blocked"** (F1 result §9, L234-258) means the
following:
- The Block C replication declaration required both co-primary instruments.
- mitoPPS does not exist in the three neoadjuvant cohorts.
- The declaration's own text says "if it is not honoured, delete it rather than
  amend it".
- Three options were recorded and none was taken. Section 5 below follows the
  outcome.

---

## 5. The BIM (BCL2L11) replication in SCAN-B

**Where it came from.** In Block C, BIM was a designated endpoint negative.
Its TCGA limb fired: log2(BCL2L11) "+0.051 (p 0.005) +0.068 (p 3e-4)"
(block_c:L226-228).
- Then `d7f4121` (2026-08-29 02:01:01) and `b2f9fb4` (02:01:22) reclassified it,
  about 40 minutes after the Block C result `f2def97` (01:21:48). "`BCL2L11` is
  NOT a negative control and the sentence above is wrong." (plan:L207-225)
- The grounds cited are a mouse-repo pre-specification. **Not verifiable here.**
- The replication commitment was written into Block C note §9 in the result
  commit `f2def97`. It is "Declared **before METABRIC, SCAN-B or any neoadjuvant
  cohort has been opened**" (block_c:L325-326). That is pre-data for the
  replication cohort, but after the TCGA number was seen.

**What was tested** (`2026-08-31_scanb_bim_replication_declaration.md`,
`a0c5519`, 2026-08-31 11:47:46):
- Model: `log2(BCL2L11) ~ M_a * OXPHOS_subunits + PROLIF_DISJOINT + PAM50`
  (L136), on GSVA and mitoPPS.
- Cohort: SCAN-B GSE202203.
- Companions: BCL2L1, predicted negative, and BBC3, predicted null.
- Comparison: the TCGA Block C limb.

**What it said in advance:**
- "`BCL2L11` replicates if and only if:" the interaction is **positive** at
  p < 0.05, on **both** instruments, at S1, and beats its expression-matched
  null (L255-260).
- "One instrument alone is not a replication." (L262)
- On failure: "reported as a TCGA-specific observation and dropped. No further
  variants are tried." (L278-282)
- **"Reversed" was not a pre-declared category.** Only replicates / fails.

**Order of commits:**

| commit | time | what |
|---|---|---|
| `a0c5519` | 11:47:46 | declaration |
| `a5d892c` | 11:59:44 | fetch script |
| `23c6748` | 12:09:36 | symbol addendum, "Added after the three files were downloaded and their symbol column read ... before any expression value entered any model" (L360-361) |
| `5a0bd8b` | 12:17:58 | its correction, "CORRECTED 2026-08-31, after script 16 ran" (L417) |
| `37816ea` | 12:28:43 | script 17 |
| `7cc3737` | 13:05:55 | result |

The declaration also records having read "the first ~3 MB" of the counts file
before declaring, but no expression value for a named gene.

**Result** (`2026-08-31_scanb_bim_replication_result.md`):
- TCGA: **+0.051 (+0.015, +0.086) p 0.005** (L24).
- SCAN-B, n = 3,143: GSVA **-0.036 (-0.066, -0.006) p 0.019**; mitoPPS
  -0.033 p 0.037 (L25).
- Matched null: 10 / 2000 on GSVA (L101), 1 / 2000 on mitoPPS.
- "**`BCL2L11` does not replicate in SCAN-B.**" (L18), and "the effect is there,
  at n = 3,143, with the opposite sign." (L112-113)
- Dropped as declared.
- It does not close criterion 4 (msg `7cc3737`).

---

## 6. Pre-declaration discipline, hypothesis by hypothesis and gate by gate

"Before" below is commit order. **Git records commits, not executions.** That
a model was not fitted before its declaration rests on the documents' own
statements.

| item | declaration (commit, date) | fitting script | before? | caveat |
|---|---|---|---|---|
| H1-H4, falsification criteria | plan `b542e70`, 2026-08-27 | first script `6382e77`, 08-28 | **yes** | - |
| G1 | plan `b542e70` | `6382e77` (script 04) | yes | correlation criterion closed only in a commit message |
| G2 | question `b542e70`; numeric rule `77968ac` (08-28 11:11:54) | script 05 in the same `77968ac` ("written, NOT run") | rule with the script | crude ORs seen before the threshold was fixed (disclosed) |
| G3 | plan `b542e70` | data snapshot `1fb703e` | yes | - |
| MYC estimator threshold rho >= 0.30 | script 06 `f9d3861` | same | in-script | "not formally signed off" (deleted handoff) |
| H1 clause 1, H2, H3 (Block C) | plan `b542e70` + `2026-08-29_script09_build_spec.md` (`48c5ccc`, 08-29 00:26:18) | `0103525` (01:00:00) | **yes** | the spec's status still reads "DRAFT - awaiting sign-off before any code is written" (L3), and **no sign-off is recorded** |
| H1 clause 2 (escape reading, Block B) | `2026-08-29_escape_reading_declaration.md` (`4d68efb`, 02:26:51) | `7ae2a3a` (02:30:22) | **yes** | declared after the Block C result; msg `4d68efb` says Block B's BCL2L1 model "has ALREADY been fitted", so E3 declares no BCL2L1 prediction |
| Block G (DepMap) | prediction in plan `b542e70`; operationalisation and script in `cbeed0d` | `cbeed0d` | prediction yes; operationalisation with the script | - |
| BCL2L11 pan-cancer check | declared in `5cd7111` together with the Block G result and its models | `5cd7111` | **same commit** | a check, not a hypothesis |
| Block D | model in plan `b542e70`; readings in script 15 (`cbeed0d`) | `cbeed0d` | model yes | "NO DIRECTION IS PRE-SPECIFIED" |
| F3-pre | G3 note §4b, same commit as script 15 (`cbeed0d`) | `cbeed0d` | same commit | - |
| F3 proper | plan | never run | - | - |
| **H4** | `b542e70` (categorical), `67f144a` (continuous primary, 08-28), `f1f0c0a` (STATE, BUFFER_c, 08-30 23:23:06) | `259efa9` (08-31 00:20:21) | **yes**, for the direction and primary model | amendment `f299c59` came after the pCR columns were read but before association; A4/A5 arrived in the script commit; `cc34f74` changed the script after its first run |
| BIM replication (SCAN-B) | block_c §9 (`f2def97`, 08-29, written after the TCGA BIM number); SCAN-B declaration `a0c5519` (08-31 11:47:46) | `37816ea` | **yes** | symbol addenda came after download, before modelling |

**Where no declaration exists:**
- The BCL2L11 pan-cancer check: its declaration and models share a commit.
- Block D's direction: none was ever specified.
- G2's numeric threshold: set after crude ORs were seen.
- The F3 gate decision: never recorded.

---

## 7. The open decisions the build specification left

Section 14 of the plan lists them. Its frontmatter at HEAD reads "resolved: D1,
D2, D3, D4, D5, D7" and "open: D6, G1's correlation criterion, script 07's
specificity panel source" (plan:L11-12).

| decision | how it was settled | source |
|---|---|---|
| **Primary MYC estimator (D2)** | **settled: M-a, the Felsher signature stripped of MitoCarta (61 genes)**, "**Confirmed, not flipped.**" (`2026-08-28_G1_result_and_decisions.md:L85`, `eb06ee6`, 08-28), then "on evidence" by the estimator validation (M_a vs M_b rho 0.613). **plan:L1089-1091 still words D2 as an open question**; only the frontmatter marks it resolved | G1 note; myc_estimator_validation |
| **Lee et al. 2017 scoping (D4)** | **settled: "Option 1, confirmed 2026-08-28"** (`2026-08-28_D4_scoping_and_G2_design.md:L138`). The verdict is "SUPPORTS, with H1 clause 2 partially pre-empted for MCL1 only". **It blocked G2**: "**D4 - Lee et al. 2017 scoping. BLOCKS G2.**" (G1 note:L216), and plan:L1098 "**Do this before G2.**" D4 was resolved in `77968ac` (11:11:54), before the G2 result `402b164` (11:35:49). The confirmer is not named. **plan:L1095-1098 is not updated** | D4 note |
| **Primary neoadjuvant cohort (D5)** | **settled: I-SPY2-990, GSE194040, n = 988**, with BrighTNess as replication and GSE25066 third ("**D5 - RESOLVED 2026-08-28.**", plan:L1099; `67f144a`). The original plan recommended GSE25066 as primary (`plan@b542e70:L803-806`) | D5 cohort selection note |
| **Fourth, unresolved: D6** | **NOT SETTLED.** "**D6 - NEW. Does the mouse arm have a metastasis phenotype to anchor F4 to?** ... If not, drop F4." (plan:L1123-1125). Open in the frontmatter; the final handoff lists D6 among outstanding, non-gating items. No D6 note exists | plan; handoff |
| also open at HEAD | G1's correlation criterion (frontmatter only); "script 07's specificity panel source" (frontmatter only, contradicted by plan:L168, "RESOLVED"); the script 09 spec sign-off; the F3 gate decision | plan:L12 |

D1, D3, D7, D8 and D9 are resolved. D8 and D9 are recorded in the plan body
but are missing from the frontmatter's list.

---

## 8. What the manuscript's human section can rely on, and what it does not use

**The manuscript itself is not in this repository**, so this section reports
what the repo says it supplies. It does not check that against the manuscript.

**What the repo says it feeds.** "Feeds the manuscript "Mitochondria integrate
oncogenic and metabolic transcriptional programs to shape breast cancer
progression" (Nature Metabolism Letter), as 2-3 main panels plus Extended Data"
(CLAUDE.md:L13-16).

**What it now offers** (`2026-08-31_human_arm_writeup_draft.md`, `7707aeb`):
- One Extended Data figure and a single main-text paragraph: "one Extended Data
  figure" (L19). The handoff calls this "**Not the three main panels plan
  section 4 allocated.**" and leaves it as "the author's call" (L141-147).
- The content is a well-formed negative:
  - H1-H4 not supported;
  - G2's co-amplification (adj OR 1.69-1.77), framed as replication of Lee et
    al. for MCL1 and as novel for BCL2L1;
  - three of four falsification criteria met;
  - BIM positive in TCGA and inverted in SCAN-B.
- The draft's own summary: a pre-registered test "**did not find the axis in
  established human breast cancer**, and said so" (L181-186).

**What the repo contains but says is not to be used:**
- **Five results "deliberately NOT pursued"** (`2026-08-31_handoff.md:L112-116`;
  CLAUDE.md:L60-63):
  - `Glycine metabolism` (Block C);
  - `Folate and 1-C` (Block B);
  - `D MB2 DESCRIPTIVE log2_BCL2L1` (Block D);
  - MCL1 in BrighTNess (F1);
  - the negative BIM association in SCAN-B, with the explicit instruction
    "**Never test whether the sign returns in the basal or ER-negative
    subset.**" (handoff:L118).
- The curated metabolic gene list: "reported not used".
- TCGA survival: "Do not run the survival analysis in TCGA" (plan:L248).
- The mammary gene-set library's human GMT tree: must not be loaded.
- The mouse result tables, which are never analysis inputs.
- Planned but not run: F2, F3 proper, F4, Block D2, Block E (RPPA). The plan's
  script list (plan:L1012-1013) names scripts 16 and 17 differently from the
  ones that exist.

---

## 9. What contradicts, or is more recent than, an outside account of the arm

**No project brief was supplied with this task, and none exists in the
repository.** `git grep -i brief HEAD` matches only three uses of the word
"brief" in two vendor README files under `docs/`. So no direct comparison is
possible here. The points below are what any brief must agree with to be
correct about this repo as frozen.

**Latest state.** The latest dated content is from 2026-08-31 13:24:15
(`d3ac60e`). Any account written from the `6dc6a47` handoff (00:48:49 the same
day) is out of date on three points:
- what survives: `6dc6a47` says "**lower BCL-XL and higher BIM**";
- the SCAN-B replication: blocked there, run and failed by `d3ac60e`;
- the write-up framing: undecided there, decided by `d3ac60e`.

**Contradictions inside the frozen repo:**
1. **H4 instrument count.**
   - The write-up draft calls the test "two-instrument" (L182).
   - The amendment says H4 is "GSVA alone" (amendment:L88-90).
   - It also says a positive or null H4 "must be reported as single-instrument"
     (L102-105).
2. **H4 sample size.**
   - The draft says 1,978 neoadjuvant patients (L50).
   - The F1 result fitted 988 + 482 + 470 = 1,940 (F1 result:L27).
3. **Stale plan text.**
   - "only H4 is untested" (plan:L951).
   - The F1 model still lists stage and grade (plan:L825), which amendment A4
     replaced with subtype and treatment.
   - D2 and D4 are still worded as open questions (L1089-1098).
   - G1's correlation criterion is still listed as open (L12).
   - The plan's own status is still "draft-for-review".
4. **Moving lock on the H4 declaration.**
   - The STATE note's lock was "after any of those outcome columns has been
     read" (L17-18).
   - The amendment moved the boundary to "score-versus-pCR association".
5. **The script 09 spec** is still "DRAFT - awaiting sign-off" (L3), although
   script 09 implements it and its result is recorded.
6. **The F3 stop gate was bypassed by a decision that is not recorded** (section
   3). So H4 ran without the guard the plan placed before it.
7. **Small numerical inconsistencies:**
   - SCAN-B mitoPPS p is 0.037 at result:L25 and 0.036 at L82.
   - Block C coefficient counts are 116 at L20 and 114 at L276.
   - Block G pan-cancer n is 1,130 at L37 and 1,127 at L170.
8. **"Frozen" is not stated anywhere in the repo.** It ends expecting a
   write-up and a script 18.

---

## 10. What I could not verify

- When anything was actually run. Git orders commits, not executions;
  `results/` is not in git and no R was run, so **no number was reproduced**.
  Every number is as the notes print it.
- The external record that names `b8ef43f` (no file in the repo cites it).
- Mouse-repo facts the repo relies on: the BIM pre-specification in
  `myc_mouse/scripts/44`, and the mouse percentiles. Other repositories were out
  of scope.
- Who confirmed D4, D5, the script 09 spec, or the F3-gate bypass. The notes
  name no decider for most of them.
- The G1 correlation-criterion numbers: whether they predate or postdate the
  averaged-aliquot fix (msg `ccc2ef1`).
- Whether the declared-but-unreported H4 sensitivities were run.
- Scripts 13, 16 and 17 were not read in full.
