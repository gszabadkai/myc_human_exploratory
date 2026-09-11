# The level-1 category comparison — what it is for and why it is designed this way

**2026-09-11.** Background for the Claude Code sessions that build it. Read this
before the prompt.

---

## 1. Where this came from

Phase 2 of the human arm compared mitoPPS **pathway rank orderings** across
species. Its conclusion: human breast tumours match the mouse **MYC** move and
not the mouse **developmental** move, with MYC-low tumours ordering the
mitochondrial compartment like a mature wild-type gland and MYC-high tumours
ordering it like a young MYC-driven one.

While reading `E25_group_orderings.pdf`, the author noticed something the
whole-panel correlation had obscured. Two readouts in the same human samples
behave differently:

| readout | 12W WT | 6W WT | 6W Myc+ | 12W Myc+ | TCGA MYC-high | SCAN-B MYC-high |
|---|---|---|---|---|---|---|
| **Mitochondrial ribosome** | 19 | 35 | **80** | 62 | **81** | **89** |
| **OXPHOS subunits** | 7 | 61 | **89** | 33 | **68** | **74** |
| mtDNA-encoded | 100 | 7 | 0.4 | 99 | 13 | 12 |

*(percentile within the 142-pathway shared panel, group-mean mitoPPS ordering)*

**Human MYC-high tumours reach or exceed the mouse mitoribosome position but
fall ~20 percentile points short on the respiratory chain.** The MYC
mitochondrial-translation programme goes all the way; the respiratory
prioritisation does not.

That is the **mass-versus-function dissociation** the paper's title is about,
appearing in human tumours, in both cohorts, as a gap between two readouts that
MYC drives together in the mouse gland.

**Caveat carried from the start:** 6W WT's interval is [20, 93] at n = 6, so
"human MYC-high sits at 6W WT level" is not supportable. What is supportable is
that it sits clearly above 12W WT and MYC-low and below 6W Myc+, with the
shortfall clearer in TCGA ([57, 77] against [71, 91]) than in SCAN-B.

---

## 2. What this analysis is meant to demonstrate

The mouse timeline's developmental claim is currently made on **selected gene
sets scored against expression-matched nulls**. The author's proposal is to make
it instead on the **native MitoCarta level-1 partition**, using three categories
with three different expected behaviours across 6W → 12W:

| category | expected across the WT window |
|---|---|
| **OXPHOS** | **down** — deprioritised |
| **Mitochondrial central dogma** | **flat** — the internal null |
| **Metabolism** | **up** |

**The null category is the point.** mitoPPS is a relative score, so "OXPHOS
down" and "Metabolism up" are not independent observations. Without a category
that does **not** move in the same samples on the same ruler, the pattern could
be an artefact of relativity rather than a reallocation. Central dogma is what
makes it a reallocation claim.

---

## 3. Three design decisions, and why

### 3.1 Report all seven level-1 categories, not three

mitoPPS is relative: if OXPHOS falls, the average of everything else must rise.
Showing that Metabolism rises is therefore only informative if the rise is
**not spread evenly** — and with three categories that cannot be demonstrated,
because the other four are unobserved.

The MitoCarta level-1 partition is exactly seven, and E25 confirmed both species
carry all seven (7 level-1, 38 level-2, 96 level-3, agreeing across species).
Reporting the complete partition converts "three chosen sets" into "the complete
top-level partition, and here is where the priority went". The extra panels cost
nothing and the figure becomes much harder to argue with.

The three named categories remain the **pre-declared readouts**; the other four
are reported as context.

### 3.2 Within a species, use mitoPPS values — not ranks

Ranks were needed in Phase 2 only to cross the species boundary. Inside one
cohort there is no such constraint, and the rank transform introduces two
problems the raw scores do not have:

- **within-sample ranks are strictly zero-sum**, so "one down, one up" is partly
  forced by construction — which is precisely the artefact the null category
  exists to exclude;
- `E24` demonstrated that the within-sample rank transform **carries a size
  bias** (rho against log set size +0.54 and +0.57) that the underlying scores do
  not (mitoPPS is centred at exactly 1 per pathway in both species, and both
  implementations take the mean per set).

So: **group-mean mitoPPS per category with bootstrap intervals.** This is E25's
estimand family, not E24's.

### 3.3 The design narrows the batch objection, and that should be stated

Batch equals timepoint in the mouse timeline and is carried as a standing
caveat. The author has ruled that timeline comparisons stay in, that DESeq2
normalisation is what handles the batch, and that no number is withdrawn on its
account. **Do not reopen that.**

But this design does something the caveat alone does not: a batch effect would
have to act **differentially across three mitochondrial categories in the same
libraries** — lowering OXPHOS, raising Metabolism, leaving Central Dogma
untouched. That is a far stronger requirement than a global batch effect. It
does not eliminate the confound; it narrows what the confound would have to be,
and the figure makes that visible rather than requiring a paragraph.

---

## 4. The statistical point about "flat"

At n = 6, an interval covering zero is **not** evidence that Central Dogma does
not move. It may simply be low power, and if that is all it says, the null panel
is doing no work.

What is wanted is an **equivalence-style reading**: Central Dogma's interval
excluding the *magnitude* of the OXPHOS change. State it that way explicitly.

If the interval is too wide for that, the honest form is *"Central Dogma's change
is not resolvable at this n"* — weaker, but still a useful contrast against a
clearly resolved OXPHOS fall, and better than implying a demonstrated null.

---

## 5. Why the human companion belongs with it

The **mitochondrial ribosome sits inside Mitochondrial central dogma**. So the
human dissociation in §1 is already a two-category version of this comparison,
run on the wrong ontology level.

Running the same seven categories in both species gives one coherent statement
across two axes:

> **Development** lowers OXPHOS priority, leaves central dogma alone, and raises
> metabolism. **MYC** raises central dogma and OXPHOS together in the mouse
> gland — but in human tumours central dogma goes all the way and OXPHOS stops
> short.

That makes the human panel a **companion to the mouse panel on the same
ontology**, rather than a separate cross-species argument carrying its own
caveats. It also connects directly to the paper's existing finding that
biogenesis and OXPHOS are largely different gene programmes (~188 against ~249
genes, overlap 17) and that OXPHOS retains coupling to outcomes where biogenesis
is a MYC-dose bystander.

---

## 6. Checks to run, not obstacles to clear

- **Category overlap.** The author's expectation is no overlap between
  Metabolism and OXPHOS. **Check and report it; do not stop on it.** If there is
  overlap, note its size and proceed — the reallocation reading weakens
  proportionally and the note says so.
- **Category names and sizes** read from the source workbook, not from memory.
  Metabolism is expected to be much the largest, and size should be reported
  beside every category.
- **Containment structure** is already characterised: E25 built a 114-member
  antichain from each species' own hierarchy and the two were identical.

---

## 7. What this is, and is not

**It is** a cleaner presentation of a claim the paper already makes, on the
native ontology, with an internal null.

**It is not** a new result, a hypothesis test, or a replacement for the
expression-matched-null analysis that established the developmental change in
the first place. Both should be reportable; if they disagree, that disagreement
is the finding and must not be resolved by choosing the friendlier one.

**N3 throughout.** These are transcript associations. "Primed" appears nowhere.
