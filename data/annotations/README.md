# data/annotations

Small tables written by hand from a cited source. They are not snapshots of an
external data release, so they carry a citation instead of a pinned commit SHA.
Edit them by hand only, and record who verified each row and when.

## bh3_binding_code.csv

**What it is.** It records which BH3-only proteins bind BCL-xL (`BCL2L1`) and
MCL1. It feeds the binding-code columns of Fig. 4g:
- `docs/2026-11-03_signature_and_fig4g_brief.md`, section 6.6;
- the E49 declaration, section 10.

**Columns:**
- `gene`, the HGNC symbol;
- `binds_BCL2L1` and `binds_MCL1`: `TRUE`, `FALSE`, or blank;
- `source`, the citation and the sentence of the abstract that supports the row;
- `verified_by`, who checked it and when.

**Source.** Chen L et al., *Mol Cell* 2005;17:393-403. PMID 15694340,
doi:10.1016/j.molcel.2004.12.030. The paper measured the affinities of eight BH3
peptides for five prosurvival proteins.

**The verified rows are `BBC3`, `BCL2L11`, `BAD` and `PMAIP1`.** Each one rests
on the abstract, which its own `source` cell quotes. They were verified as
follows:
- the author checked them on 2026-11-03;
- Claude Code re-checked them against the PubMed abstract the same day.

**The blank rows are `BID`, `BIK` and `BMF`.**
- The abstract does not state their classes.
- The author supplies them from the paper's full binding data.
- Until then their cells stay blank, and nothing is drawn for them.

**A binding class is never inferred.** It is never taken from a review, a
homologue or a pattern in the data. A blank cell is drawn blank.
