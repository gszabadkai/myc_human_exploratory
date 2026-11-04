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
- `source`, the citation and the sentence that supports the row: the abstract's,
  or the author's summary where the abstract gives none;
- `verified_by`, who checked it and when.

**Sources.**
- **Chen L et al.**, *Mol Cell* 2005;17:393-403. PMID 15694340,
  doi:10.1016/j.molcel.2004.12.030. The paper measured the affinities of eight
  BH3 peptides for five prosurvival proteins. It is the source of `BBC3`,
  `BCL2L11`, `BAD` and `PMAIP1`.
- **Ku B et al.**, *Cell Res* 2011;21:627-641. PMID 21060336,
  doi:10.1038/cr.2010.149.
- **Wang H et al.**, *Comput Struct Biotechnol J* 2023;21:3760-3767.
  PMID 37560128, doi:10.1016/j.csbj.2023.07.017.

  Ku et al. and Wang et al. are the sources of `BID`, `BIK` and `BMF`.

**All seven rows are verified.**
- **`BBC3`, `BCL2L11`, `BAD` and `PMAIP1`** rest on Chen et al.'s abstract,
  which each row's `source` cell quotes.
  - The author checked them on 2026-11-03.
  - Claude Code re-checked them against the PubMed abstract the same day.
- **`BID`, `BIK` and `BMF`** rest on the author's summary of Ku et al. and Wang
  et al. (`docs/2026-11-04_binding_classes_ruling.md`, 2026-11-04).
  - **`BID` binds both proteins.**
  - **`BIK` and `BMF` bind BCL-xL.** Against MCL1 they are coded `FALSE` by the
    author's decision, because they strongly favour BCL-xL. **Their MCL1 binding
    is weak, not absent.** Wang et al. measured a much lower affinity of the Bmf
    BH3 peptide for Mcl-1 and solved that complex.
  - **Claude Code checked both abstracts in PubMed on 2026-11-04.** Wang et al.'s
    is quoted in the `BMF` row. Ku et al.'s gives no values for these rows, and
    its full text was not retrievable from PMC.
  - **The rulings' section 5 had named Chen et al.'s full data.** The author's
    ruling of 2026-11-04 replaces it for these rows. It also adds Ku et al. and
    Wang et al. to the 4g legend and to Methods note N2.

**Where the code is drawn.** On 4g, only through the author's Illustrator
annotations (the STOP 3 rulings, 1.3 and 5). E50 does not read this table.
- **E49 read it when it ran** (`a63ef6c`), with the `BID`, `BIK` and `BMF` cells
  blank. Its saved `binding_code` and its digests record that state.
- **E49 is not re-run for this.** Re-sourcing it now would change that element
  and its displays, so they would no longer match E49's verified digests.

**A binding class is never inferred.** It is never taken from a review, a
homologue or a pattern in the data. A blank cell is drawn blank. Every class
here is the author's, from the cited papers.
