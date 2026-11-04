---
date: 2026-11-04
status: >
  The author's ruling on the binding classes of BID, BIK and BMF (the STOP 3
  rulings, section 5), and on the citations that follow from it. Recorded
  verbatim and committed alone, before it is applied.
relates-to:
  - docs/2026-11-03_stop3_rulings_consolidated.md (cb9b7ce; section 5, and legend edit 6 in section 4)
  - docs/2026-11-03_stop3_rulings_addendum.md (3f9ffbe)
  - data/annotations/bh3_binding_code.csv and its README
---

# The binding classes of BID, BIK and BMF (the author, 2026-11-04)

## 1. What this is

**Sections 2 and 3 are the author's, recorded unchanged.**
- **Section 2** is the author's message on resuming the session on 2026-11-04.
- **Section 3** gives the author's answers to two questions Claude Code asked about it. Each question and the option the author chose are quoted as they were shown.

**What it changes:**
- **The source.** The rulings' section 5 named Chen et al. 2005's full binding data. The three rows come instead from Ku et al. 2010 (doi:10.1038/cr.2010.149) and Wang et al. 2023 (doi:10.1016/j.csbj.2023.07.017).
- **BIK and BMF against MCL1 are coded "does not bind".** Their binding is weak, not absent.
- **The 4g legend and Methods note N2 cite all three papers.** This is a text change to the rulings' legend and to N2.

**Nothing in sections 2 and 3 is edited here.** It is applied in a separate commit: the binding-code table, `pack_text.py`, the E50 review pack and E49's result note, Part 4.

## 2. The author's message, verbatim

> Start with the handoff: docs/2026-11-03_e50_handoff.md  1. I will do the illustratir work. 2. Summary: Together, published binding metrics establish that BID binds both proteins potently (favoring BCL-xL functionally), whereas BIK and BMF strongly favor BCL-xL over MCL-1. https://doi.org/10.1038/cr.2010.149, https://doi.org/10.1016/j.csbj.2023.07.017 4. I will do the rest later.

## 3. The author's answers, verbatim

**Question 1:**

> The binding-code table has two classes per protein: binds (filled) or does not (open). You say BIK and BMF strongly favour BCL-xL over MCL-1, and Wang 2023 measured weak Bmf binding to Mcl-1 and solved that complex. How should BIK and BMF be coded against MCL1? (BID will be 'binds both'; BIK and BMF 'binds BCL-xL'.)

**The author's choice:**

> Does not bind (open): binds_MCL1 = FALSE for BIK and BMF, as for BAD. The table's note records that the binding is weak, not absent.

**Question 2:**

> The 4g legend says the binding code is 'from Chen et al. (2005)', and Methods note N2 cites only [Chen]. With BID, BIK and BMF now from Ku 2010 and Wang 2023, what should those two places cite?

**The author's choice:**

> Add both papers: The legend reads 'from Chen et al. (2005), Ku et al. (2010) and Wang et al. (2023)', and N2 cites all three. This is a text change to your ruled legend and N2, recorded as your ruling.
