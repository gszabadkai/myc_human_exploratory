---
date: 2026-11-03
status: >
  The author's addendum to the STOP 3 rulings, answering the four observations
  left after the rulings were applied (da548bd). Recorded verbatim and
  committed alone, before it is applied.
relates-to:
  - docs/2026-11-03_stop3_rulings_consolidated.md (cb9b7ce, the rulings it amends)
  - docs/2026-11-03_e49_result.md (Part 4, as of da548bd, where the observations are)
  - docs/2026-10-22_E47_mutual_adjustment_declaration.md (4a.7, item 5, the Methods definition it adopts)
  - docs/2026-11-03_signature_and_fig4g_brief.md (eb45a9b; rule 5.4, which it rewords)
---

# Addendum to the STOP 3 rulings (the author, 2026-11-03)

## 1. What this is

**Section 2 is the author's text, pasted into the session on 2026-11-03 and recorded unchanged.** It answers the four observations on the front page of the E50 review pack and in E49's result note, 4.1 (`da548bd`).

**What it amends:**
- **the rulings' section 3:**
  - R2's sample sizes;
  - the Methods notes N4 to N6, which are dropped, with E47 declaration 4a.7, item 5, adopted as the Methods definition, with two changes;
- **the brief's rule 5.4**, reworded. The brief itself is not edited, as with earlier amendments;
- **nothing drawn.** ED-g keeps its ±0.9 scale, and the author deletes its "Proliferation minus OXPHOS" column in Illustrator.

**Nothing in section 2 is edited here.** It is applied in a separate commit: `pack_text.py`, the E50 review pack and E49's result note, Part 4.

## 2. The addendum, verbatim

> Addendum to my STOP 3 rulings (the author, 2026-11-03), answering your four observations. Record it verbatim in its own dated doc, committed alone. Then extend check_against_rulings() to read it, apply it in pack_text.py, regenerate the E50 pack, and update E49's Part 4 (4.1 and 4.3).
>
> 1. R2: replace the sample sizes with the cohorts the search used:
>    "We therefore searched the twelve BH3-only and anti-apoptotic genes of the BCL-2 family in TCGA-BRCA (n = 1,095) and SCAN-B (n = 3,207) for those whose expression tracked OXPHOS at fixed proliferation."
>    The model sets (938 and 3,143) stay in the 4g legend, and the abstract's 4,081 stands.
>
> 2. Brief rule 5.4 now reads: "BAD is reported with its explanation in the sentence that sets its rise against BCL-xL's, or the next." R7 meets it. R4 stands unchanged.
>
> 3. Methods: N4, N5 and N6 are dropped. The definition is E47 declaration 4a.7, item 5, verbatim, with two changes:
>    a. Its sentence "Because genes and signs were chosen by their association with OXPHOS in the same cohorts, the signature's association with OXPHOS is expected by construction; the informative quantity is its association with MYC activity at fixed OXPHOS." becomes:
>       "Because genes and signs were chosen by their association with OXPHOS at fixed proliferation in the same cohorts, the signature's association with OXPHOS, and in part with proliferation, is expected by construction. The informative quantity is its association with MYC activity at fixed OXPHOS and proliferation."
>    b. Its variants sentence ends with "(Extended Data Fig. Xg)".
>    N1, N2 and N3 stay. Measure the paragraph and report any sentence outside STYLE.md's range for me to split; do not drop the construction sentence. N6 as ruled misdescribed S5, whose rule uses MYC activity; E47's paragraph states it correctly.
>
> 4. ED-g keeps its ±0.9 scale. I delete the "Proliferation minus OXPHOS" column in Illustrator; no legend edit and no new script.
