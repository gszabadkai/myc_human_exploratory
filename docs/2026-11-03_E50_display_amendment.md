---
date: 2026-11-03
status: >
  The author's further amendment to the STOP 3 rulings on E49, recorded
  verbatim. It is the specification of E50, a display-only script, and was
  committed alone, before E50 was written.
relates-to:
  - docs/2026-11-03_E49_signature_declaration.md (7bd7865)
  - docs/2026-11-03_e49_result.md (4fd661b)
  - docs/2026-11-03_e49_handoff.md (4070759)
---

# E50's specification: the author's further amendment to the STOP 3 rulings

## 1. What this is

**The text in section 2 is the author's, pasted into the session on 2026-11-03 and recorded unchanged.** It is E50's specification. E50 is display-only: it reads E49's verified saved object and fits nothing.

**It amends rulings that are not in the repo.** The author's STOP 3 rulings on E49 (items 1 to 6, and amendments A to E) did not reach the session that recorded this. Section 2 replaces three of them:
- the version choice (item 1);
- the grey cells (item 2 and amendment E);
- the E50 specification (item 6).

**Everything else in those rulings stands, including the text.** Applying the rest (the text in the review pack and in the result note's Part 4) waits until the rulings themselves are in hand.

**How E50 reads the specification, and every choice it had to make, is recorded in the script's header and in the pre-run digest note.** Nothing in section 2 is edited here.

## 2. The amendment, verbatim

> Further amendment to my STOP 3 rulings (the author, 2026-11-03). It replaces the version choice (item 1), the grey cells (item 2 and amendment E) and the E50 specification (item 6). Everything else stands, including the text.
>
> 1. No grey anywhere. Every cell is drawn in its actual colour: 4g's signature row, ED-e, ED-f, ED-g and the signature rows of ED-h. In each legend, replace the grey sentence with: "The signature was selected for its association with OXPHOS at fixed proliferation (Methods)."
>
> 2. Palette: Fig. 1g's. If Fig. 1g's source script is in the mouse repo, reuse its palette and drawing code. Otherwise use the values I measured from the Fig. 1 PDF:
>    - brown #462915 at the negative limit, #F9F9F9 at zero, green #318866 at the positive limit, linear in RGB on each side;
>    - limits symmetric, rounded outward to 0.1 (±0.6 on present values).
>
>    Use this palette for every human coefficient or correlation heatmap: 4g, ED-e, ED-f, ED-g and ED-h. ED-a keeps its band shading.
>
> 3. 4g, ED-e and ED-f become heatmaps in Fig. 1g's style, with no annotation columns. I will add the binding code, the signature signs, the § marks and any protein names in Illustrator.
>    - Rows: the 15 genes in Fig. 1g's groups and order, gene symbols in Helvetica oblique on the left. Group names (Effectors, Activators, Sensitisers, Anti-apoptotic) on the right, as in Fig. 1g. The signature as a final row, labelled "Signature", after a larger gap.
>    - Columns: OXPHOS, MYC activity and proliferation, each split into TCGA and SCAN-B; a gap of about 2.3 pt between the three blocks. Two header tiers: the score, then the cohort.
>    - Cells: 26.71 x 8.46 pt, no borders, no gaps within a group, about 2.3 pt between groups.
>    - Values printed in every cell as in Fig. 1g: signed, one decimal, "0" when the value rounds to zero.
>    - Type: Helvetica 5 pt for values and labels, 6 pt for the score headers. These panels follow Fig. 1g's 5 pt rather than the 5.5 pt floor used so far.
>    - A vertical colour bar on the right, titled "SD per SD", with signed ticks.
>    - No dots. The intervals stay in ED-i.
>
> 4. Output: a vector PDF with live Helvetica text (not outlined) and an SVG, both editable in Illustrator, plus a PNG preview, for 4g, ED-e and ED-f. ED-g and ED-h are re-coloured only.
>
> 5. E50, display-only, reads E49's saved object and implements items 1 to 4. It no longer draws the binding code. The table stays the record behind my Illustrator annotations, and I will fill in BID, BIK and BMF there.
>    - Assertions: every printed value equals its stored value, rounded; the palette's three anchors equal the colours in item 2; the cell size, gaps and fonts match item 3.
>    - The usual sequence: dry runs, digests committed before my run, then verification. Then regenerate the review pack with the amended text, and update CLAUDE.md and the handoff.

## 3. What this supersedes

- **E49 declaration section 9.1's colour rule**, "not Fig. 1g's brown to green, because the estimand differs". The author now uses Fig. 1g's palette for every human coefficient or correlation heatmap.
- **E49 declaration section 8.2's grey.** The signature's OXPHOS coefficient, and the difference that contains it, are now drawn in their actual colour. **Nothing in the text reads them**, and each legend states how the signature was selected.
- **The two versions of 4g, ED-e, ED-f and ED-h.** There are no dots. The intervals stay in ED-i.
- **The binding-code, signature-sign and § columns in the drawn panels.** The author adds them in Illustrator. `data/annotations/bh3_binding_code.csv` stays the record behind them.
