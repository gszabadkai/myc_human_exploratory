# Writes docs/2026-11-03_e50_review_pack.qmd. It reads:
#   - the text, from E49's single source, ../2026-11-03_e49_review_pack_src/pack_text.py
#     (the result note's Part 4, the proposed legend additions, the metrics);
#   - this folder's legends.json, E50's $legend_draft verbatim (4g, ED-e to ED-h);
#   - E49's and E48's legends.json, verbatim (ED-i to ED-k; ED-a to ED-c).
import sys, json, importlib.util
S, OUT = sys.argv[1], sys.argv[2]
sys.dont_write_bytecode = True
spec = importlib.util.spec_from_file_location("pt", S + "/../2026-11-03_e49_review_pack_src/pack_text.py")
pt = importlib.util.module_from_spec(spec); spec.loader.exec_module(pt)
LEG50 = json.load(open(S + "/legends.json", encoding="utf-8"))
LEG49 = json.load(open(S + "/../2026-11-03_e49_review_pack_src/legends.json", encoding="utf-8"))
LEG48 = json.load(open(S + "/../2026-11-03_e48_review_pack_src/legends.json", encoding="utf-8"))

W_HEAT = "87.6mm"   # E50's heatmaps at final size (87.6 x 58.2 mm)
W_PAGE = "178mm"

def nwords(s):
    return len([t for t in s.split() if any(c.isalnum() for c in t)])

def metrics_table():
    rows = ["| paragraph | sentences | words | mean | median | range | under 20 | under 12 | em-dashes |",
            "|:------------------------------------------|--------:|-----:|-----:|------:|--------:|---------:|-------:|-------:|"]
    for short, m in pt.metrics_rows():
        rows.append(f"| {short} | {m['sentences']} | {m['words']} | {m['mean']:.1f} | {m['median']:g} | "
                    f"{m['min']}–{m['max']} | {m['under20']} of {m['sentences']} | {m['under12']} | {m['emdash']} |")
    return "\n".join(rows)

def results_summary():
    res = [t for (title, short, ss) in pt.PARAS
           if short in ("Signature paragraph", "Balance paragraph", "MCL1 paragraph") for (_, t, _) in ss]
    m = pt.metrics(res)
    pct = 100 * m["under20"] / m["sentences"]
    return (f"**The three Results paragraphs together:** {m['sentences']} sentences, {m['words']} words, "
            f"mean {m['mean']:.1f}, median {m['median']:g}, range {m['min']}–{m['max']}, "
            f"{m['under20']} of {m['sentences']} ({pct:.0f}%) under 20 words. "
            f"**That share is below STYLE.md's 35% floor**; no sentence was shortened for the metric." if pct < 35 else "")

# (title, [(figure, width)], legend, cited as)
ED = [
  ("ED-a (E48): the scores' correlations", [("figures/2026-11-03_E48_score_correlations.png", W_PAGE)],
   LEG48["ed_a"], "Extended Data Fig. Xa, in sentences S2b and S6"),
  ("ED-b (E48): *BCL2L1* against *BBC3*", [("figures/2026-11-03_E48_bcl2l1_against_bbc3.png", W_PAGE)],
   LEG48["ed_b"], "Extended Data Fig. Xb, in sentence B1"),
  ("ED-c (E48): MYC activity on four scores, three genes", [("figures/2026-11-03_E48_myc_score_ladder.png", W_PAGE)],
   LEG48["ed_c"], "not cited in the text. ED-f sits beside it"),
  ("ED-e (E50): Fig. 4g on OXPHOS priority", [("figures/2026-11-03_E50_bcl2_family_heatmap_priority.png", W_HEAT)],
   LEG50["ed_e"], "not cited in the text; the 4g legend points to it"),
  ("ED-f (E50): Fig. 4g from the reference model", [("figures/2026-11-03_E50_bcl2_family_heatmap_reference_score.png", W_HEAT)],
   LEG50["ed_f"], "not cited in the text; the 4g legend points to it"),
  ("ED-g (E50): the signature and its two variants, re-coloured", [("figures/2026-11-03_E50_signature_strips.png", W_PAGE)],
   LEG50["ed_g"], "not cited in the text"),
  ("ED-h (E50): raw correlations beside the adjusted coefficients, re-coloured", [("figures/2026-11-03_E50_raw_beside_adjusted.png", W_PAGE)],
   LEG50["ed_h"], "Extended Data Fig. Xh, in sentence S6"),
  ("ED-i (E49): all fifteen genes, three coefficients and the difference", [("figures/2026-11-03_E49_fifteen_genes_intervals.png", W_PAGE)],
   LEG49["ed_i"], "not cited in the text"),
  ("ED-j (E49): MYC activity on four scores, the signatures and the effectors", [("figures/2026-11-03_E49_myc_ladder_signatures_effectors.png", W_PAGE)],
   LEG49["ed_j"], "Extended Data Fig. Xj, in sentence S7"),
  ("ED-k (E49): PUMA relative to BCL-xL against OXPHOS priority", [("figures/2026-11-03_E49_puma_bclxl_against_priority.png", W_PAGE)],
   LEG49["ed_k"], "Extended Data Fig. Xk, in sentence B2"),
]

q = []
q.append("""---
title: "E50 review pack: Fig. 4g in Fig. 1g's style"
subtitle: "For the author. Exploratory and post hoc. The text still awaits the STOP 3 rulings."
date: 2026-11-03
format:
  typst:
    papersize: a4
    margin:
      x: 16mm
      y: 15mm
    fontsize: 9pt
---

**What this is.** Fig. 4g as E50 draws it, at final size (87.6 mm wide), with the text, the legend draft and the STYLE.md metrics, followed by the Extended Data panels.

- **E50 follows the author's amendment** (`docs/2026-11-03_E50_display_amendment.md`):
  - no grey anywhere;
  - Fig. 1g's palette;
  - 4g, ED-e and ED-f in Fig. 1g's style, with no annotation columns and no dots;
  - ED-g and ED-h re-coloured.
- **The binding code, the signature signs, the § marks and any protein names** are the author's, added in Illustrator.
- **The panels shown are the PNG previews.** The PDF and SVG beside them in `docs/figures/` are the files for Illustrator, with live Helvetica text.

**THE TEXT IS NOT YET AMENDED.**

- **What it is:** E49's result note, Part 4, unchanged (`4fd661b`).
- **Why:** the amendment says the STOP 3 rulings on the text still stand. Those rulings (items 1 to 6, amendments A to E) are not in the repo and did not reach this session.
- **When they arrive:** the sentences are edited in their one source, E49's `pack_text.py`, and this pack is regenerated.
- **What has already changed:** only the legend drafts, as E50 changed them.

**The record.**

- **E50 is verified:** 116 of 116 assertions, and 11 of 12 saved elements by digest. The twelfth differs only in one label's wording.
  - **All comparable figure files match the checksums committed before the run** (`eae8ef5`).
  - **The details:** `docs/2026-11-03_e50_result.md`.
- **Nothing new is computed.** Every number drawn was stored by E49 (`4fd661b`).

**The palette is Fig. 1g's source, not the measured values.** Item 2 puts the source first.

| anchor | Fig. 1g's source (used) | measured from the Fig. 1 PDF (not used) |
|---|---|---|
| negative limit | `#4A3525` | `#462915` |
| zero | `#FAFAFA` | `#F9F9F9` |
| positive limit | `#2A8A6D` | `#318866` |

Fig. 1g's scale interpolates in CIELAB, and E50 follows it.

**Choices the amendment left open** (E50's pre-run note, section 3):

- **One scale for 4g, ED-e and ED-f, ±0.6.** ED-e alone would be ±0.5.
- **ED-g's scale is ±0.9,** set by its "proliferation minus OXPHOS" cells. ED-h stays at ±0.7 for the raw correlations and ±0.6 for the adjusted coefficients.
- **The gap before the Signature row** is twice the group gap.
- **Ink is black, or white where white reads better:** 9 of 288 cells.
- **The group names and the bar's title are 5 pt** (item 3's two sizes). Fig. 1g sets them at 6 pt.

**One qualification found after E49 ran** (E49 result, Part 2 H2; AD HOC arithmetic, not declared).

- **The issue.** Part of the signature's proliferation coefficient can be inherited from its selection rule. The rule chose genes by their correlation with OXPHOS at fixed proliferation.
- **How much.** At the extreme, up to about half of the coefficient on 4g's model in TCGA, and a third in SCAN-B.
- **What remains.** The remainder is negative in every cell, but it has no interval.
- **MYC activity's coefficient is unaffected.**
- **The test** would take the next free script number, `E51` or later.

**How the metrics are counted.** Words are whitespace tokens containing a letter or digit. Figure citations in parentheses and reference numbers are not counted.

**Decisions for the author.**

1. **Supply the STOP 3 rulings** so that the text can be amended.
2. **The palette:** Fig. 1g's source, as drawn, or the measured values.
3. **The open choices above.**
4. **Whether to declare the H2 test.**
5. **The proposed legend additions** (E49 result note, 4.3), still marked as proposed below.
""")

q.append("{{< pagebreak >}}\n")
q.append("# Fig. 4g (E50)\n")
q.append(f"![](figures/2026-11-03_E50_bcl2_family_heatmap_level.png){{width={W_HEAT}}}\n")
q.append("## The text, as committed in E49's result note, Part 4 (not yet amended)\n")
for (ptitle, short, ss) in pt.PARAS:
    q.append(f"**{ptitle}.** " + " ".join(t for (_, t, _) in ss) + "\n")
q.append("## Legend draft\n")
q.append(f"*Verbatim from E50* ({nwords(LEG50['fig4g'])} words). {LEG50['fig4g']}\n")
q.append(f"*Proposed additions, not in the script's draft (E49 result note, 4.3)* "
         f"({nwords(pt.LEGEND_ADDITIONS)} words). {pt.LEGEND_ADDITIONS}\n")
q.append("STYLE.md does not yet cover legends (its section 9), so only their word counts are given.\n")
q.append("## STYLE.md metrics\n")
q.append(metrics_table() + "\n")
q.append(results_summary() + "\n")
q.append("No sentence is under 12 words or over 38. There are no em-dashes and no serial commas: "
         "each \", and\" joins two clauses. The clinical sentences C1 to C3 are the brief's own.\n")

q.append("{{< pagebreak >}}\n")
q.append("# The Extended Data panels\n")
q.append("Each is shown at final size, with its legend draft verbatim from the script that drew it. "
         "The letters are placeholders: a to c are E48's; e to h are E50's, which redraws E49's; i to k are E49's. "
         "E48's ED-d is superseded for the paper by ED-i, its extension.\n")
for k, (title, figs, leg, cite) in enumerate(ED):
    if k >= 2:  # ED-a and ED-b share the first page; every later panel starts its own
        q.append("{{< pagebreak >}}\n")
    q.append(f"## {title}\n")
    q.append(f"Cited as: {cite}.\n")
    for (f, w) in figs:
        q.append(f"![]({f}){{width={w}}}\n")
    q.append(f"*Legend draft* ({nwords(leg)} words). {leg}\n")

open(OUT, "w", encoding="utf-8").write("\n".join(q))
print("wrote", OUT)
