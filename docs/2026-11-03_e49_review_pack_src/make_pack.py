# Writes docs/2026-11-03_e49_review_pack.qmd from pack_text.py (the text, the
# proposed legend additions and the metrics) and two legends files, verbatim:
# this folder's legends.json (E49's $legend_draft) and E48's, for ED-a to ED-c.
import sys, json, importlib.util
S, OUT = sys.argv[1], sys.argv[2]
sys.dont_write_bytecode = True
spec = importlib.util.spec_from_file_location("pt", S + "/pack_text.py")
pt = importlib.util.module_from_spec(spec); spec.loader.exec_module(pt)
LEG = json.load(open(S + "/legends.json", encoding="utf-8"))
LEG48 = json.load(open(S + "/../2026-11-03_e48_review_pack_src/legends.json", encoding="utf-8"))

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
            f"**That share is below STYLE.md's 35% floor**; no sentence was shortened for the metric, and whether to "
            f"split one is the author's call." if pct < 35 else "")

TITLES = {short: title for (title, short, ss) in pt.PARAS}

VERSIONS = [
  ("a", "Version (a): unmarked", "figures/2026-11-03_E49_fig4g_heatmap_unmarked.png",
   [("Legend draft, verbatim from the script", LEG["fig4g_a"])]),
  ("b", "Version (b): a dot where the 95% CI excludes zero", "figures/2026-11-03_E49_fig4g_heatmap_marked.png",
   [("Legend draft, verbatim from the script", LEG["fig4g_a"]), ("then version (b)'s addition", LEG["fig4g_b"])]),
]

# (title, [figures], legend, cited as)
ED = [
  ("ED-a (E48): the scores' correlations", ["figures/2026-11-03_E48_score_correlations.png"], LEG48["ed_a"],
   "Extended Data Fig. Xa, in sentences S2b and S6"),
  ("ED-b (E48): *BCL2L1* against *BBC3*", ["figures/2026-11-03_E48_bcl2l1_against_bbc3.png"], LEG48["ed_b"],
   "Extended Data Fig. Xb, in sentence B1"),
  ("ED-c (E48): MYC activity on four scores, three genes", ["figures/2026-11-03_E48_myc_score_ladder.png"], LEG48["ed_c"],
   "not cited in the text. ED-f sits beside it (amendment item 3)"),
  ("ED-e: Fig. 4g on OXPHOS priority, versions (a) and (b)",
   ["figures/2026-11-03_E49_heatmap_priority_unmarked.png", "figures/2026-11-03_E49_heatmap_priority_marked.png"], LEG["ed_e"],
   "not cited in the text; the 4g legend points to it"),
  ("ED-f: Fig. 4g from the reference model, versions (a) and (b)",
   ["figures/2026-11-03_E49_heatmap_reference_score_unmarked.png", "figures/2026-11-03_E49_heatmap_reference_score_marked.png"], LEG["ed_f"],
   "not cited in the text; the 4g legend points to it"),
  ("ED-g: the signature and its two variants", ["figures/2026-11-03_E49_signature_strips.png"], LEG["ed_g"],
   "not cited in the text"),
  ("ED-h: raw correlations beside the adjusted coefficients, versions (a) and (b)",
   ["figures/2026-11-03_E49_raw_beside_adjusted_unmarked.png", "figures/2026-11-03_E49_raw_beside_adjusted_marked.png"], LEG["ed_h"],
   "Extended Data Fig. Xh, in sentence S6"),
  ("ED-i: all fifteen genes, three coefficients and the difference", ["figures/2026-11-03_E49_fifteen_genes_intervals.png"], LEG["ed_i"],
   "not cited in the text"),
  ("ED-j: MYC activity on four scores, the signatures and the effectors", ["figures/2026-11-03_E49_myc_ladder_signatures_effectors.png"], LEG["ed_j"],
   "Extended Data Fig. Xj, in sentence S7"),
  ("ED-k: PUMA relative to BCL-xL against OXPHOS priority", ["figures/2026-11-03_E49_puma_bclxl_against_priority.png"], LEG["ed_k"],
   "Extended Data Fig. Xk, in sentence B2"),
]

q = []
q.append("""---
title: "E49 review pack: Fig. 4g in two versions"
subtitle: "For the author at STOP 3. Exploratory and post hoc."
date: 2026-11-03
format:
  typst:
    papersize: a4
    margin:
      x: 16mm
      y: 15mm
    fontsize: 9pt
---

**What this is.** One section per version of Fig. 4g, each starting on a new page:

- **the panel at final size**, 178 mm wide;
- **the text**: brief section 8, finalised against the verified numbers, with the amendment's two sentences and the effector sentence;
- **the legend draft**, verbatim from the script, with the § note, followed by the proposed additions (brief 5.6), marked as such;
- **STYLE.md metrics** for every paragraph.

The Extended Data panels follow the two versions.

**What is the same in both versions.**

- **The values.** Every cell draws the same coefficient; the script asserts it. Version (b) adds a dot on each gene cell whose 95% CI excludes zero.
- **The text and the metrics.** Only the legend gains one line in version (b).
- **The choice carries** to ED-e, ED-f and ED-h, which are shown here in both versions.

**The record.**

- **Source:** `docs/2026-11-03_e49_result.md`. Part 4 is this text, and its section 4.1 says what changed from the brief's drafts and why.
- **The run:** 24 of 24 digests and 192 of 192 assertions; 33 of 33 PNGs and tables byte-identical to the dry run. Built 2026-11-03 13:05:54 GMT.

**One qualification found after the run** (result note, Part 2 H2; AD HOC arithmetic, not declared).

- **The issue.** Part of the signature's proliferation coefficient can be inherited from its selection rule, which chose genes by their correlation with OXPHOS at fixed proliferation.
- **How much.** If all of the signature's OXPHOS coefficient came from that selection, the inherited part would be up to about half of the coefficient on 4g's model in TCGA, and a third in SCAN-B.
- **What remains.** The remainder is negative in every cell, but it has no interval.
- **MYC activity's coefficient is unaffected.**
- **Where it bears on the text:** S5 ("played no part in its definition"), S6, S8, A2 and Z1. The test would be E50.

**How the metrics are counted.** Words are whitespace tokens containing a letter or digit. Figure citations in parentheses and reference numbers are not counted.

**Decisions for the author.**

1. **The version:** (a) unmarked, or (b) marked.
2. **The choices marked "the author's call"** in the result note, 4.1:
   - S3 and S4b, the six-gene sentence;
   - S5;
   - the abstract (A1);
   - the closing sentence (Z1);
   - "rely on" (M3);
   - "consequences".
3. **Whether to declare E50**, the test of H2, before the proliferation sentences are final.
4. **The proposed legend additions.**
5. **The binding code for BID, BIK and BMF**, which is drawn blank until the author supplies it.
""")

for (v, title, fig, legends) in VERSIONS:
    q.append("{{< pagebreak >}}\n")
    q.append(f"# {title}\n")
    q.append(f"![]({fig}){{width=178mm}}\n")
    q.append("## The text\n")
    for (ptitle, short, ss) in pt.PARAS:
        q.append(f"**{ptitle}.** " + " ".join(t for (_, t, _) in ss) + "\n")
    q.append("## Legend draft\n")
    for (lab, txt) in legends:
        q.append(f"*{lab}* ({nwords(txt)} words). {txt}\n")
    q.append(f"*Proposed additions, not in the script's draft (brief 5.6; result note 4.3)* "
             f"({nwords(pt.LEGEND_ADDITIONS)} words). {pt.LEGEND_ADDITIONS}\n")
    q.append("STYLE.md does not yet cover legends (its section 9), so only their word counts are given.\n")
    q.append("## STYLE.md metrics\n")
    q.append(metrics_table() + "\n")
    q.append(results_summary() + "\n")
    q.append("No sentence is under 12 words or over 38. There are no em-dashes and no serial commas: "
             "each \", and\" joins two clauses. The clinical sentences C1 to C3 are the brief's own.\n")

q.append("{{< pagebreak >}}\n")
q.append("# The Extended Data panels\n")
q.append("These are reproduced at final size, each with its legend draft verbatim from the script that drew it. "
         "The letters are placeholders: a to c are E48's, e to k are E49's. E48's ED-d is superseded for the paper by ED-i, its extension; E48's S-1 is kept as it was and is not repeated here.\n")
for k, (title, figs, leg, cite) in enumerate(ED):
    if k >= 2:  # ED-a and ED-b share the first page; every later panel starts its own
        q.append("{{< pagebreak >}}\n")
    q.append(f"## {title}\n")
    q.append(f"Cited as: {cite}.\n")
    for f in figs:
        q.append(f"![]({f}){{width=178mm}}\n")
    q.append(f"*Legend draft* ({nwords(leg)} words). {leg}\n")

open(OUT, "w", encoding="utf-8").write("\n".join(q))
print("wrote", OUT)
