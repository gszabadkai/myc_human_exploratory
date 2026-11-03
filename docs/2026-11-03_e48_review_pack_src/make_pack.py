# Writes docs/2026-11-03_e48_review_pack.qmd from pack_text.py (the text and
# its citations) and legends.json (the script's legend drafts, verbatim).
import sys, json, importlib.util
S, OUT = sys.argv[1], sys.argv[2]
spec = importlib.util.spec_from_file_location("pt", S + "/pack_text.py")
pt = importlib.util.module_from_spec(spec); spec.loader.exec_module(pt)
LEG = json.load(open(S + "/legends.json", encoding="utf-8"))

def render_marked(text, layout):
    """Citations that differ from layout 1 are set in bold."""
    vals = {}
    for k, v in pt.CIT[layout].items():
        vals[k] = v if v == pt.CIT["L1"][k] else "**" + v + "**"
    return text.format(**vals)

def prose(ss, layout):
    parts = []
    for (i, t, o) in ss:
        s = render_marked(t, layout)
        parts.append("[Optional: " + s + "]" if o else s)
    return " ".join(parts)

def metrics_table():
    rows = ["| paragraph | sentences | words | mean | median | range | under 20 | under 12 | em-dashes |",
            "|:------------------------------------------|--------:|-----:|-----:|------:|--------:|---------:|-------:|-------:|"]
    for title, ss in pt.PARAS:
        for opt in (False, True):
            if opt and not any(o for (_, _, o) in ss):
                continue
            m = pt.metrics([pt.render(t, "L1") for (i, t, o) in ss if (not o or opt)])
            name = pt.SHORT[title] + (", with optional" if opt else "")
            rows.append(f"| {name} | {m['sentences']} | {m['words']} | {m['mean']:.1f} | {m['median']:g} | {m['min']}–{m['max']} | {m['under20']} of {m['sentences']} | {m['under12']} | {m['emdash']} |")
    return "\n".join(rows)

def nwords(s):
    return len([t for t in s.split() if any(c.isalnum() for c in t)])

TITLES = {pt.PARAS[0][0]: "Abstract, sentence 8, with the *MCL1* sentence the addendum requires",
          pt.PARAS[1][0]: "Abstract, sentence 11",
          pt.PARAS[2][0]: "Mouse Fig. 4 paragraph, last sentence (unchanged)",
          pt.PARAS[3][0]: "The separation paragraph",
          pt.PARAS[4][0]: "The OXPHOS-side paragraph. The two co-amplification sentences wait on brief section 9 and are not shown",
          pt.PARAS[5][0]: "The clinical paragraph, edited sentences only",
          pt.PARAS[6][0]: "The closing paragraph, its human sentences",
          pt.PARAS[7][0]: "Methods"}

LAYOUTS = [
  ("L1", "Layout 1: two scores", "figures/2026-11-03_E48_separation_layout1_two_scores.png",
   [("Legend draft", LEG["layout1"])]),
  ("L2", "Layout 2: three scores", "figures/2026-11-03_E48_separation_layout2_three_scores.png",
   [("Legend draft, layout 1's legend", LEG["layout1"]), ("then layout 2's addition", LEG["layout2"])]),
  ("L3", "Layout 3: the whole argument in one panel", "figures/2026-11-03_E48_separation_layout3_with_outcomes.png",
   [("Legend draft, verbatim from the script", LEG["layout3"])]),
]

L3_NOTE = ("**A flag for layout 3.** The facet title reads \"BCL-2 family, OXPHOS level\", and this legend uses "
           "\"the BCL-2 family\" twice for *BBC3*, *BCL2L1* and *MCL1*. The section 8 addendum keeps that phrase out of the "
           "text for the three read genes. A legend variant that follows it: \"Top, *BBC3*, *BCL2L1* and *MCL1* as in "
           "layout 1, OXPHOS level only.\" and \"The top panel holds MYC activity fixed\". Changing the facet title is a "
           "re-render, which is a display-only task and the author's call.")

ED = [
  ("ED-a: the scores' correlations", "figures/2026-11-03_E48_score_correlations.png", LEG["ed_a"], "Extended Data Fig. Xa"),
  ("ED-b: *BCL2L1* against *BBC3*", "figures/2026-11-03_E48_bcl2l1_against_bbc3.png", LEG["ed_b"], "Extended Data Fig. Xb"),
  ("ED-c: MYC activity on four scores", "figures/2026-11-03_E48_myc_score_ladder.png", LEG["ed_c"], "Extended Data Fig. Xc"),
  ("ED-d: all twelve genes, three coefficients", "figures/2026-11-03_E48_twelve_genes_three_coefficients.png", LEG["ed_d"], "Extended Data Fig. Xd"),
  ("S-1: the sensitivities C1 and C4", "figures/2026-11-03_E48_sensitivities_tp53_composition.png", LEG["s1"], "none. Whether S-1 enters the paper is the author's call (declaration 12.3)"),
]

q = []
q.append("""---
title: "E48 review pack: Fig. 4g in three layouts"
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

**What this is.** One section per layout of Fig. 4g, each starting on a new page, as brief 6.5 asks:

- **the panel at final size**, 178 mm wide;
- **the text** from brief section 8, finalised against the verified numbers, with that layout's citations;
- **the legend draft**, descriptive only;
- **STYLE.md metrics** for every paragraph.

The Extended Data panels that the text cites follow the three layouts.

**What is the same in every layout.**

- **The declared reading** (declaration section 7):
  - *BBC3* OPPOSITE and *BCL2L1* OPPOSITE;
  - *MCL1* NOT SEPARATED;
  - the agreement check holds in 12 of 12 cells.
- **The text rules.** Outcome row 3 applies, with the section 8 addendum. C3 holds at *m3*, so the text says "with the other two held fixed".
- **The words.** They are identical in all three layouts. **Only the citations change, and a citation that differs from layout 1's is set in bold.**
- **The metrics.** They do not depend on the citations, so the table is the same on every page.

**The record.**

- **Source:** `docs/2026-11-03_e48_result.md`. Part 4 is this text, and its section 4.1 says what changed from the brief's drafts and why.
- **The run:** verified 33 of 33 digests and 175 of 175 assertions, built 2026-11-03 01:14:04 GMT.

**How the metrics are counted.** Words are whitespace tokens containing a letter or digit. Figure citations in parentheses and reference numbers are not counted.

**Decisions for the author.**

1. The layout. If it is layout 3, decide whether 4h and 4i stay in the main figure.
2. The choices marked "the author's call" in the result note, 4.1:
   - abstract sentence 11;
   - the clinical opening sentence;
   - the 9-word kept sentence;
   - "consequences" in the mouse paragraph;
   - the abstract and closing *MCL1* sentences.
3. The optional sentences.
4. Whether S-1 enters the paper.
""")

for (lay, title, fig, legends) in LAYOUTS:
    q.append("{{< pagebreak >}}\n")
    q.append(f"# {title}\n")
    q.append(f"![]({fig}){{width=178mm}}\n")
    q.append(f"## The text, with {title.split(':')[0].lower()}'s citations\n")
    for (ptitle, ss) in pt.PARAS:
        q.append(f"**{TITLES[ptitle]}.** " + prose(ss, lay) + "\n")
    q.append("## Legend draft\n")
    for (lab, txt) in legends:
        q.append(f"*{lab}* ({nwords(txt)} words). {txt}\n")
    if lay == "L3":
        q.append(L3_NOTE + "\n")
    q.append("STYLE.md does not yet cover legends (its section 9), so only their word counts are given.\n")
    q.append("## STYLE.md metrics\n")
    q.append(metrics_table() + "\n")
    q.append("No sentence is under 12 words. There are no em-dashes, no serial commas and no colons in the text.\n")

q.append("{{< pagebreak >}}\n")
q.append("# The Extended Data panels the text cites\n")
q.append("These are reproduced at final size, each with its legend draft verbatim from the script. Xa to Xd are placeholder letters.\n")
for (title, fig, leg, cite) in ED:
    q.append(f"## {title}\n")
    q.append(f"Cited as: {cite}.\n")
    q.append(f"![]({fig}){{width=178mm}}\n")
    q.append(f"*Legend draft* ({nwords(leg)} words). {leg}\n")

open(OUT, "w", encoding="utf-8").write("\n".join(q))
print("wrote", OUT)
