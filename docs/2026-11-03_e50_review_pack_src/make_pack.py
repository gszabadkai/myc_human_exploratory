# Writes docs/2026-11-03_e50_review_pack.qmd. It reads:
#   - the text, from its single source, ../2026-11-03_e49_review_pack_src/pack_text.py:
#     the author's STOP 3 rulings, section 3 (N4 aligned), and section 4's edits to
#     the 4g legend, with the STYLE.md counter. pack_text checks itself against the
#     rulings file first;
#   - this folder's legends.json, E50's $legend_draft verbatim (4g, ED-e to ED-h);
#   - E49's and E48's legends.json, verbatim (ED-i to ED-k; ED-a and ED-b);
#   - data/annotations/bh3_binding_code.csv, only to check legend edit 6's flag.
import sys, json, csv, re, importlib.util
S, OUT = sys.argv[1], sys.argv[2]
sys.dont_write_bytecode = True
spec = importlib.util.spec_from_file_location("pt", S + "/../2026-11-03_e49_review_pack_src/pack_text.py")
pt = importlib.util.module_from_spec(spec); spec.loader.exec_module(pt)
N_CHECKED = pt.check_against_rulings(S + "/../2026-11-03_stop3_rulings_consolidated.md")
LEG50 = json.load(open(S + "/legends.json", encoding="utf-8"))
LEG49 = json.load(open(S + "/../2026-11-03_e49_review_pack_src/legends.json", encoding="utf-8"))
LEG48 = json.load(open(S + "/../2026-11-03_e48_review_pack_src/legends.json", encoding="utf-8"))

W_HEAT = "87.6mm"   # E50's heatmaps at final size (87.6 x 58.2 mm)
W_PAGE = "178mm"

# Legend edit 6 waits on the binding classes. Stop if the flag and the table disagree.
rows = {r["gene"]: r for r in csv.DictReader(open(S + "/../../data/annotations/bh3_binding_code.csv", encoding="utf-8"))}
entered = all(rows[g]["binds_BCL2L1"] and rows[g]["binds_MCL1"] and rows[g]["verified_by"] for g in ("BID", "BIK", "BMF"))
if entered != pt.BINDING_CLASSES_ENTERED:
    raise SystemExit(f"bh3_binding_code.csv has BID, BIK and BMF {'entered' if entered else 'not entered'}, "
                     f"but pack_text.BINDING_CLASSES_ENTERED is {pt.BINDING_CLASSES_ENTERED}")
LEG4G = pt.legend_4g(LEG50["fig4g"])

def nwords(s):
    return len([t for t in s.split() if any(c.isalnum() for c in t)])

# Where each Extended Data panel is cited, read from the text and the legend.
CITE = re.compile(r"Extended Data Fig\. X([a-z](?:,[a-z])*)")
def cited_in():
    where = {}
    for (title, short, ss) in pt.PARAS:
        for (i, t, _) in ss:
            for m in CITE.finditer(t):
                for k in m.group(1).split(","):
                    where.setdefault(k, []).append(("Methods" if short == "Methods notes" else
                                                    "Results" if short in pt.RESULTS_SHORTS else short, i))
    for m in CITE.finditer(LEG4G):
        for k in m.group(1).split(","):
            where.setdefault(k, []).append(("4g legend", "4g legend"))
    return where
WHERE = cited_in()
# The rulings' 1.4: a, b, h, j and k in the Results; g in the Methods; e, f and i from the 4g legend.
RULED = {"a": "Results", "b": "Results", "h": "Results", "j": "Results", "k": "Results",
         "g": "Methods", "e": "4g legend", "f": "4g legend", "i": "4g legend"}
for k, place in RULED.items():
    if place not in [p for (p, _) in WHERE.get(k, [])]:
        raise SystemExit(f"ED-{k}: not cited in the {place}, as the rulings' 1.4 requires")
for k in WHERE:
    if k not in RULED:
        raise SystemExit(f"ED-{k} is cited, but the rulings' 1.4 keeps it out of the paper")

def cited_as(k):
    ids = [i for (_, i) in WHERE[k]]
    sents = [i for i in ids if i != "4g legend"]
    parts = []
    if sents:
        noun = "Methods note " if all(i.startswith("N") for i in sents) else "sentence "
        if len(sents) > 1:
            noun = noun.replace("note ", "notes ").replace("sentence ", "sentences ")
        parts.append(noun + " and ".join(dict.fromkeys(sents)))
    if "4g legend" in ids:
        parts.append("the 4g legend")
    return f"Extended Data Fig. X{k}, in " + ", and in ".join(parts)

def metrics_table():
    rows = ["| paragraph | sentences | words | mean | median | range | under 20 | under 12 | em-dashes |",
            "|:------------------------------------------|--------:|-----:|-----:|------:|--------:|---------:|-------:|-------:|"]
    for short, m in pt.metrics_rows():
        rows.append(f"| {short} | {m['sentences']} | {m['words']} | {m['mean']:.1f} | {m['median']:g} | "
                    f"{m['min']}–{m['max']} | {m['under20']} of {m['sentences']} | {m['under12']} | {m['emdash']} |")
    return "\n".join(rows)

R = pt.metrics(pt.results_sentences())
R_PCT = 100 * R["under20"] / R["sentences"]
ALL = pt.metrics([t for (_, _, ss) in pt.PARAS for (_, t, _) in ss])
if ALL["under12"] or ALL["over38"] or ALL["emdash"]:
    raise SystemExit("a sentence is under 12 words or over 38, or there is an em-dash")
N4R = pt.N4_AS_RULED.replace(pt.ALIGNMENT["ruled"], "**" + pt.ALIGNMENT["ruled"] + "**")
N4A = pt.sentence("N4").replace(pt.ALIGNMENT["aligned"], "**" + pt.ALIGNMENT["aligned"] + "**")

# (letter, title, [(figure, width)], legend, note)
ED = [
  ("a", "ED-a (E48): the scores' correlations", [("figures/2026-11-03_E48_score_correlations.png", W_PAGE)], LEG48["ed_a"], ""),
  ("b", "ED-b (E48): *BCL2L1* against *BBC3*", [("figures/2026-11-03_E48_bcl2l1_against_bbc3.png", W_PAGE)], LEG48["ed_b"], ""),
  ("e", "ED-e (E50): Fig. 4g on OXPHOS priority", [("figures/2026-11-03_E50_bcl2_family_heatmap_priority.png", W_HEAT)], LEG50["ed_e"], ""),
  ("f", "ED-f (E50): Fig. 4g from the reference model", [("figures/2026-11-03_E50_bcl2_family_heatmap_reference_score.png", W_HEAT)], LEG50["ed_f"], ""),
  ("g", "ED-g (E50): the signature and its two variants", [("figures/2026-11-03_E50_signature_strips.png", W_PAGE)], LEG50["ed_g"],
   "**The author removes the \"Proliferation minus OXPHOS\" column in Illustrator** (rulings 1.3). It contains the coefficient that "
   "defines the signature, and it remains in E49's tables. The legend draft does not describe that column, so no legend edit follows."),
  ("h", "ED-h (E50): raw correlations beside the adjusted coefficients", [("figures/2026-11-03_E50_raw_beside_adjusted.png", W_PAGE)], LEG50["ed_h"], ""),
  ("i", "ED-i (E49): all fifteen genes, three coefficients and the difference", [("figures/2026-11-03_E49_fifteen_genes_intervals.png", W_PAGE)], LEG49["ed_i"], ""),
  ("j", "ED-j (E49): MYC activity on four scores, the signatures and the effectors", [("figures/2026-11-03_E49_myc_ladder_signatures_effectors.png", W_PAGE)], LEG49["ed_j"], ""),
  ("k", "ED-k (E49): PUMA relative to BCL-xL against OXPHOS priority", [("figures/2026-11-03_E49_puma_bclxl_against_priority.png", W_PAGE)], LEG49["ed_k"], ""),
]
assert [e[0] for e in ED] == sorted(RULED), "the pack's Extended Data must be the nine the rulings name"

EDIT_WHAT = {
  "1": "adds the sentence on cells and intervals, after the modelling sentence",
  "2": "the MYC activity sentence names its place among the four scores, and cites ED-f and ED-j",
  "3": "the signature sentence names the bottom row (\"Bottom, the six-gene signature\")",
  "4": "cites ED-e",
  "5": "adds the two sentences on patterns outside the text, replacing E49's proposed additions. The rulings do not place them; they are appended, where those additions stood",
  "6": "deletes \"; blank, not yet verified\" once BID, BIK and BMF are entered",
}

q = []
q.append(f"""---
title: "E50 review pack: Fig. 4g and the text, after the STOP 3 rulings"
subtitle: "For the author. Exploratory and post hoc."
date: 2026-11-03
format:
  typst:
    papersize: a4
    margin:
      x: 16mm
      y: 15mm
    fontsize: 9pt
---

**What this is.** Fig. 4g as E50 draws it, at final size (87.6 mm wide), with the text and the 4g legend as the author's STOP 3 rulings fix them, the STYLE.md metrics, and the nine Extended Data panels the rulings put in the paper.

- **The rulings:** `docs/2026-11-03_stop3_rulings_consolidated.md` (`cb9b7ce`), committed verbatim and alone.
- **Nothing new is computed.** Every number drawn was stored by E49 (`4fd661b`) and drawn by E50, which was verified (`docs/2026-11-03_e50_result.md`).

**The text** (rulings section 3).

- **Every sentence is the rulings' own, verbatim, except Methods note N4.** {N_CHECKED} comparisons against the committed file, made each time this pack is built, found no other difference.
- **N4 is aligned to E47 declaration 4a.7, item 5, as section 3 asks.**
  - *As ruled:* {N4R}
  - *Aligned:* {N4A}
  - *Why:* E47 4a.7 applies the rule among the twelve only. Human MitoCarta 3.0 also lists BAX, BAK1, BOK, BCL2L10, BCL2L13 and BNIP3, none of which was a candidate. Read literally, the ruled wording would make the effectors eligible, against N1.
- **The 4g legend is E50's draft with section 4's edits**, items 1 to 5. Item 6 waits on the binding classes.
- **It supersedes E49's result note, Part 4, as committed** (`4fd661b`). The note now carries this text.

**The displays** (rulings section 1).

- **Accepted as E50 drew them:**
  - the palette from Fig. 1g's source, interpolated in CIELAB;
  - one ±0.6 scale shared by 4g, ED-e and ED-f;
  - the 1.6 mm gap before the Signature row;
  - white ink on 9 cells.
- **The author's, in Illustrator, outside the repo's scripts.** No new display script is needed.
  - **Type sizes:** the group names and the colour bar's title at 6 pt, as in Fig. 1g.
  - **4g's annotations:** the binding code, the signature signs, the § marks and any protein names.
  - **ED-g:** the "Proliferation minus OXPHOS" column is removed.
- **So the panels below are E50's PNG previews, before that work.** The PDF and SVG beside each in `docs/figures/` are the files for Illustrator, with live Helvetica text.
- **Extended Data in the paper:** ED-a, ED-b, ED-e, ED-f, ED-g, ED-h, ED-i, ED-j and ED-k, each cited.
  - **In the repo, out of the paper:** ED-c (`docs/figures/2026-11-03_E48_myc_score_ladder`), superseded by ED-i and ED-j, and S-1 (`2026-11-03_E48_sensitivities_tp53_composition`), kept for referee responses.
  - **ED-d** stays superseded by ED-i.

**H2** (rulings section 2).

- **Accepted.** The text says first that the signature was selected for tracking OXPHOS at fixed proliferation.
- **The signature's associations with OXPHOS and proliferation are what the selection intended.** The independent evidence is *BBC3* and *BCL2L1*, MYC activity, and the effectors.
- **The H2 test is not declared.** It is the answer if a referee asks how much of the signature's proliferation coefficient the selection explains.

**For the author: four observations, not acted on.** E49's result note, section 4.1, gives the details of the first three.

1. **R2 gives the model sets, n = 938 and 3,143.** The 2026-09 selection table was computed on all 1,095 and 3,207 tumours (E10's `component_cor`).
2. **R4 names BAD among the genes that rose.** Its explanation is in R7, three sentences later, and brief 5.4 asks for the same sentence or the next. E49's committed text had the same gap.
3. **N4 to N6 restate part of E47 4a.7, item 5,** which stays the Methods' definition. If they replace that paragraph rather than sit beside it, three of its disclosures need to travel:
   - BID's stronger correlation with MYC activity (+0.30) and its split between cohorts (+0.33 and +0.09);
   - that S5's rule uses MYC activity, so S5 is selected against it;
   - the construction sentence, which E47 says may not be dropped.
4. **ED-g's colour scale, ±0.9, was set by the column the author removes.** Without that column, the amendment's rule (symmetric, rounded outward to 0.1) gives ±0.7: the largest remaining value is +0.62 (S7, TCGA, OXPHOS level, reference score).
   - The rulings accept three of E50's open choices and do not mention this one.
   - Kept at ±0.9, the remaining cells are paler than the rule would make them. Changing it means redrawing ED-g, and the rulings say no new display script is needed.

**How the metrics are counted.** Words are whitespace tokens containing a letter or digit. Figure citations in parentheses and reference numbers are not counted.
""")

q.append("{{< pagebreak >}}\n")
q.append("# Fig. 4g (E50)\n")
q.append(f"![](figures/2026-11-03_E50_bcl2_family_heatmap_level.png){{width={W_HEAT}}}\n")
q.append("*E50's preview, before the author's Illustrator work: no annotations, and the group names and the bar's title at 5 pt, which become 6 pt.*\n")
q.append("## The text: the rulings' section 3\n")
q.append("Sentence ids are for reference. N4 is aligned, as above.\n")
for (ptitle, short, ss) in pt.PARAS:
    q.append(f"**{ptitle}.** " + " ".join(f"**[{i}]** {t}" for (i, t, _) in ss) + "\n")
q.append("## The legend: E50's draft, with the rulings' section 4 edits\n")
q.append(f"*{nwords(LEG4G)} words; the draft has {nwords(LEG50['fig4g'])}.* {LEG4G}\n")
q.append("| item | the edit | status |\n|:-|:---------------------------|:----------|")
for (k, kind, old, new) in pt.LEGEND_4G_EDITS + [pt.LEGEND_4G_EDIT_6]:
    status = ("applied" if k != "6" or pt.BINDING_CLASSES_ENTERED else
              "waits on the binding classes (rulings section 5)")
    q.append(f"| {k} | {EDIT_WHAT[k]} | {status} |")
q.append("")
q.append("STYLE.md does not yet cover legends (its section 9), so only the word count is given.\n")
q.append("## STYLE.md metrics\n")
q.append(metrics_table() + "\n")
q.append(f"**The four Results paragraphs together:** {R['sentences']} sentences, {R['words']} words, mean {R['mean']:.1f}, "
         f"median {R['median']:g}, range {R['min']}–{R['max']}, {R['under20']} of {R['sentences']} ({R_PCT:.0f}%) under 20 words. "
         f"These are the rulings' section 6 figures, reproduced.\n")
q.append(f"- **Against STYLE.md's targets:** the range and the share under 20 are inside them (12 to 38; 35% to 45%). "
         f"The mean and the median sit just below them (22 to 25; 21 to 24).\n"
         f"- **No sentence is under 12 words or over 38, and there are no em-dashes.**\n"
         f"- **There are no serial commas.** Each \", and\" that the counter flags joins two clauses or two predicates.\n"
         f"- **The new Methods sentences, as ruled, are 30, 18 and 24 words; N4, aligned, is {len(pt.words(pt.sentence('N4')))}.**\n")

q.append("{{< pagebreak >}}\n")
q.append("# The Extended Data panels in the paper\n")
q.append("Each is shown at final size, with its legend draft verbatim from the script that drew it. "
         "The letters are placeholders: a and b are E48's; e to h are E50's, which redraws E49's; i to k are E49's. "
         "Where each is cited is read from the text and the legend above, and checked against the rulings' 1.4.\n")
for n, (k, title, figs, leg, note) in enumerate(ED):
    if n >= 2:  # ED-a and ED-b share the first page; every later panel starts its own
        q.append("{{< pagebreak >}}\n")
    q.append(f"## {title}\n")
    q.append(f"Cited as: {cited_as(k)}.\n")
    for (f, w) in figs:
        q.append(f"![]({f}){{width={w}}}\n")
    if note:
        q.append(note + "\n")
    q.append(f"*Legend draft* ({nwords(leg)} words). {leg}\n")

open(OUT, "w", encoding="utf-8").write("\n".join(q))
print("wrote", OUT, f"({N_CHECKED} comparisons with the rulings: N4 aligned, everything else verbatim)")
