# Writes docs/2026-11-03_e50_review_pack.qmd. It reads:
#   - the text, from its single source, ../2026-11-03_e49_review_pack_src/pack_text.py:
#     the author's STOP 3 rulings, section 3, as their addendum amends it, with
#     E47 4a.7, item 5, as the Methods definition; section 4's edits to the 4g
#     legend; the STYLE.md counter. pack_text first checks itself against the
#     rulings, the addendum and E47's declaration;
#   - this folder's legends.json, E50's $legend_draft verbatim (4g, ED-e to ED-h);
#   - E49's and E48's legends.json, verbatim (ED-i to ED-k; ED-a and ED-b);
#   - data/annotations/bh3_binding_code.csv, only to check legend edit 6's flag.
import sys, json, csv, re, importlib.util
S, OUT = sys.argv[1], sys.argv[2]
sys.dont_write_bytecode = True
spec = importlib.util.spec_from_file_location("pt", S + "/../2026-11-03_e49_review_pack_src/pack_text.py")
pt = importlib.util.module_from_spec(spec); spec.loader.exec_module(pt)
N_CHECKED = pt.check_against_rulings(S + "/..")
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
# Addendum 1: the model sets stay in the 4g legend.
if not ("(n = 938)" in LEG4G and "(n = 3,143)" in LEG4G):
    raise SystemExit("the 4g legend no longer gives the model sets, 938 and 3,143 (addendum 1)")

def nwords(s):
    return len([t for t in s.split() if any(c.isalnum() for c in t)])

def place(short):
    return "Methods" if short.startswith("Methods") else "Results" if short in pt.RESULTS_SHORTS else short

# Where each Extended Data panel is cited, read from the text and the legend.
CITE = re.compile(r"Extended Data Fig\. X([a-z](?:,[a-z])*)")
def cited_in():
    where = {}
    for (title, short, ss) in pt.PARAS:
        for (i, t, _) in ss:
            for m in CITE.finditer(t):
                for k in m.group(1).split(","):
                    where.setdefault(k, []).append((place(short), i))
    for m in CITE.finditer(LEG4G):
        for k in m.group(1).split(","):
            where.setdefault(k, []).append(("4g legend", "4g legend"))
    return where
WHERE = cited_in()
# The rulings' 1.4: a, b, h, j and k in the Results; g in the Methods; e, f and i from the 4g legend.
RULED = {"a": "Results", "b": "Results", "h": "Results", "j": "Results", "k": "Results",
         "g": "Methods", "e": "4g legend", "f": "4g legend", "i": "4g legend"}
for k, pl in RULED.items():
    if pl not in [p for (p, _) in WHERE.get(k, [])]:
        raise SystemExit(f"ED-{k}: not cited in the {pl}, as the rulings' 1.4 requires")
for k in WHERE:
    if k not in RULED:
        raise SystemExit(f"ED-{k} is cited, but the rulings' 1.4 keeps it out of the paper")

def cited_as(k):
    parts = []
    for pl in ("Results", "Methods", "4g legend"):
        ids = list(dict.fromkeys(i for (p, i) in WHERE[k] if p == pl))
        if not ids:
            continue
        if pl == "4g legend":
            parts.append("the 4g legend")
        elif pl == "Methods":
            parts.append("the Methods (" + " and ".join(ids) + ")")
        else:
            parts.append(("sentence " if len(ids) == 1 else "sentences ") + " and ".join(ids))
    return f"Extended Data Fig. X{k}, in " + ", and in ".join(parts)

def metrics_table():
    rows = ["| paragraph | sentences | words | mean | median | range | under 20 | under 12 | over 38 | em-dashes |",
            "|:-----------------------------|-------:|-----:|-----:|------:|-------:|--------:|-------:|-------:|---------:|"]
    for short, m in pt.metrics_rows():
        rows.append(f"| {short} | {m['sentences']} | {m['words']} | {m['mean']:.1f} | {m['median']:g} | "
                    f"{m['min']}–{m['max']} | {m['under20']} of {m['sentences']} | {m['under12']} | "
                    f"{m['over38']} | {m['emdash']} |")
    return "\n".join(rows)

R = pt.metrics(pt.results_sentences())
R_PCT = 100 * R["under20"] / R["sentences"]
ALL = pt.metrics([t for (_, _, ss) in pt.PARAS for (_, t, _) in ss])
if ALL["emdash"]:
    raise SystemExit("there is an em-dash")
OUT_OF_RANGE = pt.outside_range()
LO, HI = pt.RANGE
LONG = [(i, w) for (i, w, _) in OUT_OF_RANGE if w > HI]
SHORT = [(i, w) for (i, w, _) in OUT_OF_RANGE if w < LO]
OOR_PARAS = sorted(set(s for (_, _, s) in OUT_OF_RANGE))
fmt = lambda xs: " and ".join(f"{i} ({w} words)" for (i, w) in xs)
def oor_lines():
    if not OUT_OF_RANGE:
        return "- **Every sentence is inside STYLE.md's range, 12 to 38 words.**\n"
    s = [f"- **{len(OUT_OF_RANGE)} sentences are outside STYLE.md's range, 12 to 38 words, all in: "
         f"{', '.join(OOR_PARAS)}.** Each is E47's wording, kept verbatim, for the author to revise."]
    if LONG:
        s.append(f"  - **Over 38, to split:** {fmt(LONG)}.")
    if SHORT:
        s.append(f"  - **Under STYLE.md's hard floor of 12:** {fmt(SHORT)}. Splitting cannot fix these; joining can.")
    s.append("  - **The construction sentences, D8 and D9, are inside the range** "
             f"({len(pt.words(pt.sentence('D8')))} and {len(pt.words(pt.sentence('D9')))} words).")
    return "\n".join(s) + "\n"

# (letter, title, [(figure, width)], legend, note)
ED = [
  ("a", "ED-a (E48): the scores' correlations", [("figures/2026-11-03_E48_score_correlations.png", W_PAGE)], LEG48["ed_a"], ""),
  ("b", "ED-b (E48): *BCL2L1* against *BBC3*", [("figures/2026-11-03_E48_bcl2l1_against_bbc3.png", W_PAGE)], LEG48["ed_b"], ""),
  ("e", "ED-e (E50): Fig. 4g on OXPHOS priority", [("figures/2026-11-03_E50_bcl2_family_heatmap_priority.png", W_HEAT)], LEG50["ed_e"], ""),
  ("f", "ED-f (E50): Fig. 4g from the reference model", [("figures/2026-11-03_E50_bcl2_family_heatmap_reference_score.png", W_HEAT)], LEG50["ed_f"], ""),
  ("g", "ED-g (E50): the signature and its two variants", [("figures/2026-11-03_E50_signature_strips.png", W_PAGE)], LEG50["ed_g"],
   "**The author deletes the \"Proliferation minus OXPHOS\" column in Illustrator** (rulings 1.3; addendum 4). It contains the "
   "coefficient that defines the signature, and it remains in E49's tables. **ED-g keeps its ±0.9 scale.** No legend edit and "
   "no new script follow."),
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
title: "E50 review pack: Fig. 4g and the text, after the STOP 3 rulings and their addendum"
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

**What this is.** Fig. 4g as E50 draws it, at final size (87.6 mm wide). With it come the text and the 4g legend as the author's STOP 3 rulings and their addendum fix them, the STYLE.md metrics, and the nine Extended Data panels the rulings put in the paper.

- **The rulings:** `docs/2026-11-03_stop3_rulings_consolidated.md` (`cb9b7ce`), committed verbatim and alone.
- **The addendum:** `docs/2026-11-03_stop3_rulings_addendum.md` (`3f9ffbe`), committed verbatim and alone.
- **Nothing new is computed.** Every number drawn was stored by E49 (`4fd661b`) and drawn by E50, which was verified (`docs/2026-11-03_e50_result.md`).

**The text** (rulings section 3, as the addendum amends it).

- **No sentence is Claude Code's.**
  - The sentences are the rulings' own.
  - R2 carries the search's cohorts, from the addendum.
  - The Methods definition is E47 declaration 4a.7, item 5, verbatim, with the addendum's two changes.
- **{N_CHECKED} comparisons, made each time this pack is built, found no difference.** They check the text against the rulings, the addendum and E47's declaration.
- **The 4g legend is E50's draft with the rulings' section 4 edits**, items 1 to 5. Item 6 waits on the binding classes.
- **It supersedes E49's result note, Part 4, as committed** (`4fd661b`). The note now carries this text.

**The addendum's answers to the four observations.**

1. **R2 now gives the cohorts the search used,** n = 1,095 and 3,207. The model sets, 938 and 3,143, stay in the 4g legend, and the abstract's 4,081 stands.
2. **The brief's rule 5.4 now reads:** "{pt.RULE_5_4}" R7 meets it, and R4 stands unchanged.
3. **N4, N5 and N6 are dropped.** The Methods definition is E47 4a.7, item 5, verbatim, with two changes:
   - the construction sentence now says that selection was at fixed proliferation. It says the association with proliferation is in part expected by construction, and that the informative quantity is MYC activity's at fixed OXPHOS and proliferation (D8 and D9);
   - the variants sentence ends with "(Extended Data Fig. Xg)" (D10).

   N1 to N3 stay. N6, as ruled, misdescribed S5, whose rule uses MYC activity; E47's paragraph states it correctly (D10, D11).
4. **ED-g keeps its ±0.9 scale.** The author deletes the "Proliferation minus OXPHOS" column in Illustrator. No legend edit and no new script follow.

**For the author: the Methods definition, measured** (addendum 3: any sentence outside the range is reported for the author to split).

{oor_lines()}- **D10 ends "(omitting BID) (Extended Data Fig. Xg)."**, as addendum 3b reads. The two parentheses now sit side by side.

**The displays** (rulings section 1; addendum 4).

- **Accepted as E50 drew them:**
  - the palette from Fig. 1g's source, interpolated in CIELAB;
  - one ±0.6 scale shared by 4g, ED-e and ED-f;
  - the 1.6 mm gap before the Signature row;
  - white ink on 9 cells;
  - ED-g's ±0.9 scale.
- **The author's, in Illustrator, outside the repo's scripts.** No new display script is needed.
  - **Type sizes:** the group names and the colour bar's title at 6 pt, as in Fig. 1g.
  - **4g's annotations:** the binding code, the signature signs, the § marks and any protein names.
  - **ED-g:** the "Proliferation minus OXPHOS" column is deleted.
- **So the panels below are E50's PNG previews, before that work.** The PDF and SVG beside each in `docs/figures/` are the files for Illustrator, with live Helvetica text.
- **Extended Data in the paper:** ED-a, ED-b, ED-e, ED-f, ED-g, ED-h, ED-i, ED-j and ED-k, each cited.
  - **In the repo, out of the paper:** ED-c (`docs/figures/2026-11-03_E48_myc_score_ladder`), superseded by ED-i and ED-j, and S-1 (`2026-11-03_E48_sensitivities_tp53_composition`), kept for referee responses.
  - **ED-d** stays superseded by ED-i.

**H2** (rulings section 2).

- **Accepted.** The text says first that the signature was selected for tracking OXPHOS at fixed proliferation.
- **The signature's associations with OXPHOS and proliferation are what the selection intended.** The independent evidence is *BBC3* and *BCL2L1*, MYC activity, and the effectors.
- **The H2 test is not declared.** It is the answer if a referee asks how much of the signature's proliferation coefficient the selection explains.

**How the metrics are counted.** Words are whitespace tokens containing a letter or digit. Figure citations in parentheses and reference numbers are not counted.
""")

q.append("{{< pagebreak >}}\n")
q.append("# Fig. 4g (E50)\n")
q.append(f"![](figures/2026-11-03_E50_bcl2_family_heatmap_level.png){{width={W_HEAT}}}\n")
q.append("*E50's preview, before the author's Illustrator work: no annotations, and the group names and the bar's title at 5 pt, which become 6 pt.*\n")
q.append("## The text: the rulings' section 3, as the addendum amends it\n")
q.append("Sentence ids are for reference.\n")
for (ptitle, short, ss) in pt.PARAS:
    run_in = f"*{pt.RUN_IN[short]}* " if short in pt.RUN_IN else ""
    q.append(f"**{ptitle}.** {run_in}" + " ".join(f"**[{i}]** {t}" for (i, t, _) in ss) + "\n")
q.append("## The legend: E50's draft, with the rulings' section 4 edits\n")
q.append(f"*{nwords(LEG4G)} words; the draft has {nwords(LEG50['fig4g'])}.* {LEG4G}\n")
q.append("| item | the edit | status |\n|:-|:---------------------------|:----------|")
for (k, kind, old, new) in pt.LEGEND_4G_EDITS + [pt.LEGEND_4G_EDIT_6]:
    status = ("applied" if k != "6" or pt.BINDING_CLASSES_ENTERED else
              "waits on the binding classes (rulings section 5)")
    q.append(f"| {k} | {EDIT_WHAT[k]} | {status} |")
q.append("")
q.append("The model sets, n = 938 and 3,143, stay in this legend (addendum 1). STYLE.md does not yet cover legends (its section 9), so only the word count is given.\n")
q.append("## STYLE.md metrics\n")
q.append(metrics_table() + "\n")
q.append(f"**The four Results paragraphs together:** {R['sentences']} sentences, {R['words']} words, mean {R['mean']:.1f}, "
         f"median {R['median']:g}, range {R['min']}–{R['max']}, {R['under20']} of {R['sentences']} ({R_PCT:.0f}%) under 20 words. "
         f"These are the rulings' section 6 figures, reproduced. The addendum changes no word count there.\n")
q.append(f"- **Against STYLE.md's targets:** the range and the share under 20 are inside them (12 to 38; 35% to 45%). "
         f"The mean and the median sit just below them (22 to 25; 21 to 24).\n"
         + oor_lines() +
         f"- **There are no em-dashes and no serial commas.** Each \", and\" that the counter flags joins two clauses, two predicates or two groups, or opens a parenthesis.\n")

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
print("wrote", OUT, f"({N_CHECKED} comparisons with the rulings, the addendum and E47's declaration, all exact)")
