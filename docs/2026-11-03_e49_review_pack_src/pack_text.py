# The paper's human text, as the author's STOP 3 rulings fix it
# (docs/2026-11-03_stop3_rulings_consolidated.md, cb9b7ce) and their addendum
# amends it (docs/2026-11-03_stop3_rulings_addendum.md, 3f9ffbe):
#   - the rulings' section 3 sentences, with R2's sample sizes from addendum 1;
#   - the Methods definition, E47 declaration 4a.7, item 5, verbatim, with the
#     addendum's two changes (3a, 3b); N4 to N6 are dropped (addendum 3);
#   - the rulings' section 4 edits to E50's 4g legend;
#   - the binding ruling (docs/2026-11-04_binding_classes_ruling.md, 66c205a):
#     legend edit 6 applies, and the legend and N2 cite Ku 2010 and Wang 2023.
# No sentence is Claude Code's. This supersedes the E49 text that this file held
# at c117f8f and the version of da548bd. Single source for E49's result note,
# Part 4, and the E50 review pack. Run directly to check the text against the
# three source files and print the STYLE.md metrics.
import re, statistics, sys

XA, XB, XG, XH, XK = ("Extended Data Fig. Xa", "Extended Data Fig. Xb",
                      "Extended Data Fig. Xg", "Extended Data Fig. Xh",
                      "Extended Data Fig. Xk")
RULINGS = "2026-11-03_stop3_rulings_consolidated.md"
ADDENDUM = "2026-11-03_stop3_rulings_addendum.md"
E47_DECL = "2026-10-22_E47_mutual_adjustment_declaration.md"
BINDING_RULING = "2026-11-04_binding_classes_ruling.md"

# The Methods definition's run-in heading, from E47 4a.7, item 5.
METHODS_HEADING = "BCL-2-family genes and signature."
RUN_IN = {"Methods definition": METHODS_HEADING}

# (title, short name, [(id, text, optional)])
PARAS = [
 ("Abstract: the human sentences", "Abstract", [
  ("A1", "Extending this relationship to 4,081 human breast tumours[8,9], a BCL-2 family signature that tracked OXPHOS at fixed proliferation showed BCL-xL rising more steeply than PUMA.", False),
  ("A2", "MYC activity, which played no part in defining it, went with the same signature when measured with scores least correlated with proliferation.", False),
 ]),
 ("Mouse Fig. 4 paragraph: last sentence", "Mouse Fig. 4, last sentence", [
  ("ML", "Consequently, established human tumours with high OXPHOS and MYC activity were predicted to carry raised pro-apoptotic signalling restrained by BCL-xL, a configuration distinct from that accompanying proliferation.", False),
 ]),
 ("Results, paragraph 1: the selection", "Results 1, selection", [
  ("R1", "In tumours, however, OXPHOS rises with both MYC activity and proliferation (" + XA + "), so testing this prediction required OXPHOS to be separated from proliferation.", False),
  ("R2", "We therefore searched the twelve BH3-only and anti-apoptotic genes of the BCL-2 family in TCGA-BRCA (n = 1,095) and SCAN-B (n = 3,207) for those whose expression tracked OXPHOS at fixed proliferation.", False),
  ("R3", "Six of the eleven mitochondrial members did, and together they formed a signature (Fig. 4g).", False),
  ("R4", "BCL2L1, encoding BCL-xL, rose with the activators BBC3 and BID and the sensitisers BAD and BIK, whilst MCL1 fell.", False),
 ]),
 ("Results, paragraph 2: the balance", "Results 2, balance", [
  ("R5", "Its structure fitted the prediction, BCL2L1 rising more steeply with OXPHOS than BBC3, encoding PUMA, in both cohorts (" + XB + ").", False),
  ("R6", "This is the reverse of the MYC-expressing gland, in which PUMA rose relative to BCL-xL with OXPHOS priority (Fig. 1i; " + XK + "), and resembles the mouse tumours that sustained high OXPHOS only alongside Bcl-xL (Fig. 3h,i).", False),
  ("R7", "BAD rose more steeply than BCL2L1 (Fig. 4g), but BAD neutralises BCL-xL only when unphosphorylated, and survival signals phosphorylate it and hold it away from BCL-xL[Zha].", False),
  ("R8", "Its expression therefore does not measure its activity, unlike that of PUMA, which is controlled mainly through transcription[Nakano].", False),
 ]),
 ("Results, paragraph 3: MCL1", "Results 3, MCL1", [
  ("R9", "MCL1 fell as BCL-xL rose, so BCL-xL was the only one of the two main guardians to rise alongside PUMA, which binds both[Chen].", False),
  ("R10", "MCL1's reported support of respiration in breast cancer stem cells depends on a mitochondrial function separable from its anti-apoptotic role[27].", False),
  ("R11", "This suggests that tumours with high OXPHOS rely on BCL-xL rather than MCL1, like persistent OXPHOS-high triple-negative breast cancer cells, which depend on BCL-xL and are insensitive to MCL1 inhibition[26].", False),
 ]),
 ("Results, paragraph 4: the joint model", "Results 4, joint model", [
  ("R12", "With OXPHOS, MYC activity and proliferation modelled together and intrinsic subtype adjusted for, the signature followed OXPHOS but not proliferation, as its selection intended (Fig. 4g).", False),
  ("R13", "So did BBC3 and BCL2L1, the two genes named by the mouse before any selection (" + XH + ").", False),
  ("R14", "MYC activity played no part in the selection, yet it went with the signature in the same direction as OXPHOS (Fig. 4g).", False),
  ("R15", "This held on the two MYC scores least correlated with proliferation, and was weaker on scores that track proliferation closely (Extended Data Fig. Xa,j).", False),
  ("R16", "The effectors BAX, BAK1 and BOK, outside the signature, also rose with OXPHOS and MYC activity in both cohorts (Fig. 4g).", False),
 ]),
 ("The clinical paragraph: the new opening sentence, then its two unchanged sentences", "Clinical paragraph", [
  ("C1", "Having separated OXPHOS from proliferation on PUMA and BCL-xL, we next asked whether the two are also separable in their association with chemotherapy response.", False),
  ("C2", "Additional adjustment for MYC activity, and then for the mean expression of BCL2L1 and MCL1, did not materially alter the OXPHOS–response association (pooled odds ratios 0.87 and 0.88).", False),
  ("C3", "The transcriptional coupling of OXPHOS and proliferation therefore concealed a clinically relevant distinction.", False),
 ]),
 ("The closing paragraph: its human sentences", "Closing, human sentences", [
  ("Z1", "In human breast tumours, BCL-xL rose more steeply than PUMA as OXPHOS rose, and MYC activity went with the OXPHOS-associated signature.", False),
  ("Z2", "Proliferation, although correlated with OXPHOS, was not part of this configuration and went the other way on chemotherapy response.", False),
 ]),
 ("Methods: the signature's definition (E47 declaration 4a.7, item 5, verbatim, with the addendum's two changes)", "Methods definition", [
  ("D1", "We examined twelve BCL-2-family genes that set the threshold for mitochondrial apoptosis: seven BH3-only proteins (BCL2L11/BIM, BID, BBC3/PUMA, BAD, BIK, BMF and PMAIP1/NOXA) and five anti-apoptotic proteins (BCL2, BCL2L1/BCL-xL, MCL1, BCL2L2/BCL-w and BCL2A1/BFL-1).", False),
  ("D2", "For each gene we computed the Spearman correlation of its expression with the OXPHOS subunit score (GSVA), partial on a proliferation score, in TCGA-BRCA and SCAN-B, and averaged the two cohorts.", False),
  ("D3", "The signature comprised every gene that is a mitochondrial protein (Human MitoCarta 3.0: eleven of the twelve, all outer-membrane; BMF is not listed) and whose mean correlation with OXPHOS was at least 0.20 in absolute value: BAD, BCL2L1, BIK, BBC3 and BID, which rise with OXPHOS (mean rho 0.21 to 0.50), and MCL1, which falls (mean rho -0.27).", False),
  ("D4", "Selection did not use correlation with MYC activity.", False),
  ("D5", "BID qualified narrowly (0.21) and is the one member more strongly correlated with MYC activity (0.30) than with OXPHOS.", False),
  ("D6", "Each gene was assigned a signature sign, +1 if it rises with OXPHOS and -1 if it falls, recording direction, not pro- or anti-apoptotic function.", False),
  ("D7", "A tumour's signature score was the mean of the six signed z-scores (log2-normalised expression, standardised across the cohort), rank-normalised within each cohort.", False),
  ("D8", "Because genes and signs were chosen by their association with OXPHOS at fixed proliferation in the same cohorts, the signature's association with OXPHOS, and in part with proliferation, is expected by construction.", False),
  ("D9", "The informative quantity is its association with MYC activity at fixed OXPHOS and proliferation.", False),
  ("D10", "All analyses were repeated with two rule-based variants: the same rule without the mitochondrial restriction (adding BMF, sign -1); and the five genes whose correlation with OXPHOS exceeded their correlation with MYC by at least 0.20 in at least one cohort (omitting BID) (" + XG + ").", False),
  ("D11", "The second variant is selected against association with MYC and is reported for transparency, not as an independent test.", False),
  ("D12", "All twelve genes are also reported individually.", False),
 ]),
 ("Methods notes (E49's three, kept)", "Methods notes", [
  ("N1", "The effectors BAX, BAK1 and BOK were fitted for comparison with the groups of Fig. 1g and were not added to the signature.", False),
  ("N2", "Binding of BCL-xL and MCL1 by each BH3-only protein was taken from published peptide affinities and is shown only where the source states it[Chen,Ku,Wang].", False),
  ("N3", "This analysis was post hoc, and the signature's proliferation coefficient had been seen in an earlier analysis before it was declared; no reading rule was applied, and every coefficient is reported.", False),
 ]),
]

# The paragraphs STYLE.md's Results targets apply to; its range, 12 to 38
# words, and its hard floor of 12 apply to every sentence.
RESULTS_SHORTS = ("Results 1, selection", "Results 2, balance", "Results 3, MCL1",
                  "Results 4, joint model")
RANGE = (12, 38)

# The record: what the rulings' section 3 said where the addendum changes it.
R2_AS_RULED = ("We therefore searched the twelve BH3-only and anti-apoptotic genes of the BCL-2 family "
               "in TCGA-BRCA (n = 938) and SCAN-B (n = 3,143) for those whose expression tracked OXPHOS "
               "at fixed proliferation.")
N4_AS_RULED = ("The signature comprised the BCL-2 family members listed in MitoCarta 3.0 "
               "whose partial Spearman correlation with OXPHOS at fixed proliferation, averaged "
               "over both cohorts, was at least 0.20 in magnitude.")
N5_AS_RULED = "Each was entered with the sign of that correlation, and MYC activity played no part in the selection."
N6_AS_RULED = ("Variants without BID, which met the threshold by the smallest margin, and with BMF, which "
               "met it but is not mitochondrial, gave similar associations (" + XG + ").")
N2_AS_RULED = ("Binding of BCL-xL and MCL1 by each BH3-only protein was taken from published "
               "peptide affinities and is shown only where the source states it[Chen].")
# The brief's rule 5.4, as the addendum (item 2) rewords it.
RULE_5_4 = ("BAD is reported with its explanation in the sentence that sets its rise against "
            "BCL-xL's, or the next.")

# The 4g legend: the rulings' section 4, as edits to E50's draft
# ($legend_draft$fig4g; docs/2026-11-03_e50_review_pack_src/legends.json).
# (item, kind, anchor in the draft, text). Each anchor must occur exactly once.
LEGEND_4G_EDITS = [
 ("1", "after", "every variable is rank-normalised within each cohort.",
  "Cells give the coefficient, rounded to one decimal; intervals are in Extended Data Fig. Xi."),
 ("2", "replace", "MYC activity is scored on MYC_UP.V1_UP, a score with few proliferation genes; "
  "the other MYC activity scores are shown in Extended Data.",
  "MYC activity is scored on MYC_UP.V1_UP, one of the two MYC activity scores least correlated "
  "with proliferation; the reference score is used in Extended Data Fig. Xf, and all four in "
  "Extended Data Fig. Xj."),
 ("3", "replace", "The signature was selected for its association with OXPHOS at fixed proliferation (Methods).",
  "Bottom, the six-gene signature, selected for its association with OXPHOS at fixed proliferation (Methods)."),
 ("4", "replace", "OXPHOS priority is shown in Extended Data.",
  "OXPHOS priority is shown in Extended Data Fig. Xe."),
 # Item 5 replaces E49's proposed additions (result note 4.3). The rulings do
 # not place it; it is appended, where those additions stood.
 ("5", "append", None,
  "Outside the signature, BCL2L11, BMF and BCL2 also fell with OXPHOS. At fixed OXPHOS and MYC "
  "activity, proliferation went with lower BAX and BOK and higher PMAIP1 and BCL2L11 in both "
  "cohorts; its association with BAK1 differed in sign between cohorts, and MCL1 was higher only "
  "in SCAN-B."),
]
# Item 6 waited on the binding classes (the rulings' section 5). The author
# ruled them on 2026-11-04 and they are entered in
# data/annotations/bh3_binding_code.csv; the E50 pack's make_pack.py stops if
# this flag and the table disagree.
BINDING_CLASSES_ENTERED = True
LEGEND_4G_EDIT_6 = ("6", "delete", "; blank, not yet verified", "")
# Item 7, from the binding ruling (question 2): the legend cites all three papers.
LEGEND_4G_EDIT_7 = ("7", "replace", "from Chen et al. (2005):",
                    "from Chen et al. (2005), Ku et al. (2010) and Wang et al. (2023):")

def legend_4g(draft):
    """E50's 4g legend draft with the rulings' section 4 edits applied."""
    s = draft
    edits = LEGEND_4G_EDITS + ([LEGEND_4G_EDIT_6] if BINDING_CLASSES_ENTERED else []) + [LEGEND_4G_EDIT_7]
    for (k, kind, old, new) in edits:
        if kind == "append":
            s = s + " " + new
            continue
        n = s.count(old)
        if n != 1:
            raise SystemExit(f"legend edit {k}: its anchor occurs {n} times in the draft, not once")
        s = s.replace(old, old + " " + new) if kind == "after" else s.replace(old, new)
    return s

CITE_RE = re.compile(r"\((?:Fig\.|Extended Data Fig\.)[^)]*\)")
SEMI_RE = re.compile(r";\s*(?:Fig\.|Extended Data Fig\.)[^)]*(?=\))")
REF_RE = re.compile(r"\[[^\]]*\]")

def words(s):
    s = SEMI_RE.sub("", s)
    s = CITE_RE.sub("", s)
    s = REF_RE.sub("", s)
    s = s.replace("*", "")
    return [t for t in s.split() if any(ch.isalnum() for ch in t)]

def metrics(sents):
    n = [len(words(s)) for s in sents]
    return {"sentences": len(n), "words": sum(n), "mean": sum(n) / len(n),
            "median": statistics.median(n), "min": min(n), "max": max(n),
            "under20": sum(1 for k in n if k < 20), "under12": sum(1 for k in n if k < 12),
            "over38": sum(1 for k in n if k > 38),
            "emdash": sum(s.count("—") for s in sents), "each": n}

def metrics_rows():
    rows = []
    for title, short, ss in PARAS:
        m = metrics([t for (_, t, _) in ss])
        rows.append((short, m))
    return rows

def results_sentences():
    return [t for (title, short, ss) in PARAS if short in RESULTS_SHORTS for (_, t, _) in ss]

def sentence(i):
    return {k: t for (_, _, ss) in PARAS for (k, t, _) in ss}[i]

def outside_range():
    """(id, words, paragraph) for every sentence outside STYLE.md's 12 to 38."""
    lo, hi = RANGE
    return [(i, len(words(t)), short) for (_, short, ss) in PARAS for (i, t, _) in ss
            if not lo <= len(words(t)) <= hi]

def _e47_item5(path):
    """E47 declaration 4a.7, item 5's Methods paragraph, its lines joined."""
    blk = open(path, encoding="utf-8").read().split(
        "**5. Methods text, revised. This SUPERSEDES the paragraph in 4a.5.**", 1)[1]
    q = []
    for ln in blk.splitlines()[1:]:
        if ln.startswith(">"):
            q.append(ln[1:].strip())
        elif q:
            break
    return " ".join(q)

def check_against_rulings(docs):
    """Every sentence and legend edit here, against the rulings, their addendum
    and E47's declaration, in the directory `docs`. Returns the number of
    comparisons; stops on any mismatch."""
    n = 0
    def same(a, b, what):
        nonlocal n
        if a != b:
            raise SystemExit(f"{what} differs:\n  {a}\n  {b}")
        n += 1
    j = lambda ids: " ".join(sentence(i) for i in ids)
    # 1. The rulings' section 3, as ruled: R2 and N4 to N6 from the record.
    txt = open(f"{docs}/{RULINGS}", encoding="utf-8").read()
    s3 = txt.split("## 3. The text", 1)[1].split("## 4.", 1)[0]
    s4 = txt.split("## 4. The 4g legend", 1)[1].split("## 5.", 1)[0]
    quotes = [ln[2:].strip() for ln in s3.splitlines() if ln.startswith("> ")]
    want = [j(["A1", "A2"]), j(["ML"]), " ".join([sentence("R1"), R2_AS_RULED, j(["R3", "R4"])]),
            j(["R5", "R6", "R7", "R8"]), j(["R9", "R10", "R11"]), j(["R12", "R13", "R14", "R15", "R16"]),
            j(["C1"]), j(["Z1", "Z2"]), " ".join([N4_AS_RULED, N5_AS_RULED, N6_AS_RULED])]
    same(len(quotes), len(want), "the number of quoted paragraphs in the rulings' section 3")
    for k, (q, w) in enumerate(zip(quotes, want)):
        same(q, w, f"the rulings' section 3, quoted paragraph {k + 1}")
    m = re.search(r'Its other two sentences are unchanged: "(.+?)" and "(.+?)"', s3)
    same((m.group(1), m.group(2)) if m else None, (sentence("C2"), sentence("C3")),
         "the clinical paragraph's two unchanged sentences")
    # 2. The rulings' section 4: every legend edit, verbatim.
    q4 = re.findall(r'"([^"]+)"', s4)
    for (k, kind, old, new) in LEGEND_4G_EDITS + [LEGEND_4G_EDIT_6]:
        for part in [p for p in (old, new) if p]:
            same(any(part == q or (q.startswith("...") and part.endswith(q[3:])) for q in q4), True,
                 f"legend edit {k}, found verbatim in the rulings' section 4")
    # 3. The addendum, read from its verbatim section.
    a = open(f"{docs}/{ADDENDUM}", encoding="utf-8").read().split("## 2. The addendum, verbatim", 1)[1]
    av = "\n".join(ln[2:] if ln.startswith("> ") else ln[1:] for ln in a.splitlines() if ln.startswith(">"))
    qa = re.findall(r'"([^"]+)"', av)
    pick = lambda start: [q for q in qa if q.startswith(start)]
    same(pick("We therefore searched"), [sentence("R2")], "R2, against addendum 1")
    same(R2_AS_RULED.replace("(n = 938)", "(n = 1,095)").replace("(n = 3,143)", "(n = 3,207)"),
         sentence("R2"), "R2 against the rulings' R2, apart from the sample sizes")
    same("4,081" in sentence("A1"), True, "the abstract's 4,081 (addendum 1)")
    same(pick("BAD is reported"), [RULE_5_4], "the reworded rule 5.4 (addendum 2)")
    old_c = pick("Because genes and signs were chosen by their association with OXPHOS in")
    new_c = pick("Because genes and signs were chosen by their association with OXPHOS at")
    same(new_c, [j(["D8", "D9"])], "the construction sentence, against addendum 3a")
    same(pick("(Extended Data Fig. Xg)"), ["(" + XG + ")"], "the variants' citation, addendum 3b")
    ids = [k for (_, _, ss) in PARAS for (k, _, _) in ss]
    same([k for k in ("N1", "N2", "N3", "N4", "N5", "N6") if k in ids], ["N1", "N2", "N3"],
         "N1 to N3 kept and N4 to N6 dropped (addendum 3)")
    # 4. E47 4a.7, item 5, verbatim, with the addendum's two changes.
    para = _e47_item5(f"{docs}/{E47_DECL}")
    head = "**" + METHODS_HEADING + "** "
    same(para.startswith(head), True, "the Methods definition's heading, against E47 4a.7, item 5")
    body = para[len(head):]
    same(len(old_c) == 1 and body.count(old_c[0]) == 1, True,
         "addendum 3a's quoted sentence, found once in E47 4a.7, item 5")
    same(body.count("(omitting BID)."), 1, "the end of E47's variants sentence")
    body = body.replace(old_c[0], new_c[0]).replace("(omitting BID).", "(omitting BID) (" + XG + ").")
    same(body, j(["D" + str(k) for k in range(1, 13)]), "the Methods definition, against E47 4a.7, item 5")
    # 5. The binding ruling: legend edit 7, N2's citations and the MCL1 call.
    b = open(f"{docs}/{BINDING_RULING}", encoding="utf-8").read().split("## 3. The author's answers, verbatim", 1)[1]
    k, kind, old, new = LEGEND_4G_EDIT_7
    same(f"The legend reads '{new.rstrip(':')}'" in b, True, "legend edit 7, against the binding ruling")
    same("and N2 cites all three" in b, True, "N2's three citations, against the binding ruling")
    same(sentence("N2").replace("[Chen,Ku,Wang]", "[Chen]"), N2_AS_RULED, "N2, apart from its citations")
    same("> Does not bind (open): binds_MCL1 = FALSE for BIK and BMF" in b, True,
         "the MCL1 call for BIK and BMF, in the binding ruling")
    return n

if __name__ == "__main__":
    here = sys.argv[0].rsplit("/", 1)[0] if "/" in sys.argv[0] else "."
    n = check_against_rulings(here + "/..")
    print(f"against the rulings, the addendum, E47 4a.7, item 5, and the binding ruling: {n} comparisons, all exact\n")
    for short, m in metrics_rows():
        print(f"{short:30s} n={m['sentences']:2d} words={m['words']:3d} mean={m['mean']:.1f} median={m['median']} range={m['min']}-{m['max']} <20={m['under20']} <12={m['under12']} >38={m['over38']} emdash={m['emdash']} each={m['each']}")
    m = metrics(results_sentences())
    print(f"\nResults, four paragraphs: n={m['sentences']} words={m['words']} mean={m['mean']:.1f} "
          f"median={m['median']} range={m['min']}-{m['max']} <20={m['under20']} "
          f"({100 * m['under20'] / m['sentences']:.0f}%)")
    print("outside STYLE.md's range (12 to 38):", outside_range())
    print("\nserial-comma candidates:")
    for title, short, ss in PARAS:
        for (i, t, o) in ss:
            for mm in re.finditer(r",\s+and\b", t):
                print(" ", i, "...", t[max(0, mm.start() - 50): mm.end() + 10])
