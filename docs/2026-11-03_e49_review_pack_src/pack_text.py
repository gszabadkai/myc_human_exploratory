# The paper's human text, as the author's STOP 3 rulings fix it
# (docs/2026-11-03_stop3_rulings_consolidated.md, cb9b7ce): section 3's
# sentences and section 4's edits to E50's 4g legend draft.
#
# - Every sentence is the rulings' own, verbatim, except N4, which is aligned to
#   E47 declaration 4a.7, item 5, as the rulings' section 3 asks (N4_AS_RULED
#   keeps the rulings' wording; ALIGNMENT says why).
# - This supersedes the E49 text that this file held at c117f8f (E49's result
#   note, Part 4, as committed at 4fd661b).
# - Single source for E49's result note, Part 4, and the E50 review pack.
# - Run directly to check the text against the rulings file and print the
#   STYLE.md metrics.
import re, statistics, sys

XA, XB, XG, XH, XK = ("Extended Data Fig. Xa", "Extended Data Fig. Xb",
                      "Extended Data Fig. Xg", "Extended Data Fig. Xh",
                      "Extended Data Fig. Xk")
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
  ("R2", "We therefore searched the twelve BH3-only and anti-apoptotic genes of the BCL-2 family in TCGA-BRCA (n = 938) and SCAN-B (n = 3,143) for those whose expression tracked OXPHOS at fixed proliferation.", False),
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
 ("Methods notes (the signature's definition is E47 declaration 4a.7, item 5, unchanged and not re-measured here)", "Methods notes", [
  ("N1", "The effectors BAX, BAK1 and BOK were fitted for comparison with the groups of Fig. 1g and were not added to the signature.", False),
  ("N2", "Binding of BCL-xL and MCL1 by each BH3-only protein was taken from published peptide affinities and is shown only where the source states it[Chen].", False),
  ("N3", "This analysis was post hoc, and the signature's proliferation coefficient had been seen in an earlier analysis before it was declared; no reading rule was applied, and every coefficient is reported.", False),
  ("N4", "The signature comprised those of the twelve BH3-only and anti-apoptotic genes listed in MitoCarta 3.0 whose partial Spearman correlation with OXPHOS at fixed proliferation, averaged over both cohorts, was at least 0.20 in magnitude.", False),
  ("N5", "Each was entered with the sign of that correlation, and MYC activity played no part in the selection.", False),
  ("N6", "Variants without BID, which met the threshold by the smallest margin, and with BMF, which met it but is not mitochondrial, gave similar associations (" + XG + ").", False),
 ]),
]

# The paragraphs STYLE.md's Results targets apply to.
RESULTS_SHORTS = ("Results 1, selection", "Results 2, balance", "Results 3, MCL1",
                  "Results 4, joint model")

# N4 as the rulings wrote it, and the one change made to align it with E47
# declaration 4a.7, item 5 (the rulings' section 3, "Methods notes").
N4_AS_RULED = ("The signature comprised the BCL-2 family members listed in MitoCarta 3.0 "
               "whose partial Spearman correlation with OXPHOS at fixed proliferation, averaged "
               "over both cohorts, was at least 0.20 in magnitude.")
ALIGNMENT = {
    "id": "N4",
    "ruled": "the BCL-2 family members listed in MitoCarta 3.0",
    "aligned": "those of the twelve BH3-only and anti-apoptotic genes listed in MitoCarta 3.0",
    "why": ("E47 4a.7 applies the rule among the twelve only (item 2: \"among the twelve, take "
            "every gene that is in Human MitoCarta 3.0\"; item 5: \"every gene that is a "
            "mitochondrial protein (Human MitoCarta 3.0: eleven of the twelve ...)\"). Human "
            "MitoCarta 3.0 also lists BAX, BAK1, BOK, BCL2L10, BCL2L13 and BNIP3, none of which "
            "was a candidate; read literally, the ruled wording would make the effectors "
            "eligible, against N1."),
}

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
# Item 6 waits on the binding classes (the rulings' section 5). Set this to True
# once the author has entered BID, BIK and BMF, with verified_by, in
# data/annotations/bh3_binding_code.csv; the E50 pack's make_pack.py stops if
# the flag and the table disagree.
BINDING_CLASSES_ENTERED = False
LEGEND_4G_EDIT_6 = ("6", "delete", "; blank, not yet verified", "")

def legend_4g(draft):
    """E50's 4g legend draft with the rulings' section 4 edits applied."""
    s = draft
    edits = LEGEND_4G_EDITS + ([LEGEND_4G_EDIT_6] if BINDING_CLASSES_ENTERED else [])
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

def check_against_rulings(path):
    """Every sentence and legend edit here, against the rulings file, verbatim
    (N4 as ruled). Returns the number of comparisons; stops on any mismatch."""
    txt = open(path, encoding="utf-8").read()
    s3 = txt.split("## 3. The text", 1)[1].split("## 4.", 1)[0]
    s4 = txt.split("## 4. The 4g legend", 1)[1].split("## 5.", 1)[0]
    quotes = [ln[2:].strip() for ln in s3.splitlines() if ln.startswith("> ")]
    j = lambda ids: " ".join(sentence(i) for i in ids)
    want = [j(["A1", "A2"]), j(["ML"]), j(["R1", "R2", "R3", "R4"]), j(["R5", "R6", "R7", "R8"]),
            j(["R9", "R10", "R11"]), j(["R12", "R13", "R14", "R15", "R16"]), j(["C1"]),
            j(["Z1", "Z2"]), " ".join([N4_AS_RULED, sentence("N5"), sentence("N6")])]
    if len(quotes) != len(want):
        raise SystemExit(f"rulings section 3: {len(quotes)} quoted paragraphs, expected {len(want)}")
    for k, (q, w) in enumerate(zip(quotes, want)):
        if q != w:
            raise SystemExit(f"rulings section 3, quoted paragraph {k + 1} differs:\n  {q}\n  {w}")
    unchanged = re.search(r'Its other two sentences are unchanged: "(.+?)" and "(.+?)"', s3)
    if not unchanged or unchanged.group(1) != sentence("C2") or unchanged.group(2) != sentence("C3"):
        raise SystemExit("rulings section 3: the clinical paragraph's unchanged sentences differ")
    if not sentence("N4").replace(ALIGNMENT["aligned"], ALIGNMENT["ruled"]) == N4_AS_RULED:
        raise SystemExit("N4 differs from the rulings' wording by more than the alignment")
    q4 = re.findall(r'"([^"]+)"', s4)
    n = len(want) + 3
    for (k, kind, old, new) in LEGEND_4G_EDITS + [LEGEND_4G_EDIT_6]:
        for part in [p for p in (old, new) if p]:
            if not any(part == q or (q.startswith("...") and part.endswith(q[3:])) for q in q4):
                raise SystemExit(f"legend edit {k}: not found verbatim in the rulings' section 4: {part}")
            n += 1
    return n

if __name__ == "__main__":
    here = sys.argv[0].rsplit("/", 1)[0] if "/" in sys.argv[0] else "."
    n = check_against_rulings(here + "/../2026-11-03_stop3_rulings_consolidated.md")
    print(f"against the rulings file: {n} comparisons, all verbatim (N4 as ruled; aligned here)\n")
    for short, m in metrics_rows():
        print(f"{short:30s} n={m['sentences']:2d} words={m['words']:3d} mean={m['mean']:.1f} median={m['median']} range={m['min']}-{m['max']} <20={m['under20']} <12={m['under12']} >38={m['over38']} emdash={m['emdash']} each={m['each']}")
    m = metrics(results_sentences())
    print(f"\nResults, four paragraphs: n={m['sentences']} words={m['words']} mean={m['mean']:.1f} "
          f"median={m['median']} range={m['min']}-{m['max']} <20={m['under20']} "
          f"({100 * m['under20'] / m['sentences']:.0f}%)")
    print(f"N4 as ruled: {len(words(N4_AS_RULED))} words; aligned: {len(words(sentence('N4')))}")
    print("\nserial-comma candidates:")
    for title, short, ss in PARAS:
        for (i, t, o) in ss:
            for mm in re.finditer(r",\s+and\b", t):
                print(" ", i, "...", t[max(0, mm.start() - 50): mm.end() + 10])
