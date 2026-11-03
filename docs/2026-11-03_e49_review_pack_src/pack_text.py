# The E49 result note's Part 4 text: brief section 8 finalised against the
# verified numbers, with the amendment's two sentences (declaration 14.1) and
# the effector sentence. Single source for the result note and the review
# pack. Run directly to print the STYLE.md metrics.
import re, statistics, sys

XA, XB, XH, XJ, XK = ("Extended Data Fig. Xa", "Extended Data Fig. Xb",
                      "Extended Data Fig. Xh", "Extended Data Fig. Xj",
                      "Extended Data Fig. Xk")
# (title, short name, [(id, text, optional)])
PARAS = [
 ("Abstract: the human sentences", "Abstract", [
  ("A1", "Extending this relationship to 4,081 human breast tumours[8,9], an OXPHOS-associated BCL-2 family signature showed BCL-xL rising more steeply than PUMA.", False),
  ("A2", "MYC activity, once separated from proliferation, was associated with the same signature, whereas proliferation was associated with the opposite when OXPHOS was held fixed.", False),
 ]),
 ("Mouse Fig. 4 paragraph: last sentence (unchanged)", "Mouse Fig. 4, last sentence", [
  ("ML", "Consequently, established human tumours with high OXPHOS and MYC activity were predicted to carry raised pro-apoptotic signalling restrained by BCL-xL, distinct from the consequences of proliferation.", False),
 ]),
 ("The signature paragraph", "Signature paragraph", [
  ("S1", "In tumours, however, OXPHOS rises with both MYC activity and proliferation, so their associations can only be interpreted once the three are separated.", False),
  ("S2", "In TCGA-BRCA (n = 938) and SCAN-B (n = 3,143) we therefore modelled each BCL-2 family gene on all three together, adjusting for intrinsic subtype (Fig. 4g).", False),
  ("S2b", "Two of four MYC activity scores, including the widely used Hallmark set, correlated with proliferation too closely to be separated from it, so we used one containing few proliferation genes (" + XA + ").", False),
  ("S3", "Six of the twelve BH3-only and anti-apoptotic genes, the mitochondrial members most closely associated with OXPHOS, formed a signature (Fig. 4g).", False),
  ("S4", "BCL2L1, encoding BCL-xL, rose with the activators BBC3 and BID and the sensitisers BAD and BIK, whilst MCL1 fell.", False),
  ("S4b", "Outside the signature, BCL2L11, BMF and BCL2 also fell with OXPHOS in both cohorts (Fig. 4g).", False),
  ("S5", "Because the signature was defined by OXPHOS, that association is descriptive, whereas its associations with proliferation and MYC activity played no part in its definition.", False),
  ("S6", "Although proliferation correlated with OXPHOS (ρ = 0.23–0.42; " + XA + "), the signature did not follow it (" + XH + "), and at fixed OXPHOS and MYC activity higher proliferation went with lower signature expression (Fig. 4g).", False),
  ("S7", "On the score with few proliferation genes, and on a second, equally separable score, MYC activity went with the signature in the same direction as OXPHOS, whereas the other two scores showed this less (Fig. 4g and " + XJ + ").", False),
  ("S8", "OXPHOS and MYC activity therefore went with the signature in one direction, and proliferation in the other.", False),
  ("S9", "The effectors BAX, BAK1 and BOK, outside the signature, rose with OXPHOS and MYC activity in both cohorts, whilst proliferation went with lower BAX and BOK (Fig. 4g).", False),
 ]),
 ("The balance paragraph", "Balance paragraph", [
  ("B1", "Within the signature, BCL2L1 rose more steeply with OXPHOS than BBC3, encoding PUMA, in both cohorts (" + XB + ").", False),
  ("B2", "This is the reverse of the MYC-expressing gland, in which PUMA rose relative to BCL-xL with OXPHOS priority (Fig. 1i; " + XK + "), and resembles the mouse tumours that sustained high OXPHOS only alongside Bcl-xL (Fig. 3h,i).", False),
  ("B3", "BAD rose most steeply of all (Fig. 4g), but BAD neutralises BCL-xL only when unphosphorylated, and survival signals phosphorylate it and hold it away from BCL-xL[Zha].", False),
  ("B4", "Its expression therefore does not measure its activity, unlike that of PUMA, which is controlled mainly through transcription[Nakano].", False),
 ]),
 ("The MCL1 paragraph", "MCL1 paragraph", [
  ("M1", "MCL1 fell as BCL-xL rose, so BCL-xL was the only one of the two main guardians to rise alongside PUMA, which binds both[Chen].", False),
  ("M2", "MCL1's reported support of respiration in breast cancer stem cells depends on a mitochondrial function separable from its anti-apoptotic role[27].", False),
  ("M3", "This suggests that tumours with high OXPHOS rely on BCL-xL rather than MCL1, like persistent OXPHOS-high triple-negative breast cancer cells, which depend on BCL-xL and are insensitive to MCL1 inhibition[26].", False),
 ]),
 ("The clinical paragraph: edited sentences only", "Clinical, edited sentences", [
  ("C1", "Having separated OXPHOS from proliferation on the BCL-2 family, we next asked whether the two are also separable in their association with chemotherapy response.", False),
  ("C2", "Additional adjustment for MYC activity, and then for the mean expression of BCL2L1 and MCL1, did not materially alter the OXPHOS–response association (pooled odds ratios 0.87 and 0.88).", False),
  ("C3", "The transcriptional coupling of OXPHOS and proliferation therefore concealed a clinically relevant distinction.", False),
 ]),
 ("The closing paragraph: its human sentence", "Closing, human sentence", [
  ("Z1", "In human breast tumours, OXPHOS and MYC activity were associated with a BCL-2 family signature in which BCL-xL outpaced PUMA, whereas proliferation ran the other way, on this signature as on chemotherapy response.", False),
 ]),
 ("Methods notes (the signature's definition is E47 declaration 4a.7, item 5, unchanged and not re-measured here)", "Methods notes", [
  ("N1", "The effectors BAX, BAK1 and BOK were fitted for comparison with the groups of Fig. 1g and were not added to the signature.", False),
  ("N2", "Binding of BCL-xL and MCL1 by each BH3-only protein was taken from published peptide affinities and is shown only where the source states it[Chen].", False),
  ("N3", "This analysis was post hoc, and the signature's proliferation coefficient had been seen in an earlier analysis before it was declared; no reading rule was applied, and every coefficient is reported.", False),
 ]),
]

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

if __name__ == "__main__":
    for short, m in metrics_rows():
        print(f"{short:30s} n={m['sentences']:2d} words={m['words']:3d} mean={m['mean']:.1f} median={m['median']} range={m['min']}-{m['max']} <20={m['under20']} <12={m['under12']} >38={m['over38']} emdash={m['emdash']} each={m['each']}")
    print("\nserial-comma candidates:")
    for title, short, ss in PARAS:
        for (i, t, o) in ss:
            for mm in re.finditer(r",\s+and\b", t):
                print(" ", i, "...", t[max(0, mm.start() - 50): mm.end() + 10])
