# Part 4 text of the E48 result note, with per-layout citations and STYLE.md
# metrics. Single source for the result note and the review pack.
import re, statistics

XA, XB, XC, XD = ("Extended Data Fig. Xa", "Extended Data Fig. Xb",
                  "Extended Data Fig. Xc", "Extended Data Fig. Xd")
CIT = {
  "L1": {"c3a": "Fig. 4g and " + XA, "c3b": XA, "c6": "Fig. 4g and " + XD,
         "c7": "Fig. 4g", "c10": "Extended Data Fig. Xa,c"},
  "L2": {"c3a": "Fig. 4g and " + XA, "c3b": "Fig. 4g", "c6": "Fig. 4g and " + XD,
         "c7": "Fig. 4g", "c10": "Fig. 4g and Extended Data Fig. Xa,c"},
  "L3": {"c3a": "Fig. 4g and " + XA, "c3b": XA, "c6": "Fig. 4g and " + XD,
         "c7": "Fig. 4g and " + XD, "c10": "Extended Data Fig. Xa,c"},
}
# (id, text, optional)
PARAS = [
 ("Abstract: sentence 8, with the MCL1 sentence the addendum requires", [
  ("A8", "Extending this relationship to 4,081 human breast tumours[8,9], OXPHOS and proliferation had opposite associations with *BBC3* and *BCL2L1*, encoding PUMA and BCL-xL, when modelled together with MYC activity, OXPHOS tracking *BCL2L1* more steeply than *BBC3*.", False),
  ("A8m", "The association of proliferation with *MCL1* was negative in one cohort and positive in the other.", False),
  ("A8b", "The association of MYC activity with *BBC3* and *BCL2L1* depended on how much proliferation its score carried.", True),
 ]),
 ("Abstract: sentence 11", [
  ("A11", "Together, these findings establish an apoptotic constraint on MYC-driven transformation and reveal that, in human breast tumours, OXPHOS is linked to apoptotic regulator expression and treatment response separably from proliferation.", False),
 ]),
 ("Mouse Fig. 4 paragraph: last sentence (unchanged)", [
  ("M4L", "Consequently, human tumours sustaining high OXPHOS were predicted to carry raised BCL-xL alongside raised PUMA, distinguishable from the consequences of the proliferation that OXPHOS also supports.", False),
 ]),
 ("The separation paragraph", [
  ("P1", "In tumours, however, OXPHOS rises with both MYC activity and proliferation, so their associations can only be interpreted once the three are separated.", False),
  ("P2", "In TCGA-BRCA (n = 938) and SCAN-B (n = 3,143), OXPHOS was scored as subunit expression and as its priority within the mitochondrial programme, as in the mouse gland (Fig. 1e).", False),
  ("P3", "OXPHOS correlated moderately with proliferation and with the reference MYC activity score (Spearman ρ = 0.23–0.43; {c3a}).", False),
  ("P4", "The reference score tracked proliferation more closely, making MYC activity the axis least separable from proliferation (ρ = 0.61–0.70; {c3b}).", False),
  ("P5", "We therefore modelled *BBC3*, *BCL2L1* and *MCL1* expression on all three together, adjusting for intrinsic subtype in both cohorts and for tumour purity and leukocyte fraction in TCGA-BRCA.", False),
  ("P6", "With the other two held fixed, higher OXPHOS was associated with higher *BBC3* and *BCL2L1* and with lower *MCL1*, in both cohorts and on both measures ({c6}).", False),
  ("P7", "Proliferation had associations opposite to those of OXPHOS with *BBC3* and *BCL2L1*, being associated with lower expression of both at fixed OXPHOS and MYC activity ({c7}).", False),
  ("P8", "Proliferation was associated with lower *MCL1* expression in TCGA-BRCA and with higher *MCL1* expression in SCAN-B ({c7}).", False),
  ("P9", "Across the nine other BCL-2 family genes, OXPHOS and proliferation had opposite associations with three, *BAD*, *BCL2L11* and *BCL2*, in both cohorts and on both measures (" + XD + ").", True),
  ("P10", "The association of MYC activity with *BBC3* and *BCL2L1* accordingly depended on how much proliferation its score carried, the four scores examined correlating with proliferation at ρ = 0.41–0.84 ({c10}).", False),
  ("P11", "On the MYC activity score least entangled with proliferation, the association of MYC activity with each of the three genes had the same sign as that of OXPHOS (" + XC + ").", True),
  ("P12", "The correlation of OXPHOS with proliferation therefore concealed opposite associations of the two with *BBC3* and *BCL2L1*.", False),
 ]),
 ("The OXPHOS-side paragraph (the two co-amplification sentences wait on brief section 9 and are not shown)", [
  ("O1", "Of the genes that rose with OXPHOS, *BCL2L1* did so more steeply than *BBC3* in both cohorts and on both measures, with the other two held fixed (" + XB + ").", False),
  ("O2", "This resembled the mouse tumours in which high OXPHOS persisted only alongside Bcl-xL (Fig. 3h,i).", False),
  ("O5", "A similar state has been described in quiescent, chemotherapy-persistent triple-negative breast cancer cells, in which high OXPHOS accompanies dependence on BCL-xL[26].", False),
 ]),
 ("The clinical paragraph: edited sentences only", [
  ("C1", "Having separated OXPHOS from proliferation on *BBC3* and *BCL2L1*, we next asked whether the two are also separable in their association with chemotherapy response.", False),
  ("C2", "Additional adjustment for MYC activity, and then for the mean expression of *BCL2L1* and *MCL1*, did not materially alter the OXPHOS–response association (pooled odds ratios 0.87 and 0.88).", False),
  ("C3", "The transcriptional coupling of OXPHOS and proliferation therefore concealed a clinically relevant distinction.", False),
 ]),
 ("The closing paragraph: its human sentences", [
  ("Z1", "In human breast tumours, OXPHOS and proliferation, although correlated, had opposite associations with *BBC3*, *BCL2L1* and chemotherapy response, whilst the association of MYC activity with both genes depended on how much proliferation its score carried.", False),
  ("Z2", "The association of proliferation with *MCL1* was negative in one cohort and positive in the other.", False),
 ]),
 ("Methods", [
  ("M3", "The proliferation coefficients are a post hoc characterisation of a declared covariate, recovered from the saved model fits without refitting.", False),
  ("M4", "Each difference between the proliferation and OXPHOS coefficients was computed within one model as a linear contrast.", False),
  ("M4b", "Its variance was taken from the full coefficient covariance matrix, with a t-based 95% confidence interval on the residual degrees of freedom.", False),
  ("M7", "That *BCL2L1* rose more steeply than *BBC3* was assessed on 2,000 paired bootstrap resamples of tumours, refitting both gene models on each resample.", False),
  ("M5", "The reference MYC activity score correlated with proliferation at ρ = 0.61–0.70, a range in which the two are not cleanly separable, whereas OXPHOS correlated with proliferation at ρ = 0.23–0.42.", False),
  ("M6", "Variance inflation factors were 1.06–1.19 for OXPHOS, 1.35–1.60 for MYC activity and 1.66–1.72 for proliferation.", False),
  ("M8", "In TCGA-BRCA, the models were also refitted with TP53 mutation status added and, separately, without tumour purity and leukocyte fraction, as in SCAN-B.", True),
  ("M1", "This analysis was post hoc, and the directions of the proliferation coefficients could already be inferred from previously reported coefficients for MYC activity.", False),
  ("M2", "Their confidence intervals and the within-model differences were computed under a reading rule fixed before any of them was examined.", False),
 ]),
]

SHORT = dict(zip([p[0] for p in PARAS], ["Abstract, sentence 8", "Abstract, sentence 11", "Mouse Fig. 4, last sentence", "Separation paragraph", "OXPHOS-side paragraph", "Clinical, edited sentences", "Closing, human sentences", "Methods"]))

def render(text, layout):
    return text.format(**CIT[layout])

CITE_RE = re.compile(r"\((?:Fig\.|Extended Data Fig\.)[^)]*\)")
SEMI_RE = re.compile(r";\s*(?:Fig\.|Extended Data Fig\.)[^)]*(?=\))")
REF_RE = re.compile(r"\[\d[\d,–-]*\]")

def words(s):
    s = SEMI_RE.sub("", s)
    s = CITE_RE.sub("", s)
    s = REF_RE.sub("", s)
    s = s.replace("*", "")
    return [t for t in s.split() if any(ch.isalnum() for ch in t)]

def metrics(sents):
    n = [len(words(s)) for s in sents]
    if not n:
        return None
    return {"sentences": len(n), "words": sum(n), "mean": sum(n) / len(n),
            "median": statistics.median(n), "min": min(n), "max": max(n),
            "under20": sum(1 for k in n if k < 20), "under12": sum(1 for k in n if k < 12),
            "emdash": sum(s.count("—") for s in sents), "colons": sum(s.count(":") for s in sents),
            "each": n}

if __name__ == "__main__":
    for title, ss in PARAS:
        for opt in (False, True):
            use = [render(t, "L1") for (i, t, o) in ss if (not o or opt)]
            m = metrics(use)
            tag = "with optional" if opt else "core"
            if opt and not any(o for (_, _, o) in ss):
                continue
            print(f"{title[:40]:40s} [{tag:13s}] n={m['sentences']:2d} words={m['words']:3d} mean={m['mean']:.1f} median={m['median']} range={m['min']}-{m['max']} <20={m['under20']}/{m['sentences']} <12={m['under12']} emdash={m['emdash']} colons={m['colons']} each={m['each']}")
    # serial-comma scan: list ", <word> and" patterns for manual review
    print("\nserial-comma candidates (', X, and'):")
    for title, ss in PARAS:
        for (i, t, o) in ss:
            for mm in re.finditer(r",\s+and\b", t):
                print(" ", i, "...", t[max(0, mm.start() - 40): mm.end() + 10])

def note_block(layout="L1"):
    """Markdown for the result note's Part 4: each paragraph as a quote,
    optional sentences marked, then a metrics table."""
    out = []
    rows = []
    for title, ss in PARAS:
        out.append(f"**{title}.**\n")
        q = []
        for (i, t, o) in ss:
            s = render(t, layout)
            q.append(("[optional] " if o else "") + s)
        out.append("> " + "\n>\n> ".join(q) + "\n")
        for opt in (False, True):
            if opt and not any(o for (_, _, o) in ss):
                continue
            m = metrics([render(t, layout) for (i, t, o) in ss if (not o or opt)])
            rows.append((SHORT[title] + (", with optional" if opt else ""), m))
    out.append("**STYLE.md metrics**, counted on the text above. Words are whitespace tokens containing a letter or digit; figure citations in parentheses and reference numbers are not counted.\n")
    out.append("| paragraph | sentences | words | mean | median | range | under 20 | under 12 | em-dashes |")
    out.append("|---|---|---|---|---|---|---|---|---|")
    for name, m in rows:
        out.append(f"| {name} | {m['sentences']} | {m['words']} | {m['mean']:.1f} | {m['median']:g} | {m['min']}-{m['max']} | {m['under20']} of {m['sentences']} | {m['under12']} | {m['emdash']} |")
    return "\n".join(out) + "\n"

