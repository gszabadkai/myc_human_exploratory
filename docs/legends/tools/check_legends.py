# Checks the draft legends against the author's legend rules.
# For each figure: every panel letter present and in order; every data panel
# carries n, a centre value, error bars and a test; P written with an italic P;
# Extended Data figures have no title; no interpretive verb anywhere; ASCII-safe
# characters only where the Doc would mangle them. Prints a table and exits 1 if
# anything fails. Reads only the draft file.
import re, sys

PATH = sys.argv[1] if len(sys.argv) > 1 else "docs/legends/2026-11-06_legends_draft.md"
text = open(PATH, encoding="utf-8").read()

# The panels each figure must have, from the assembled PDFs.
WANT = {
    "Fig. 2": list("abcdefghij"),
    "Extended Data Fig. 4": list("abcdefg"),
    "Fig. 3": list("abcdefghij"),
    "Fig. 4": list("abcdefghij"),
    "Extended Data Fig. 6": list("abcd"),
}
# Panels that draw no data of their own: blots and images need no n or test.
NO_STATS = {
    ("Fig. 2", "a"), ("Fig. 2", "c"), ("Fig. 2", "f"),
    ("Fig. 3", "b"), ("Fig. 3", "d"),
}
# Panels whose statistics are given in the figure's closing sentence.
SHARED = {
    ("Fig. 2", "g"), ("Fig. 2", "h"), ("Fig. 2", "i"),
    ("Fig. 4", "a"), ("Fig. 4", "c"), ("Fig. 4", "e"),
    ("Extended Data Fig. 4", "a"), ("Extended Data Fig. 4", "b"),
    ("Extended Data Fig. 4", "c"), ("Extended Data Fig. 4", "d"),
    ("Extended Data Fig. 4", "e"), ("Extended Data Fig. 4", "f"),
    ("Extended Data Fig. 4", "g"),
    ("Extended Data Fig. 6", "a"), ("Extended Data Fig. 6", "c"),
}
INTERPRETIVE = [
    r"\bsuggest(s|ing|ed)?\b", r"\bdemonstrat(e|es|ing|ed)\b", r"\bshow(s|ing|ed)? that\b",
    r"\bindicat(e|es|ing|ed)\b", r"\bconfirm(s|ing|ed)?\b", r"\breveal(s|ing|ed)?\b",
    r"\bsignificantly\b", r"\bimportantly\b", r"\bconsistent with\b", r"\bsuppress(es|ed|ing)?\b",
    r"\bdrives?\b", r"\bcauses?\b", r"\bbecause\b", r"\btherefore\b", r"\bthus\b",
]

sections = {}
for m in re.finditer(r"(?m)^## (Fig\. \d|Extended Data Fig\. \d)\s*$", text):
    name = m.group(1)
    end = text.find("\n---", m.end())
    sections[name] = text[m.end():end if end > 0 else len(text)]

rows, failed = [], 0


def check(label, ok, detail=""):
    global failed
    rows.append((label, "pass" if ok else "FAIL", detail))
    if not ok:
        failed += 1


for fig, letters in WANT.items():
    body = sections.get(fig, "")
    check(f"{fig}: section present", bool(body))
    if not body:
        continue
    found = re.findall(r"\*\*([a-z])\*\*,", body)
    seen, order = [], True
    for f in found:
        if f not in seen:
            seen.append(f)
    check(f"{fig}: every panel present, in order", seen == letters,
          f"found {''.join(seen)} want {''.join(letters)}")
    # A figure title only for main figures.
    has_title = bool(re.search(r"(?m)^\*\*Fig\. \d\.", body))
    check(f"{fig}: title present only if a main figure",
          has_title == fig.startswith("Fig."))
    # Split the body at each panel letter so each panel can be checked alone.
    parts = re.split(r"(?<=\.\s)(?=\*\*[a-z]\*\*,)|\A(?=\*\*[a-z]\*\*,)", body)
    for p in parts:
        m = re.match(r"\*\*([a-z])\*\*,", p)
        if not m:
            continue
        L = m.group(1)
        if (fig, L) in NO_STATS or (fig, L) in SHARED:
            continue
        p = re.sub(r"\s+", " ", p)
        has_n = bool(re.search(r"\*n\*", p))
        has_test = bool(re.search(r"ANOVA|t test|\*t\* test|Welch|meta-analysis|mixed-effects|"
                                  r"linear model|logistic model|Spearman|bootstrap", p))
        check(f"{fig} {L}: gives n", has_n, p[:60].replace("\n", " "))
        check(f"{fig} {L}: names a test", has_test, p[:60].replace("\n", " "))
    # Centre and error bars: in each panel or in the figure's closing sentence.
    has_centre = bool(re.search(r"mean|median", body))
    has_err = bool(re.search(r"standard deviation|s\.e\.m\.|confidence interval", body))
    check(f"{fig}: names the centre value", has_centre)
    check(f"{fig}: names the error bars", has_err)

# P is italic wherever it is a P value, and q likewise.
bad_p = [m.start() for m in re.finditer(r"(?<![*\w(])P(?=\s*[=<,;]|\s*values)", text)]
check("every P value is italic", not bad_p, f"{len(bad_p)} upright")
# No interpretation anywhere in the legends.
hits = []
for pat in INTERPRETIVE:
    bodies = "\n".join(re.sub(r"(?m)^\*\*Fig\. \d\..*$", "", v) for v in sections.values())
    for m in re.finditer(pat, bodies, re.I):
        hits.append(m.group(0))
check("no interpretive verb in any legend", not hits, ", ".join(sorted(set(hits))[:6]))
# Every flag referenced exists in the inventory.
inv = PATH.replace("legends_draft", "panel_inventory_and_stats")
try:
    invtext = open(inv, encoding="utf-8").read()
    used = sorted(set(re.findall(r"\[(F\d+)[,\]]", text)) | set(re.findall(r"\[F\d+, (F\d+)\]", text)))
    missing = [f for f in used if f"### {f}." not in invtext]
    check("every flag cited is defined in the inventory", not missing, ", ".join(missing))
except FileNotFoundError:
    check("inventory found", False, inv)

w = max(len(r[0]) for r in rows)
for label, verdict, detail in rows:
    print(f"{label:<{w}}  {verdict}" + (f"  ({detail})" if verdict == "FAIL" and detail else ""))
print(f"\n{len(rows) - failed} of {len(rows)} checks pass")
sys.exit(1 if failed else 0)
