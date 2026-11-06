# Map every P value printed on an assembled figure PDF to the panel that carries
# it, by position on the page, so a legend's claim about a panel is checked
# against that panel and not against the reading order of the file. Read-only.
#
#   python3 docs/legends/tools/pdf_pvalue_map.py NM_Fig2.pdf
#
# Panel letters are found as the largest single lower-case letters on the page;
# where two panels share a row the letters' x positions give the column, which
# the printout carries so a borderline value can be placed by hand. Needs
# pdftotext (poppler).
import html
import re
import statistics
import subprocess
import sys

pdf = sys.argv[1]
xml = subprocess.run(["pdftotext", "-bbox-layout", pdf, "-"],
                     capture_output=True, text=True).stdout
words = []
for m in re.finditer(r'<word xMin="([\d.]+)" yMin="([\d.]+)" xMax="([\d.]+)" '
                     r'yMax="([\d.]+)">(.*?)</word>', xml):
    words.append((float(m[1]), float(m[2]), float(m[3]), float(m[4]),
                  html.unescape(m[5])))

# A panel letter is a single lower-case letter set far larger than the panels'
# own text: on these figures 17.3 pt against a median word height of 8 to 10.
med = statistics.median(w[3] - w[1] for w in words) if words else 0
letters = sorted((w for w in words
                  if re.fullmatch(r"[a-z]", w[4]) and (w[3] - w[1]) >= 1.5 * med),
                 key=lambda w: (w[1], w[0]))
print("panel letters: " + ", ".join(f"{w[4]} (x={w[0]:.0f}, y={w[1]:.0f})" for w in letters))

vals = []
for i, w in enumerate(words):
    if (w[4] in ("P", "q") and i + 2 < len(words) and words[i + 1][4] in "=<"
            and re.fullmatch(r"[\d.]+", words[i + 2][4])):
        vals.append((w[0], w[1], f"{w[4]} {words[i+1][4]} {words[i+2][4]}"))


def panel_of(px, py):
    above = [L for L in letters if L[1] <= py + 2]
    if not above:
        return "?"
    row = max(L[1] for L in above)
    same_row = [L for L in above if L[1] >= row - 1]
    left = [L for L in same_row if L[0] <= px]
    return (max(left, key=lambda L: L[0]) if left else min(same_row, key=lambda L: L[0]))[4]


print(f"\n{len(vals)} printed values:")
for px, py, s in sorted(vals, key=lambda v: (v[1], v[0])):
    print(f"  {panel_of(px, py)}   x={px:7.1f} y={py:7.1f}   {s}")
