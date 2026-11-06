# Extract every drawn two-point segment from a Prism/Illustrator PDF, in page
# points, so a panel's error bars can be read against its own axis. Read-only.
import re, sys, zlib

raw = open(sys.argv[1], "rb").read()
blocks = []
for m in re.finditer(rb"stream\r?\n", raw):
    s = m.end()
    e = raw.find(b"endstream", s)
    try:
        t = zlib.decompress(raw[s:e]).decode("latin-1")
    except Exception:
        continue
    if " cm\n" not in t and " cm " not in t:
        continue
    for b in re.finditer(r"q 1 0 0 1 (-?[\d.]+) (-?[\d.]+) cm\s*(.*?)\sQ", t, re.S):
        X, Y, body = float(b.group(1)), float(b.group(2)), b.group(3)
        pts = re.findall(r"(-?[\d.]+) (-?[\d.]+) (?:m|l)\b", body)
        if len(pts) == 2 and re.search(r"\bS\b", body):
            (x1, y1), (x2, y2) = [(float(a) + X, float(c) + Y) for a, c in pts]
            blocks.append((x1, y1, x2, y2))

x0, x1b, y0, y1b = map(float, sys.argv[2:6])
sel = [s for s in blocks if x0 <= s[0] <= x1b and y0 <= s[1] <= y1b]
vert = sorted([s for s in sel if abs(s[0] - s[2]) < 0.01 and abs(s[3] - s[1]) > 0.3], key=lambda s: (s[0], s[1]))
horiz = sorted([s for s in sel if abs(s[1] - s[3]) < 0.01], key=lambda s: s[1])
print(f"{len(blocks)} segments in file; {len(sel)} in window; {len(vert)} vertical, {len(horiz)} horizontal")
print("VERTICAL (x, ylow, yhigh, length pt):")
for s in vert:
    lo, hi = sorted((s[1], s[3]))
    print(f"   {s[0]:8.2f}  {lo:8.2f} {hi:8.2f}  {hi - lo:7.3f}")
print("HORIZONTAL (y, xlow, xhigh, length pt)  - ticks and axis:")
for s in horiz[:25]:
    lo, hi = sorted((s[0], s[2]))
    print(f"   {s[1]:8.2f}  {lo:8.2f} {hi:8.2f}  {hi - lo:7.3f}")
