# Prism .prism archive -> per-panel inventory for figure legends.
# For each data sheet: groups (data sets), rows (row titles / x), per group and
# row the n, mean, s.d. and s.e.m. of the values; for each analysis on it: the
# test, design, post-test method and every adjusted P with the two samples.
# Read-only; the archive is data and nothing in it is executed.
import csv, io, json, math, sys, zipfile

z = zipfile.ZipFile(sys.argv[1])
names = set(z.namelist())
J = lambda n: json.loads(z.read(n))
doc = J("document.json")
attr = doc["sheetAttributesMap"]
want = sys.argv[2:]


def txt(v):
    if isinstance(v, dict):
        return v.get("string", "")
    return v or ""


def dataset(uid):
    p = f"data/sets/{uid}.json"
    return J(p) if p in names else {}


def n_subcols(ds):
    if "replicate ranges" in ds:
        return max(int(r["range"].split("~")[-1]) + 1 for r in ds["replicate ranges"])
    return len(ds.get("replicates", [])) or 1


def num(v):
    try:
        return float(v)
    except ValueError:
        return None


def stats(vals):
    v = [x for x in vals if x is not None]
    n = len(v)
    if n == 0:
        return (0, None, None, None)
    m = sum(v) / n
    sd = math.sqrt(sum((x - m) ** 2 for x in v) / (n - 1)) if n > 1 else None
    return (n, m, sd, sd / math.sqrt(n) if sd is not None else None)


def describe_data(uid):
    s = J(f"data/sheets/{uid}/sheet.json")
    t = s["table"]
    p = f"data/tables/{t['uid']}/data.csv"
    rows = list(csv.reader(io.StringIO(z.read(p).decode("utf-8", "replace")))) if p in names else []
    col = 0
    rtitles = xs = None
    if t.get("rowTitlesDataSet"):
        rtitles = [r[0] if r else "" for r in rows]
        col += 1
    if t.get("xDataSet"):
        xs = [r[col] if len(r) > col else "" for r in rows]
        col += 1
    groups = []
    for i, ds_uid in enumerate(t.get("dataSets", [])):
        ds = dataset(ds_uid)
        k = n_subcols(ds)
        cells = [[num(r[c]) if len(r) > c and r[c].strip() else None for c in range(col, col + k)] for r in rows]
        col += k
        title = txt(ds.get("title")) or f"[{chr(65 + i)}]"
        if t.get("format") == "column":
            flat = [v for row in cells for v in row]
            per = [stats(flat)]
        else:
            per = [stats(row) for row in cells]
        groups.append({"uid": ds_uid, "title": title, "per_row": per})
    return {"title": s["title"], "format": t.get("format"), "dataFormat": t.get("dataFormat"),
            "rows": rtitles, "x": xs, "groups": groups, "nrows": len(rows)}


def comparisons(res):
    out = []
    def walk(o):
        if isinstance(o, dict):
            if ("pValueAdjusted" in o or "qValue" in o) and "sample1" in o:
                s1, s2 = o["sample1"], o["sample2"]
                def lab(s):
                    i = s["id"]
                    c = i.get("col", {}); r = i.get("row", {})
                    return (txt(c.get("title")) or f"[{chr(65 + c.get('idx', 0))}]", r.get("idx"), s["value"].get("n"))
                pv = o.get("pValueAdjusted") or o.get("qValue")
                kind = "adjP" if "pValueAdjusted" in o else "q"
                out.append((lab(s1), lab(s2), pv.get("value") if isinstance(pv, dict) else pv, kind))
            for v in o.values():
                walk(v)
        elif isinstance(o, list):
            for v in o:
                walk(v)
    walk(res)
    return out


def describe_analysis(uid):
    base = f"analyses/{uid}/"
    sh = J(base + "sheet.json")
    par = (J(base + "parameters.json") if base + "parameters.json" in names else {}).get("content", {})
    res = (J(base + "results.json") if base + "results.json" in names else {}).get("content", {})
    method = None
    def find_method(o):
        nonlocal method
        if isinstance(o, dict):
            if "method" in o and isinstance(o["method"], str) and method is None:
                method = o["method"]
            for v in o.values():
                find_method(v)
    find_method(res)
    pt = par.get("postTests") or {}
    pw = pt.get("pairwiseTests") or {}
    design = (pw.get("meansDesignType"), (pw.get("pairwiseDesign") or {}).get("type"))
    models = [k for k in res.keys() if k.endswith("Model")]
    return {"title": sh["title"], "class": sh.get("analysisClass"),
            "inputs": [txt((i or {}).get("title")) for i in (sh.get("inputSheets") or [])],
            "design": par.get("experimentalDesign"), "post": design, "method": method,
            "models": models, "comparisons": comparisons(res),
            "main": [(k, (res.get(models[0], {}).get("anovaTable", {}).get(k, {}).get("fTest", {}).get("pValue", {}) or {}).get("value"))
                     for k in ["rowFactor", "colFactor", "rowFactorVsColFactor"]] if models else []}


analyses = [describe_analysis(u) for u in doc["sheets"]["analyses"]]
by_input = {}
for a in analyses:
    for t in a["inputs"]:
        by_input.setdefault(t, []).append(a)

fmt = lambda v, d=3: "-" if v is None else f"{v:.{d}g}"
for uid in doc["sheets"]["data"]:
    d = describe_data(uid)
    if want and not any(d["title"].startswith(w) for w in want):
        continue
    print(f"\n### {d['title']}   [{d['format']}/{d['dataFormat']}]")
    labels = None
    if d["format"] != "column":
        labels = [(d["rows"][i] if d["rows"] and d["rows"][i] else "") + (f" x={d['x'][i]}" if d["x"] and d["x"][i] else "") for i in range(d["nrows"])]
    for g in d["groups"]:
        if d["format"] == "column":
            n, m, sd, se = g["per_row"][0]
            if n:
                print(f"   {g['title']:<22} n={n} mean={fmt(m)} sd={fmt(sd)} sem={fmt(se)}")
        else:
            ns = [p[0] for p in g["per_row"]]
            if any(ns):
                nn = sorted(set(x for x in ns if x))
                print(f"   {g['title']:<22} n per row {nn} over {sum(1 for x in ns if x)} rows")
                for i, (n, m, sd, se) in enumerate(g["per_row"]):
                    if n:
                        print(f"       {labels[i]:<18} n={n} mean={fmt(m)} sd={fmt(sd)} sem={fmt(se)}")
    for a in by_input.get(d["title"], []):
        print(f"   >> {a['title']} | {a['class']} | design={a['design']} | post={a['post']} | method={a['method']}")
        if a["main"]:
            print("      main effects:", ", ".join(f"{k} P={fmt(v)}" for k, v in a["main"] if v is not None))
        for s1, s2, p, st in a["comparisons"]:
            r1 = labels[s1[1]] if labels and s1[1] is not None and s1[1] < len(labels) else s1[1]
            r2 = labels[s2[1]] if labels and s2[1] is not None and s2[1] < len(labels) else s2[1]
            print(f"      {s1[0]}[{r1}] (n={s1[2]}) vs {s2[0]}[{r2}] (n={s2[2]}): {st}={fmt(p, 4)}")
