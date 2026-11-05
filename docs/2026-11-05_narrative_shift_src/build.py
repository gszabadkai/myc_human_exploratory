import json, math, sys
D = sys.argv[1]
a = json.load(open(f"{D}/data.json")); c = json.load(open(f"{D}/data_c.json"))
raw = {(r["cohort"], r["readout"], r["exposure"]): r["rho"] for r in a["raw"]}
RAW = [("OXPHOS level", "OX_gsva", None), ("MYC activity, reference score", "MYC_M_a", "the “not MYC” reading"),
       ("MYC activity, MYC_UP.V1_UP", "MYC_clean", None), ("Proliferation", "PROLIF", None)]
raw_rows = [{"row": lab, "note": note, "TCGA": raw[("TCGA", "S6", k)], "SCAN-B": raw[("SCAN-B", "S6", k)]} for lab, k, note in RAW]
sig = {(r["cohort"], r["estimator"]): r for r in a["sig"]}
g3 = {(r["cohort"], r["estimator"]): r["rho_MYC_PROLIF"] for r in a["g3"]}
LAD = [("MYC_UP.V1_UP", "clean", "the score in Fig. 4g"), ("CollecTRI MYC regulon", "regulon", ""),
       ("Felsher", "M_a", "the reference score"), ("Hallmark MYC targets", "entangled", "most tied to proliferation")]
ladder = []
for lab, k, note in LAD:
    rows = {}
    for coh in ("TCGA", "SCAN-B"):
        r = sig[(coh, k)]
        rows[coh] = {"m1": [r["b_MYC_m1"], r["myc_lo_m1"], r["myc_hi_m1"]],
                     "m3": [r["b_MYC_m3"], r["myc_lo_m3"], r["myc_hi_m3"]], "void": bool(r["void"])}
    ladder.append({"score": lab, "note": note, "rho": {coh: g3[(coh, k)] for coh in ("TCGA", "SCAN-B")}, "rows": rows})
# Order check: by mean correlation with proliferation, ascending.
means = [sum(l["rho"].values()) / 2 for l in ladder]
assert means == sorted(means), means
ci = {(r["readout"], r["cohort"]): r for r in c["ci"]}
cells = []
for x in c["cells"]:
    r = ci[(x["readout"], x["cohort"])]
    lo, hi = {"OX": (r["ox_lo_m3"], r["ox_hi_m3"]), "MYC": (r["myc_lo_m3"], r["myc_hi_m3"]), "PROLIF": (r["pr_lo"], r["pr_hi"])}[x["term"]]
    cells.append({**{k: x[k] for k in ("readout", "cohort", "term", "value", "fill", "ink")}, "lo": lo, "hi": hi})
stops = [l.strip() for l in open(f"{D}/stops.txt") if l.strip()]
heat = {"rows": [{"key": "BCL2L1", "gene": "BCL2L1", "protein": "BCL-xL", "label": "BCL2L1"},
                 {"key": "BBC3", "gene": "BBC3", "protein": "PUMA", "label": "BBC3"},
                 {"key": "S6", "label": "Signature"}], "cells": cells, "stops": stops}
p = a["pcr"]; ox = {r["rung"]: r for r in p["ox"]}; pr = p["prolif"][0]; dl = p["delta"][0]
E = math.exp
pcr = [
 {"label": "OXPHOS", "sub": "without proliferation term", "kind": "ghost", "est": E(ox["m0"]["fe_estimate"]),
  "lo": E(ox["m0"]["fe_ci_lo"]), "hi": E(ox["m0"]["fe_ci_hi"]), "detail": "fixed effect; cohort signs + − +"},
 {"label": "OXPHOS", "sub": "with proliferation term", "kind": "ink", "est": E(ox["m1"]["fe_estimate"]),
  "lo": E(ox["m1"]["fe_ci_lo"]), "hi": E(ox["m1"]["fe_ci_hi"]), "detail": "fixed effect; signs − − −; I² 0%"},
 {"label": "Proliferation", "sub": "same model", "kind": "ink", "est": E(pr["re_estimate"]),
  "lo": E(pr["re_ci_lo"]), "hi": E(pr["re_ci_hi"]), "detail": "random effects; signs + + +; I² %d%%" % round(pr["I2"])},
 {"label": "Proliferation ÷ OXPHOS", "sub": "difference, tested", "kind": "ink", "shape": "diamond", "divider": True,
  "est": E(dl["fe_estimate"]), "lo": E(dl["fe_ci_lo"]), "hi": E(dl["fe_ci_hi"]), "detail": "ratio of odds ratios; fixed effect; I² 0%"},
]
data = {"raw": raw_rows, "ladder": ladder, "heat": heat, "pcr": pcr}
html = open(f"{D}/template.html", encoding="utf-8").read()
assert html.count("__DATA__") == 1
open(f"{D}/fig4_narrative_shift.html", "w", encoding="utf-8").write(html.replace("__DATA__", json.dumps(data, ensure_ascii=False)))
for r in raw_rows: print("a", r["row"], round(r["TCGA"], 3), round(r["SCAN-B"], 3))
for l in ladder: print("b", l["score"], {k: (round(v["m1"][0], 3), round(v["m3"][0], 3), v["void"]) for k, v in l["rows"].items()}, {k: round(v, 2) for k, v in l["rho"].items()})
for r in pcr: print("d", r["label"], r["sub"], "%.2f [%.2f, %.2f]" % (r["est"], r["lo"], r["hi"]))
print("written", len(html))
