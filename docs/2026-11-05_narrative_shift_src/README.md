# Source of the narrative-shift explainer

**What it is.** A one-page explainer for co-authors on how the human arm of Fig. 4
moved:
- **from** "the signature tracks OXPHOS but not MYC";
- **to** OXPHOS against proliferation, with MYC activity going with the signature
  in the same direction as OXPHOS.

It is EXPLORATORY and POST HOC, and it is **not a manuscript figure**.

**Where it lives.** It is published as a private Claude artifact,
"Fig. 4 Narrative Shift" (<https://claude.ai/artifact/PdASXUXr2x6nu91oXwVvXY>,
2026-11-05). Co-authors can open it only once the author shares it from the
page's Share menu.

**It computes nothing new.** Every number is read from a saved results object,
each verified against digests committed before the author's run:

| panel | what | read from |
|---|---|---|
| a | the signature's raw Spearman correlations | E49 `raw_correlations` |
| b | MYC activity's coefficient on the signature, without and with proliferation, on four MYC scores | E49 `signature_coefficients`; correlations with proliferation from E47's gate G3 |
| c | Fig. 4g's cells and colours for `BCL2L1`, `BBC3` and the signature | E50 `values`, checked equal to E49's coefficients (18 of 18); intervals from E49 `all_readouts_m3` |
| d | pooled pCR odds ratios | E40, E43 and E44 `meta` |

**The files:**
- **`fig4_narrative_shift.html`** is the page as published. It is a fragment: the
  artifact platform wraps it in a document skeleton.
- **`template.html`** is the page with a `__DATA__` placeholder.
- **`extract.R` and `extract_c.R`** read the saved objects (read-only) and write
  `data.json` and `data_c.json`.
- **`stops.txt`** holds the colour key's stops: Fig. 1g's CIELAB scale
  (`scales::pal_gradient_n`, as E50 drew it).
- **`build.py`** injects the data into the template and writes the page.

**To rebuild**, from the repo root, with `results/` present (it is gitignored):

```
Rscript docs/2026-11-05_narrative_shift_src/extract.R docs/2026-11-05_narrative_shift_src/data.json
Rscript docs/2026-11-05_narrative_shift_src/extract_c.R docs/2026-11-05_narrative_shift_src/data_c.json
Rscript -e 'f <- scales::pal_gradient_n(c("#4A3525","#FAFAFA","#2A8A6D"), values = c(0, 0.5, 1)); cat(f(seq(0, 1, length.out = 13)), sep = "\n")' > docs/2026-11-05_narrative_shift_src/stops.txt
python3 -I -B docs/2026-11-05_narrative_shift_src/build.py docs/2026-11-05_narrative_shift_src
```

From this folder's committed files, `build.py` reproduces the published page
byte for byte (checked on 2026-11-05).

**Two notes:**
- **OXPHOS's adjusted pCR interval is [0.75, 0.96].** The lower bound is
  exp(-0.2813) = 0.7548; the 2 November brief rounded it to 0.76.
- **The page was not previewed before publishing.** The local headless browser
  could not run in the session's sandbox, so only its script was syntax-checked.
  The author reviews it live.
