# CollecTRI human TF regulons

CollecTRI regulons for the M-b MYC estimator (plan section 7.1) and the FOXO3 activity
score (plan section 7.4).

**Do not edit in place. Re-snapshot and update this file if a newer retrieval is needed.**

## Provenance

- File: `collectri_human_omnipath.tsv.gz` (gzip of the raw TSV as retrieved)
- Upstream: OmniPath web service, https://omnipathdb.org
  CollecTRI: Muller-Dott et al. 2023, Nucleic Acids Research
- Query URL (verbatim):
  `https://omnipathdb.org/interactions?resources=CollecTRI&genesymbols=yes&fields=sources,references,curation_effort&format=tsv&organisms=9606`
- Retrieval date: 2026-08-28
- Method: direct `curl` of the query above; gzipped with no other modification
- SHA-256 of the uncompressed TSV: `d72660703f7ccb8068b75994a1e74986b451d0c098e36d62efd7a88e631d287d`
- Size: 17,192,959 bytes uncompressed / 2,410,971 bytes gzipped

Verified after compression: `gzcat` of the stored file reproduces the SHA-256 above.

## Why a snapshot and not a live decoupleR call

The plan names `decoupleR::get_collectri()`. That call does not work in this environment.

`OmnipathR` 4.0.0 (Bioconductor 3.23) was installed on 2026-08-28. Both
`decoupleR::get_collectri()` and `OmnipathR::collectri()` fail with:

    Caused by error in `full_join()`:
    ! Join columns in `y` must be present in the data.
    x Problem with `ncbi_tax_id`.

This is a package bug, not an environment or connectivity problem: the OmniPath server
returns HTTP 200 in under 0.1 s and the query above returns complete data. The failure
persists after `omnipath_cache_wipe()` and is unaffected by passing `organism = 9606`.

Snapshotting is in any case the better fit for this repo. A live API call at analysis
time is not reproducible; a dated file with a checksum is. This matches the standing
rule that nothing is sourced across repos or services at runtime.

`OmnipathR` and `decoupleR` remain installed and are still used for the scoring step
(`run_ulm` / `run_viper`); only the network retrieval is replaced by this file.

## Contents

64,723 interactions, 1,201 unique TFs, 6,662 unique targets.

| Column | Note |
|---|---|
| `source`, `target` | UniProt accessions. `COMPLEX:` prefix on 22,890 source rows |
| `source_genesymbol`, `target_genesymbol` | HGNC symbols - use these |
| `is_directed`, `is_stimulation`, `is_inhibition` | `True`/`False` strings, not logicals |
| `consensus_direction`, `consensus_stimulation`, `consensus_inhibition` | |
| `sources`, `references`, `curation_effort` | Provenance per interaction |

Regulon sizes for the TFs this project needs (unique targets):

| TF | Targets | Stimulation | Inhibition |
|---|---|---|---|
| `MYC` | 886 | 814 | 169 |
| `FOXO3` | 196 | 166 | 52 |
| `E2F1` | 315 | 285 | 57 |

## Reading notes

- `is_stimulation` / `is_inhibition` arrive as the strings `"True"` / `"False"`. Convert
  explicitly; do not rely on coercion.
- decoupleR's `mor` is built from the sign: +1 where `is_stimulation`, -1 where
  `is_inhibition`. A small number of interactions carry both flags and need an explicit
  rule - decide it in the script and state it, do not let it fall through silently.
- The `COMPLEX:` rows correspond to `split_complexes = FALSE`. `MYC` and `FOXO3` appear
  as plain single-protein sources, so neither is affected.
- Retrieval was filtered to `organisms=9606`. This is a human-only table.
