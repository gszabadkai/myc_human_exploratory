# R Coding Instructions
# ====================
# General rules for writing R scripts in this project.
# These address known pitfalls and namespace conflicts.

## 1. Tibble printing: avoid `print(n = ...)` after `head()`

# When a tibble passes through `head()`, the result may become a plain
# data.frame. Passing `n` to `print.data.frame()` is interpreted as
# `na.print`, causing errors.
#
# BAD:
#   df %>% head(20) %>% print(n = 20)
#
# GOOD (head already limits rows):
#   df %>% head(20) %>% print()
#
# ALSO GOOD (if you know it's a tibble, use n directly without head):
#   df %>% print(n = 20)

## 2. Namespace conflicts: always use `dplyr::count()`

# `count()` can conflict with other packages (e.g. plyr). Always
# use the explicit namespace to avoid silent errors:
#
# BAD:
#   df %>% count(group)
#
# GOOD:
#   df %>% dplyr::count(group)
#
# This also applies to other commonly conflicting functions:
#   dplyr::select()   — conflicts with MASS::select
#   dplyr::filter()   — conflicts with stats::filter
#   dplyr::lag()      — conflicts with stats::lag
#   dplyr::rename()   — conflicts with plyr::rename
#
# In this project, dplyr::select and dplyr::filter are already used
# with explicit namespaces where conflicts are likely.

## 3. SAVE THE ANALYSIS FRAME. A numbered script's `.rds` must carry the
##    per-patient frame it fitted on, not only its summary tables.

# Added 2026-10-04, from a real block. `E41` saved 26 summary tables and no
# per-patient frame. `X01` then needed to stratify METABRIC by PAM50 and could
# not: nothing on disk held the scores, survival, covariates and identifiers
# together. The only per-patient element, `er_discordant`, was 127 rows with no
# scores and no survival. Commit `51465b6` had to add `frame = DM` afterwards
# and the author had to re-source `E41`.
#
# THE RULE. Every numbered script that fits anything saves, alongside its
# results, the data frame the fits were run on - one row per analysis unit,
# carrying at minimum:
#
#   - the join key (`sample`, `sample_id`, `patient`, whatever the arm uses);
#   - every exposure and covariate that entered any model, on the scale it
#     entered on, AND the raw version if one was transformed;
#   - every outcome, with its time variable where there is one;
#   - every stratifier the script used, and every candidate stratifier it
#     audited and rejected.
#
# BOUNDARY. This does NOT mean save the expression matrix. The frame is the
# analysis-unit table - 1,981 rows for `E41`, trivial on disk - not the
# 23,043 x 1,981 input. Big inputs stay where their provenance README says.
#
# AND SAVE MORE COLUMNS THAN THE SCRIPT FITTED. `X01` ran into this twice. The
# second time was worse: asked to weigh confounding by indication behind a
# treatment-group reversal, it could not, because `grade`, `size`,
# `lymph_nodes_positive` and `age_at_diagnosis` sit in METABRIC's clinical
# table and `E41` did not carry them into its frame. They cost nothing to
# carry and their absence left a question unanswerable from disk. So the frame
# holds the covariates a downstream question might plausibly need, not only
# the ones the declaring script put in a model.
#
# WHY THIS IS NOT A LICENCE TO RE-FIT. The frame exists so a downstream script
# can ask a NEW question of the same data without rescoring. It does NOT
# license refitting what the declaring script already fitted, or re-reading its
# declared quantities under a different rule. The declaration that governed the
# original fit still governs it.
#
# ADDING A FRAME TO AN ALREADY-VERIFIED OBJECT IS SAFE AND MUST BE SHOWN TO BE.
# It appends an element and changes none of the existing ones, so recorded
# per-element digests still hold. `51465b6` demonstrated this rather than
# asserting it: all 26 of `E41`'s pre-existing elements came back
# byte-identical and the author's re-run matched them, so a digest check finds
# 26 matching plus one new. Any such amendment states that explicitly and the
# re-source requirement with it.
