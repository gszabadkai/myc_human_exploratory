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

## 4. SAVE THE FULL COEFFICIENT MATRIX, THE FULL VIF VECTOR AND THE FITTED
##    MODEL OBJECTS. Never a selected row.

# Added 2026-10-08, from the third block of this class. Declaration section 14.7
# fixes it for the respiratory-axis line; it applies to the arm.
#
# WHAT HAPPENED. `E40` and `E41` carried `PROLIF_DISJOINT` in every rung from m1
# onward and extracted coefficients through a helper taking `co["OX", ]`. Every
# other term's estimate, SE, interval, p and VIF was computed and thrown away.
# The saved tables are one row per cohort-by-rung, so the `terms` column records
# the formula as a string while the coefficients it names exist nowhere. Neither
# object saved its fitted models either. `E43` had to refit all eight models to
# recover one term that had already been estimated twice.
#
# THE RULE. Every script that fits anything saves:
#
#   - the FULL coefficient matrix, every term, as a tidy table - set, term,
#     estimate, SE, statistic, p, and the interval;
#   - the FULL VIF vector, every term, with its Df and the measure named;
#   - the FITTED MODEL OBJECTS themselves.
#
# WHY A SELECTED ROW IS NEVER ENOUGH. The coefficient you did not keep is the
# one the next question needs. That has now happened three times: script 09's
# discarded Block C main effects, script 13's F1 main effect printed once in an
# `if (FALSE)` sandbox, and `E40`/`E41`'s PROLIF row. Each cost a refit and one
# of them cost a declaration amendment.
#
# ONE CAVEAT FOR VERIFICATION, found by `E43`'s dry runs. A FITTED MODEL IS NOT
# DIGEST-STABLE ACROSS RUNS. `coxph` and `glm` objects carry `terms`, `formula`,
# `model` and `family`, each holding an ENVIRONMENT whose address differs
# between sessions, so `identical()` on two runs' fits returns FALSE even when
# every number agrees. Verified on `E43`: coefficients and variance matrices
# identical in every fit, only those four carriers differing.
#
#   SO: save a `fit_digests` table - one row per fit, with digests of `coef()`
#   and `vcov()` - and EXCLUDE the fitted objects from any digest comparison.
#   Say so in the script and in the commit, or a verification pass reads a
#   cosmetic difference as a failure.
#
# AND READ A MULTI-DF GVIF ON THE RIGHT SCALE. A generalised VIF for a term with
# Df > 1 is not comparable to a 1-df VIF threshold. Compare `GVIF^(1/(2*Df))`
# against the SQUARE ROOT of the 1-df threshold (Fox and Monette). `E43`'s first
# version compared a raw multi-df GVIF of 7.979 against a threshold of 5 and
# flagged a caveat that does not exist: adjusted, that term is 1.414 on 3 df.
