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
