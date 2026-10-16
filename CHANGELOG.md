# Changelog

## 0.3.0 — 2026-10-16

### Added
- `from_log` (CLI `--from-log`): aggregate log-normalised input,
  `log1p(count / total * scale_factor)`. Raw counts are recovered as
  `expm1(x) * total / scale_factor` from each cell's stored total
  (`total_col`; `scale_factor` inferred unless given), then summed. Without
  stored totals, scagg prints a message and falls back to the
  smallest-non-zero trick, extended to cells with no count of exactly 1. A
  warning is raised if the recovered values are not close to whole numbers.
  In R this reads the `data` layer.
- `mean_normalized` (CLI `--mean-normalized`), only with `from_log`: return
  `log1p(mean(expm1(x)))` per metacell instead of summed counts.
- `uns["scagg_values"]` (Python) records whether `X` holds counts, recovered
  counts or mean normalised values. Python warns when the matrix being summed
  has non-integer values and `from_log` is off.

### Fixed
- R: aggregation no longer goes through
  `Seurat::AggregateExpression(return.seurat = TRUE)`, which copied the whole
  object with `subset()` and ran `ScaleData()` on every gene (more than 32 GB
  on large objects). Counts are now summed directly with a sparse indicator
  product, including across split Seurat v5 `counts.*` layers. Summed counts
  are unchanged; on a 100k-cell × 20k-gene test object, peak memory fell from
  10.3 GB to 5.9 GB (most of that is the input object) and run time from
  22 s to 5 s.
- Metacells could silently pool cells from different samples or cell types
  when group values containing `_` joined to the same ID (e.g. `("A_1", "x")`
  and `("A", "1_x")` both give `A_1_x`). All functions now check that every
  metacell holds cells from a single `group_vars` combination and stop with
  an error otherwise.
- Python: summing per metacell is one sparse product instead of a loop over
  metacells (which was quadratic), and no longer copies the input. Output is
  unchanged.

### Changed
- R: the result has a log-normalised `data` layer but no `scale.data`; run
  `ScaleData()` on the genes you need. New `normalize` argument
  (default `TRUE`); `normalize = FALSE` returns counts only.
- R: column names of the result are now the original metacell IDs
  (`mc1_A_s1`) instead of Seurat's renamed form (`mc1-A-s1`), so they match
  `metacell_id` and the membership table.

## 0.2.0 — 2026-10-15

### Added
- `method="similarity"` for `make_metacells_by_size` and `make_metacells_by_num`
  (Python, R and CLI). Cells that are close in an embedding (`use_rep` in
  Python, `reduction` in R) are pooled together via balanced recursive PCA
  bisection with size-constrained 2-means at each split, followed by optional
  size-preserving pairwise refinement (`refine`). Metacell counts and sizes are
  identical to the random method (balanced to within one cell); the
  partition is deterministic and identical between Python and R.
- `n_dims` to restrict the embedding to its first N dimensions.
- `n_init_pcs` (CLI `--n-init-pcs`): try each split starting from PC1 … PCN
  and keep the cut with the lowest within-half sum of squares. Helps when a
  broad gradient dominates PC1 but sub-states separate along a later PC.
  Default `1` (PC1 only).
- `compactness` column in metacell metadata for similarity metacells, and a
  log line comparing it against size-matched random metacells.
- CLI flags `--method`, `--use-rep`, `--n-dims`, `--no-refine`, `--n-init-pcs`.
- R testthat suite and an R job in CI.

### Fixed
- R: with Seurat v5, every metadata column of the aggregated object
  (`n_cells`, group and meta variables) was `NA`, because
  `AggregateExpression` renames groups from `mc1_A_s1` to `mc1-A-s1` and the
  metadata join used the original IDs. Columns are now joined correctly, and
  `metacell_id` keeps the original underscore form that matches the
  membership table.

### Changed
- R package version aligned with Python (was 0.1.0 vs 0.1.2).
- The default `method="random"` output is unchanged from 0.1.x for a given
  seed.

## 0.1.2
- Docker image and GHCR publish workflow.

## 0.1.1
- `idents_col` parameter for R functions; PyPI publish workflow.

## 0.1.0
- Initial release.
