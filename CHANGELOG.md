# Changelog

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
