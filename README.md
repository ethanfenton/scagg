# scagg

**Pseudobulk and metacell aggregation for single-cell RNA-seq**

`scagg` provides three aggregation strategies for turning single-cell data into
sample-level or metacell-level objects, ready for bulk-style differential
expression (DESeq2, edgeR, dream) or other downstream analyses:

| Strategy | Function | When to use |
|---|---|---|
| **Pseudobulk** | `make_pseudobulk` | One aggregate per unique sample × cell-type combination |
| **Metacells by size** | `make_metacells_by_size` | ~N cells per metacell, maximises usage of all cells |
| **Metacells by count** | `make_metacells_by_num` | Exactly N metacells per group, balances replication |

Both metacell strategies can group cells **randomly** (default) or by
**similarity** in an embedding (`method="similarity"`): cells that are close in
PCA / Harmony / scVI space are pooled together, while metacell sizes stay
balanced to within one cell.

Available for **Python (AnnData/scanpy)** and **R (Seurat)**.
Full cell → metacell traceability is built in.

[![CI](https://github.com/ethanfenton/scagg/actions/workflows/ci.yml/badge.svg)](https://github.com/ethanfenton/scagg/actions)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)

---

## Installation

### Python

```bash
pip install scagg

# or the development version
pip install git+https://github.com/ethanfenton/scagg.git
```

Requirements: Python ≥ 3.10, anndata, numpy, pandas, scipy (all present in
any standard scanpy environment).

### R

```r
if (!requireNamespace("remotes")) install.packages("remotes")
remotes::install_github("ethanfenton/scagg", subdir = "R/scagg")
```

Requirements: Seurat, dplyr.

### Docker

No Python environment needed — mount your data directory and run:

```bash
docker pull ghcr.io/ethanfenton/scagg:latest

docker run --rm \
  -v /path/to/your/data:/data \
  ghcr.io/ethanfenton/scagg:latest \
  pseudobulk \
  --input  /data/input.h5ad \
  --output /data/output.h5ad \
  --group-vars cell_type sample

# All sub-commands work the same way
docker run --rm -v $(pwd):/data ghcr.io/ethanfenton/scagg:latest --help
```

The `/data` mount is where input and output files must live — paths inside the
container must start with `/data/`.  Versioned images are available as
`ghcr.io/ethanfenton/scagg:0.2.0`, etc.

---

## Example

A complete run that builds similarity metacells from public PBMC data.

### Python

```python
import scanpy as sc
import scagg

# 2,638 PBMCs: raw counts + published cell-type labels
adata = sc.datasets.pbmc3k()
labels = sc.datasets.pbmc3k_processed().obs["louvain"]
adata = adata[labels.index].copy()
adata.obs["cell_type"] = labels
adata.layers["counts"] = adata.X.copy()

# Standard preprocessing to get a PCA embedding
sc.pp.normalize_total(adata)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, n_top_genes=2000)
sc.pp.pca(adata, n_comps=30)

# ~20-cell metacells of similar cells within each cell type
mc = scagg.make_metacells_by_size(
    adata,
    group_vars = "cell_type",
    cell_size  = 20,
    layer      = "counts",      # sum raw counts
    method     = "similarity",
    use_rep    = "X_pca",
)
# [scagg] by_size_20_similarity: 2638 cells → 129 metacells (0 excluded)

mc.obs[["cell_type", "n_cells", "compactness"]].head()
mc.uns["metacell_membership"]   # which cell went into which metacell
```

### R

```r
library(Seurat)
library(scagg)

# 80 PBMCs bundled with SeuratObject, already clustered, with a PCA
pbmc <- UpdateSeuratObject(pbmc_small)

res <- make_metacells_by_size(
  pbmc,
  group_vars = "RNA_snn_res.1",
  cell_size  = 5,
  assays     = "RNA",
  method     = "similarity",
  reduction  = "pca"
)
# [scagg] 80 cells → 15 metacells (0 excluded)

res$obj                                    # Seurat object, one column per metacell
head(res$obj@meta.data[, c("metacell_id", "n_cells", "compactness")])
head(res$membership)                       # which cell went into which metacell
```

In real data, add your sample column to `group_vars`
(e.g. `c("cell_type", "sample")`) so that metacells never mix samples.
Cells are only pooled with cells that share every `group_vars` value; this
holds for both `method = "random"` and `method = "similarity"`, and scagg
stops with an error if two groups would ever map to the same metacell ID.

---

## Quick start — Python

```python
import scanpy as sc
import scagg

adata = sc.read_h5ad("my_data.h5ad")

# True pseudobulk: one obs per sample × cell_type
pb = scagg.make_pseudobulk(
    adata,
    group_vars = ["cell_type", "sample"],
    save_membership = "results/pseudobulk_membership.csv",
)

# Metacells of ~10 cells each
mc_size = scagg.make_metacells_by_size(
    adata,
    group_vars = ["cell_type", "sample"],
    cell_size  = 10,
    save_membership = "results/metacell_size10_membership.csv",
)

# Exactly 5 metacells per cell_type × sample group
mc_num = scagg.make_metacells_by_num(
    adata,
    group_vars  = ["cell_type", "sample"],
    n_metacells = 5,
    save_membership = "results/metacell_num5_membership.csv",
)
```

### Similarity-based metacells

```python
sc.pp.pca(adata)   # or use an existing obsm key, e.g. "X_harmony", "X_scVI"

mc_sim = scagg.make_metacells_by_size(
    adata,
    group_vars = ["cell_type", "sample"],
    cell_size  = 10,
    method     = "similarity",
    use_rep    = "X_pca",
    n_dims     = 30,          # optional: first 30 dims only
)
mc_sim.obs["compactness"]    # mean cell-to-centroid distance per metacell
# logs: [scagg] median metacell compactness: <x> (random metacells of the same sizes: <y>)
```

### Log-normalised input

Sums are only meaningful on counts: adding log values multiplies the
underlying expression. If you only have log-normalised values,
`log1p(count / total * scale_factor)` (e.g. `adata.X` after `sc.pp.log1p`, or
Seurat's `data` layer), pass `from_log=True`. scagg then works in linear space:

```python
# Recover raw counts, then sum them (usable for DESeq2 / edgeR)
mc = scagg.make_metacells_by_size(adata, ["cell_type", "sample"], 10,
                                  from_log=True)

# Or: average normalised profile per metacell, log1p(mean(expm1(x)))
mc_avg = scagg.make_metacells_by_size(adata, ["cell_type", "sample"], 10,
                                      from_log=True, mean_normalized=True)
```

- **Recovered counts** (`from_log=True`, default `mean_normalized=False`):
  `count = expm1(x) * total / scale_factor`, using each cell's total count
  from `total_col` (default: `total_counts`, `n_counts` or `nCount_RNA` in
  Python; `nCount_<assay>` in R). `scale_factor=None` infers the target sum,
  so Seurat's 1e4 and scanpy's default median target both work. If no totals
  are stored, scagg prints a message and falls back to the
  smallest-non-zero trick: every value in a cell is a whole multiple of
  one count, so the cell is rescaled until all its values are whole numbers.
  Either way, scagg warns if the recovered values are not close to whole
  numbers (e.g. the data were scaled or batch-corrected).
- **Mean normalised profile** (`mean_normalized=True`, only with
  `from_log=True`): every cell weighs the same regardless of depth. Good for
  plots and correlations; not counts, so not for count-based DE. In R it is
  stored in the `data` layer, with no `counts` layer.

`uns["scagg_values"]` records which one you got (`"counts"`,
`"recovered_counts"` or `"mean_normalized"`). Python warns when the matrix it
is about to sum has non-integer values and `from_log` is off.

In R, `from_log = TRUE` reads the `data` layer instead of `counts`. It is not
needed for SCT assays, whose `counts` layer already holds corrected counts.

### Access results

```python
# Result is a standard AnnData
mc_size.obs          # metacell metadata (n_cells, cell_type, sample, …)
mc_size.X            # summed counts matrix (metacells × genes)
mc_size.uns["metacell_membership"]  # DataFrame: cell_barcode → metacell_id
```

### CLI

```bash
# Pseudobulk
scagg pseudobulk \
  --input data.h5ad --output results/pb.h5ad \
  --group-vars cell_type sample \
  --save-membership results/pb_membership.csv

# Metacells by size
scagg metacells-by-size \
  --input data.h5ad --output results/mc.h5ad \
  --group-vars cell_type sample \
  --cell-size 10 \
  --save-membership results/mc_membership.csv

# Metacells by count
scagg metacells-by-num \
  --input data.h5ad --output results/mc.h5ad \
  --group-vars cell_type sample \
  --n-metacells 5

# Similarity-based metacells (either metacell sub-command)
scagg metacells-by-size \
  --input data.h5ad --output results/mc_sim.h5ad \
  --group-vars cell_type sample \
  --cell-size 10 \
  --method similarity --use-rep X_pca --n-dims 30 --n-init-pcs 3
```

---

## Quick start — R (Seurat)

```r
library(scagg)

# True pseudobulk
res <- make_pseudobulk(
  so_obj     = seurat_obj,
  group_vars = c("cell_type", "sample"),
  assays     = c("SCT", "RNA"),
  save_membership = "results/pseudobulk_membership.csv"
)
pb_so      <- res$obj         # aggregated Seurat object
membership <- res$membership  # data.frame: barcode → metacell_id

# Metacells of ~10 cells
res <- make_metacells_by_size(
  so_obj     = seurat_obj,
  group_vars = c("cell_type", "sample"),
  cell_size  = 10,
  save_membership = "results/mc_size10_membership.csv"
)

# Exactly 5 metacells per group
res <- make_metacells_by_num(
  so_obj      = seurat_obj,
  group_vars  = c("cell_type", "sample"),
  n_metacells = 5
)

# Similarity-based metacells (uses the "pca" reduction by default)
res <- make_metacells_by_size(
  so_obj     = seurat_obj,
  group_vars = c("cell_type", "sample"),
  cell_size  = 10,
  method     = "similarity",
  reduction  = "harmony",
  n_dims     = 30
)
```

### Carry custom metadata

```r
res <- make_metacells_by_size(
  so_obj     = seurat_obj,
  group_vars = c("cell_type", "sample"),
  cell_size  = 10,
  meta_vars  = c("Treatment", "Timepoint", "sex", "Prepper", "batch")
)
```

---

## API reference

### Python

```
scagg.make_pseudobulk(adata, group_vars, *, meta_vars, layer, save_membership,
                      from_log, mean_normalized, total_col, scale_factor)
scagg.make_metacells_by_size(adata, group_vars, cell_size, *, cell_min, meta_vars,
                              layer, seed, save_membership,
                              method, use_rep, n_dims, refine, n_init_pcs,
                              from_log, mean_normalized, total_col, scale_factor)
scagg.make_metacells_by_num(adata, group_vars, n_metacells, *, min_cells, meta_vars,
                             layer, seed, save_membership,
                             method, use_rep, n_dims, refine, n_init_pcs,
                             from_log, mean_normalized, total_col, scale_factor)
```

| Parameter | Default | Notes |
|---|---|---|
| `group_vars` | required | One or more obs column names |
| `cell_size` | required | Target cells per metacell |
| `cell_min` | 70 % of cell_size | Min cells to include a group |
| `n_metacells` | required | Metacells per group |
| `min_cells` | n_metacells | Min cells to include a group |
| `meta_vars` | all obs columns | Which metadata to carry forward |
| `layer` | None (use X) | AnnData layer to aggregate |
| `seed` | 42 | RNG seed for reproducibility (random method only) |
| `save_membership` | None | CSV path for traceability output |
| `method` | `"random"` | `"random"` or `"similarity"` |
| `use_rep` | `"X_pca"` | obsm embedding for `method="similarity"` (R: `reduction = "pca"`) |
| `n_dims` | None (all) | Use only the first N embedding dimensions |
| `refine` | True | Size-preserving refinement after bisection |
| `n_init_pcs` | 1 | Try each split from PC1…PCN, keep the tightest (see below) |
| `from_log` | False | Input is log-normalised; recover raw counts before summing |
| `mean_normalized` | False | With `from_log`: return `log1p(mean(expm1(x)))` instead |
| `total_col` | auto | Column with per-cell totals for `from_log` |
| `scale_factor` | None (infer) | Normalisation target sum for `from_log` |

### R

```
make_pseudobulk(so_obj, group_vars, meta_vars, assays, save_membership,
                idents_col = NULL, normalize = TRUE, from_log = FALSE,
                mean_normalized = FALSE, total_col = NULL, scale_factor = NULL)
make_metacells_by_size(so_obj, group_vars, cell_size, cell_min, meta_vars,
                       assays, seed, save_membership, idents_col = NULL,
                       method = "random", reduction = "pca", n_dims = NULL,
                       refine = TRUE, n_init_pcs = 1, normalize = TRUE,
                       from_log = FALSE, mean_normalized = FALSE,
                       total_col = NULL, scale_factor = NULL)
make_metacells_by_num(so_obj, group_vars, n_metacells, min_cells, meta_vars,
                      assays, seed, save_membership, idents_col = NULL,
                      method = "random", reduction = "pca", n_dims = NULL,
                      refine = TRUE, n_init_pcs = 1, normalize = TRUE,
                       from_log = FALSE, mean_normalized = FALSE,
                       total_col = NULL, scale_factor = NULL)
```

R results hold the summed raw counts for each assay in `assays` and, with
`normalize = TRUE` (default), a log-normalised `data` layer. No genes are
scaled, so memory stays low on large objects; run `ScaleData()` on the genes
you need.

`idents_col`: name of any column in the result's metadata to set as the active
Seurat `Idents`. `NULL` (default) leaves Idents unset. The column can be a
pre-existing metadata column (e.g. `"cell_type"`) or one you derive after
aggregation:

All R functions return a named list:
- `$obj` — aggregated Seurat object
- `$membership` — data.frame mapping every input barcode to its metacell ID

---

## Traceability output

Every function writes (optionally) a CSV with one row per input cell:

```
cell_barcode, metacell_id, cell_type, sample, Treatment, ...
AAACCCAGTCCGAACC-1, mc3_A_s1, ExN, s1, KO, ...
AAACCCATCAGCTTGC-1, mc3_A_s1, ExN, s1, KO, ...
AAACGAAGTCCTGTAG-1, mc1_B_s2, InN, s2, WT, ...
AAAGGATAGCTCCATG-1, unassigned, ExN, s3, KO, ...   ← below cell_min threshold
```

The `metacell_id` column uses the format `mc{N}_{group_vars_joined}`.
Cells below the minimum threshold are marked `unassigned`.

---

## Choosing a strategy

**Pseudobulk** is the gold standard for DE testing (recommended when you have
≥ 4 biological replicates per condition). It produces one observation per
sample × cell-type and is unambiguous.

**Metacells by size** is best when samples have variable cell counts and you
want to maximise the number of metacells from well-represented groups while
naturally excluding sparse ones (below `cell_min`).

**Metacells by count** is best when you want equal replication across all
groups regardless of cell count, e.g. for downstream tools that assume balanced
designs.

### Random vs. similarity grouping

Both metacell strategies decide *how many* metacells each group gets and
*how big* they are the same way; `method` only changes *which* cells share a
metacell.

- **`random`** pools arbitrary cells from the group. Each metacell is a
  small random sample of the group, so metacells from one sample are
  near-exchangeable.
- **`similarity`** pools cells that are neighbours in the embedding, so each
  metacell represents one local cell state (finer than your cluster labels).
  Expression profiles are less noisy, and within-group heterogeneity (e.g. a
  sub-state or gradient inside a cell type) is preserved instead of averaged
  away.

How it works, within each group of *n* cells and *k* metacells:

1. Sort the cells by their score on the group's first principal component
   (computed on the embedding).
2. Cut the sorted list at rank `floor(n · k_left / k)`, where
   `k_left = k // 2`, so each side gets cells in proportion to the number of
   metacells it will hold. Polish the cut with a size-constrained 2-means
   step, so it follows gaps in the full embedding rather than only PC1.
3. Repeat on each side, recomputing PC1 within that side, until every piece
   is one metacell. This guarantees every metacell has `floor(n/k)` or
   `ceil(n/k)` cells.
4. Optionally (`refine=True`), re-split neighbouring metacells' pooled cells
   with their sizes fixed. This never increases within-metacell spread.

**`n_init_pcs`.** The 2-means polish in step 2 only moves the cut locally, so
it can stay on PC1 when PC1 is a broad, smooth gradient and the real gap
between sub-states lies along PC2 or PC3. With `n_init_pcs=3`, each split is
run starting from PC1, PC2 and PC3, and the cut with the lowest within-half
sum of squares is kept. The default `1` uses PC1 only. Runtime grows
sub-linearly in practice (≈2× for `3`).

The result is deterministic, and the Python and R implementations give
identical partitions for the same embedding.

> **Statistical note.** Metacells are not independent biological
> replicates under either method, and similarity grouping makes metacells
> from the same sample *more* different from one another. For differential
> expression, keep sample as a random effect (e.g. `dream`) or aggregate to
> true pseudobulk.

In all cases, expression is aggregated by **summing** raw counts.  If you are
using normalised/scaled values (SCT slot), pass `layer=` pointing to the raw
counts layer for biologically meaningful aggregation.

---

## License

MIT © 2026 Ethan Fenton
