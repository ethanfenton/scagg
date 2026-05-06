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

Available for **Python (AnnData/scanpy)** and **R (Seurat)**.
Full cell → metacell traceability is built in.

[![CI](https://github.com/ethanfenton/scagg/actions/workflows/ci.yml/badge.svg)](https://github.com/ethanfenton/scagg/actions)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)

---

## Installation

### Python

```bash
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
scagg.make_pseudobulk(adata, group_vars, *, meta_vars, layer, save_membership)
scagg.make_metacells_by_size(adata, group_vars, cell_size, *, cell_min, meta_vars,
                              layer, seed, save_membership)
scagg.make_metacells_by_num(adata, group_vars, n_metacells, *, min_cells, meta_vars,
                             layer, seed, save_membership)
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
| `seed` | 42 | RNG seed for reproducibility |
| `save_membership` | None | CSV path for traceability output |

### R

```
make_pseudobulk(so_obj, group_vars, meta_vars, assays, save_membership)
make_metacells_by_size(so_obj, group_vars, cell_size, cell_min, meta_vars,
                       assays, seed, save_membership)
make_metacells_by_num(so_obj, group_vars, n_metacells, min_cells, meta_vars,
                      assays, seed, save_membership)
```

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

In all cases, expression is aggregated by **summing** raw counts.  If you are
using normalised/scaled values (SCT slot), pass `layer=` pointing to the raw
counts layer for biologically meaningful aggregation.

---

## Migrating from your existing code

If you were using `make_metacells_by_size(so_obj, cell_size=10)` with
hardcoded `ct3`/`sample` columns:

```r
# Old (hardcoded)
metacell_so <- make_metacells_by_size(so_obj, cell_size = 10)

# New (scagg, generalised)
res <- make_metacells_by_size(
  so_obj,
  group_vars = c("ct3", "sample"),
  cell_size  = 10,
  meta_vars  = c("Treatment", "Timepoint", "sex", "Prepper"),
  assays     = "SCT"
)
metacell_so <- res$obj
# Add any derived columns you need:
metacell_so$time_treat <- paste(metacell_so$Timepoint, metacell_so$Treatment, sep = "_")
```

---

## License

MIT © 2026 Ethan Fenton
