# ── Internal helpers ──────────────────────────────────────────────────────────

.default_cell_min <- function(cell_size) {
  if (cell_size == 10L) return(7L)
  if (cell_size == 25L) return(15L)
  floor(cell_size * 0.7)
}

# Assign metacell IDs into mdata$metacell_id.
# Returns updated mdata data.frame and the number of unassigned cells.
.assign_metacells <- function(mdata, group_vars, labels_fn, seed = 42L) {
  set.seed(seed)
  mdata$metacell_id <- "unassigned"

  combos <- unique(mdata[, group_vars, drop = FALSE])
  for (i in seq_len(nrow(combos))) {
    cond <- rep(TRUE, nrow(mdata))
    for (v in group_vars) cond <- cond & (mdata[[v]] == combos[i, v])
    idx <- which(cond)
    result <- labels_fn(idx)
    if (!is.null(result)) mdata$metacell_id[idx] <- result
  }
  .check_metacell_purity(mdata, group_vars)
  mdata
}

# Validate `method`; return the embedding matrix (rows = cells) or NULL.
.get_embedding <- function(so_obj, method, reduction, n_dims) {
  if (method == "random") return(NULL)
  if (!reduction %in% Seurat::Reductions(so_obj)) {
    stop(sprintf(
      "[scagg] reduction '%s' not found (available: %s). Run e.g. RunPCA() first.",
      reduction, paste(Seurat::Reductions(so_obj), collapse = ", ")))
  }
  emb <- Seurat::Embeddings(so_obj, reduction = reduction)[colnames(so_obj), , drop = FALSE]
  if (!is.null(n_dims)) {
    if (n_dims < 1L) stop("[scagg] n_dims must be >= 1")
    emb <- emb[, seq_len(min(n_dims, ncol(emb))), drop = FALSE]
  }
  if (any(!is.finite(emb))) stop("[scagg] embedding contains NA or infinite values")
  emb
}

# Metacell labels 1..n_mc for the cells `idx` of one group.
.group_labels <- function(idx, n_mc, emb, refine, n_init_pcs) {
  if (is.null(emb)) return(sample(rep(seq_len(n_mc), length.out = length(idx))))
  .balanced_partition(emb[idx, , drop = FALSE], n_mc, refine = refine,
                      n_init_pcs = n_init_pcs)
}

# Per-metacell compactness, plus a log line comparing it to random metacells
# of the same sizes.  Returns data.frame(metacell_id, compactness).
.compactness_table <- function(mdata, emb, group_vars, seed) {
  keep <- mdata$metacell_id != "unassigned"
  ids  <- mdata$metacell_id[keep]
  sub  <- emb[keep, , drop = FALSE]
  comp <- .compactness(sub, ids)

  set.seed(seed)
  strata <- interaction(mdata[keep, group_vars, drop = FALSE], drop = TRUE)
  shuffled <- ids
  for (rows in split(seq_along(ids), strata)) shuffled[rows] <- sample(ids[rows])
  message(sprintf(
    "[scagg] median metacell compactness: %.3g (random metacells of the same sizes: %.3g)",
    stats::median(comp), stats::median(.compactness(sub, shuffled))))

  data.frame(metacell_id = names(comp), compactness = unname(comp),
             stringsAsFactors = FALSE)
}

# Stop unless every metacell holds cells from a single combination of
# `group_vars` (e.g. one sample × one cell type).  Different groups can map to
# the same ID when values contain "_" (("A_1", "x") and ("A", "1_x") both give
# "A_1_x"), which would otherwise sum cells from different samples together.
.check_metacell_purity <- function(mdata, group_vars) {
  keep <- mdata$metacell_id != "unassigned"
  if (!any(keep)) return(invisible(NULL))
  key <- do.call(paste, c(lapply(mdata[keep, group_vars, drop = FALSE], as.character),
                          sep = "\r"))
  pairs <- unique(data.frame(id = mdata$metacell_id[keep], key = key,
                             stringsAsFactors = FALSE))
  bad <- unique(pairs$id[duplicated(pairs$id)])
  if (length(bad) > 0L) {
    stop(sprintf(paste0(
      "[scagg] %d metacell ID(s) would combine cells from different %s groups ",
      "(e.g. '%s'), because group values joined with '_' collide. Rename the ",
      "values in %s so that no joined combination is ambiguous."),
      length(bad), paste(group_vars, collapse = " × "), bad[1],
      paste(group_vars, collapse = ", ")))
  }
  invisible(NULL)
}

.check_value_args <- function(from_log, mean_normalized) {
  if (isTRUE(mean_normalized) && !isTRUE(from_log)) {
    stop("[scagg] mean_normalized = TRUE requires from_log = TRUE (it averages ",
         "log-normalised values in linear space).")
  }
}

# Every layer of an assay named `layer` (Seurat v5 splits it into one per
# batch: counts.s1, counts.s2, ...), as a list of features × cells matrices.
.assay_layers <- function(assay_obj, assay, layer) {
  mats <- if (inherits(assay_obj, "Assay5")) {
    lapply(SeuratObject::Layers(assay_obj, search = layer),
           function(l) SeuratObject::LayerData(assay_obj, layer = l))
  } else {
    list(methods::slot(assay_obj, layer))
  }
  mats <- Filter(function(m) prod(dim(m)) > 0, mats)
  if (length(mats) == 0L) stop(sprintf("[scagg] assay '%s' has no '%s' layer", assay, layer))
  mats
}

# expm1 of a matrix, as a dgCMatrix.
.expm1_sparse <- function(m) {
  m <- methods::as(methods::as(m, "CsparseMatrix"), "dMatrix")
  m@x <- expm1(m@x)
  m
}

# Smallest positive value in each column of a dgCMatrix (NA if none).
.col_min_positive <- function(m) {
  col <- rep.int(seq_len(ncol(m)), diff(m@p))
  pos <- m@x > 0
  mins <- tapply(m@x[pos], factor(col[pos], levels = seq_len(ncol(m))), min)
  as.vector(mins)
}

# Per-column value of a single count in linear-space data (dgCMatrix, genes
# × cells).  Every non-zero value is a whole multiple of it, so the unit is
# min / c for the smallest c in 1..max_multiple that makes all of the cell's
# values whole multiples.  c = 1 is the plain smallest-non-zero trick; larger
# c handles cells with no count of exactly 1.  NA for cells without counts.
.count_unit <- function(m, max_multiple = 10L, tol = 0.05) {
  min_pos <- .col_min_positive(m)
  unit <- min_pos
  todo <- is.finite(min_pos)
  col  <- rep.int(seq_len(ncol(m)), diff(m@p))
  colf <- factor(col, levels = seq_len(ncol(m)))
  for (c in seq_len(max_multiple)) {
    if (!any(todo)) break
    u     <- min_pos / c
    ratio <- m@x / u[col]
    dev   <- ifelse(todo[col], abs(ratio - round(ratio)), 0)
    worst <- as.vector(tapply(dev, colf, max))
    ok    <- todo & !is.na(worst) & worst < tol
    unit[ok] <- u[ok]
    todo <- todo & !ok
  }
  unit
}

# Raw counts from log-normalised layers, x = log1p(count / total * sf):
# count = expm1(x) * total / sf.  Without totals, each cell's smallest
# non-zero value is taken as a count of 1.  `scale_factor = NULL` infers sf.
.recover_counts <- function(mats, totals, scale_factor) {
  lin  <- lapply(mats, .expm1_sparse)
  mins <- lapply(lin, .count_unit)        # value of one count in each cell

  if (is.null(totals)) {
    message("[scagg] from_log: no per-cell totals found; falling back to the ",
            "smallest-non-zero trick (each cell's values are scaled so its ",
            "smallest non-zero value is the smallest whole count consistent ",
            "with the rest, usually 1). Set total_col to use stored totals.")
    sizes <- mins
  } else {
    tot <- lapply(lin, function(m) unname(totals[colnames(m)]))
    if (is.null(scale_factor)) {
      est <- unlist(Map(`*`, mins, tot))
      est <- est[is.finite(est) & est > 0]
      if (length(est) == 0L) {
        stop("[scagg] cannot infer scale_factor: no cell has both a valid total ",
             "and a non-zero value")
      }
      scale_factor <- stats::median(est)
      message(sprintf("[scagg] from_log: inferred scale_factor = %.6g", scale_factor))
    }
    n_bad <- 0L
    sizes <- Map(function(t, mn) {
      ok <- is.finite(t) & t > 0
      n_bad <<- n_bad + sum(!ok)
      ifelse(ok, scale_factor / t, mn)
    }, tot, mins)
    if (n_bad > 0L) {
      message(sprintf(paste0(
        "[scagg] from_log: %d cell(s) have a missing or zero total; falling back ",
        "to the smallest-non-zero trick for them."), n_bad))
    }
  }

  dev <- 0
  out <- Map(function(m, size) {
    size[!is.finite(size) | size <= 0] <- 1          # cells with no counts
    col <- rep.int(seq_len(ncol(m)), diff(m@p))
    m@x <- m@x / size[col]
    if (length(m@x)) dev <<- max(dev, abs(m@x - round(m@x)))
    m@x <- round(m@x)
    m
  }, lin, sizes)
  if (dev > 0.1) {
    warning(sprintf(paste0(
      "[scagg] recovered counts are up to %.3g away from whole numbers, so the ",
      "input may not be log1p(count / total * scale_factor) (e.g. it was scaled ",
      "or used a different total). Check total_col and scale_factor."), dev))
  }
  out
}

# Sum each matrix's columns per metacell: features × metacells dgCMatrix.
# Each matrix is multiplied by a sparse cells × metacells indicator, so
# memory scales with the non-zeros, not with features × cells.
.sum_by_group <- function(mats, groups) {
  total <- NULL
  for (m in mats) {
    g <- groups[colnames(m)]
    rows <- which(!is.na(g))
    ind <- Matrix::sparseMatrix(
      i = rows, j = as.integer(g[rows]), x = 1,
      dims = c(ncol(m), nlevels(groups)))
    part <- methods::as(m %*% ind, "CsparseMatrix")
    rownames(part) <- rownames(m)
    if (is.null(total)) {
      total <- part
    } else {
      if (!identical(rownames(part), rownames(total))) {
        stop("[scagg] layers of one assay have different features")
      }
      total <- total + part
    }
  }
  colnames(total) <- levels(groups)
  total
}

# Aggregate one assay.  Returns list(values, kind): kind "counts" (summed
# counts, or recovered counts with from_log) or "mean_normalized"
# (log1p of the mean of expm1(data)).
.aggregate_assay <- function(so_obj, assay, groups, from_log, mean_normalized,
                             total_col, scale_factor) {
  assay_obj <- so_obj[[assay]]
  if (!from_log) {
    return(list(values = .sum_by_group(.assay_layers(assay_obj, assay, "counts"), groups),
                kind = "counts"))
  }
  if (inherits(assay_obj, "SCTAssay")) {
    stop("[scagg] from_log = TRUE is not needed for SCT assays: their 'counts' ",
         "layer already holds the corrected counts. Use from_log = FALSE.")
  }
  mats <- .assay_layers(assay_obj, assay, "data")

  if (mean_normalized) {
    sums  <- .sum_by_group(lapply(mats, .expm1_sparse), groups)
    sizes <- as.vector(table(groups))
    sums@x <- log1p(sums@x / sizes[rep.int(seq_len(ncol(sums)), diff(sums@p))])
    return(list(values = sums, kind = "mean_normalized"))
  }

  md <- so_obj@meta.data
  if (is.null(total_col)) {
    total_col <- paste0("nCount_", assay)
    totals <- if (total_col %in% colnames(md)) md[[total_col]] else NULL
    if (!is.null(totals)) message("[scagg] from_log: using per-cell totals from '",
                                  total_col, "'")
  } else {
    if (!total_col %in% colnames(md)) {
      stop(sprintf("[scagg] total_col '%s' not found in metadata", total_col))
    }
    totals <- md[[total_col]]
  }
  if (!is.null(totals)) {
    totals <- stats::setNames(suppressWarnings(as.numeric(totals)), rownames(md))
  }
  list(values = .sum_by_group(.recover_counts(mats, totals, scale_factor), groups),
       kind = "counts")
}

# Per-metacell metadata: mean for numeric columns, most common value otherwise
# (ties go to the first value in sort order).
.metacell_meta <- function(mdata, cols, groups) {
  out <- data.frame(row.names = levels(groups))
  for (col in cols) {
    x <- mdata[[col]]
    if (is.numeric(x)) {
      ok <- !is.na(x)
      out[[col]] <- as.vector(
        rowsum(ifelse(ok, x, 0), groups, reorder = TRUE) /
          rowsum(as.numeric(ok), groups, reorder = TRUE))
    } else {
      x  <- as.character(x)
      ok <- !is.na(x)
      gi <- as.integer(groups)[ok]
      xv <- x[ok]
      o  <- order(gi, xv)
      gi <- gi[o]; xv <- xv[o]
      n  <- length(gi)
      starts <- if (n == 0L) logical(0) else
        c(TRUE, gi[-1L] != gi[-n] | xv[-1L] != xv[-n])
      run_len <- tabulate(cumsum(starts))
      run_g   <- gi[starts]
      run_x   <- xv[starts]
      best    <- order(run_g, -run_len)
      first   <- best[!duplicated(run_g[best])]
      vals    <- rep(NA_character_, nlevels(groups))
      vals[run_g[first]] <- run_x[first]
      out[[col]] <- vals
    }
  }
  out$n_cells <- as.vector(table(groups))
  out
}

# Build the aggregated Seurat object: summed counts for each assay, plus
# per-metacell metadata.  Counts are summed directly rather than through
# Seurat::AggregateExpression(return.seurat = TRUE), which also copies the
# object with subset() and runs ScaleData() on every gene; that needs more than
# 32 GB on large objects.  With `normalize = TRUE` each assay is
# log-normalised (sparse and cheap); nothing is scaled.  With `from_log` the
# `data` layer is used instead of `counts` (see .aggregate_assay).
.aggregate_and_meta <- function(so_obj, group_vars, meta_vars, assays,
                                extra_meta = NULL, normalize = TRUE,
                                from_log = FALSE, mean_normalized = FALSE,
                                total_col = NULL, scale_factor = NULL) {
  mdata <- so_obj@meta.data
  keep  <- mdata$metacell_id != "unassigned"
  if (!any(keep)) stop("[scagg] all cells were excluded; no metacells formed")

  ids    <- mdata$metacell_id[keep]
  groups <- factor(ids, levels = unique(ids))
  cell_groups <- stats::setNames(factor(rep(NA_character_, nrow(mdata)),
                                        levels = levels(groups)),
                                 rownames(mdata))
  cell_groups[keep] <- groups

  missing <- setdiff(assays, SeuratObject::Assays(so_obj))
  if (length(missing) > 0L) {
    stop(sprintf("[scagg] assay(s) not found: %s (available: %s)",
                 paste(missing, collapse = ", "),
                 paste(SeuratObject::Assays(so_obj), collapse = ", ")))
  }

  cols <- setdiff(intersect(unique(c(group_vars, meta_vars)), colnames(mdata)),
                  "metacell_id")
  mc_meta <- .metacell_meta(mdata[keep, , drop = FALSE], cols, groups)
  mc_meta <- cbind(data.frame(metacell_id = levels(groups),
                              row.names = levels(groups),
                              stringsAsFactors = FALSE),
                   mc_meta)
  if (!is.null(extra_meta)) {
    mc_meta <- dplyr::left_join(mc_meta, extra_meta, by = "metacell_id")
    rownames(mc_meta) <- mc_meta$metacell_id
  }

  mc_so <- NULL
  for (assay in assays) {
    agg <- .aggregate_assay(so_obj, assay, cell_groups, from_log, mean_normalized,
                            total_col, scale_factor)
    first <- is.null(mc_so)
    v5 <- first || inherits(mc_so[[assays[1L]]], "Assay5")
    new_assay <- suppressWarnings(if (agg$kind == "mean_normalized") {
      if (v5) SeuratObject::CreateAssay5Object(data = agg$values)
      else SeuratObject::CreateAssayObject(data = agg$values)
    } else {
      if (v5) SeuratObject::CreateAssay5Object(counts = agg$values)
      else SeuratObject::CreateAssayObject(counts = agg$values)
    })
    if (first) {
      mc_so <- SeuratObject::CreateSeuratObject(new_assay, assay = assay)
    } else {
      mc_so[[assay]] <- new_assay
    }
    if (normalize && agg$kind == "counts") {
      mc_so <- Seurat::NormalizeData(mc_so, assay = assay, verbose = FALSE)
    }
  }

  # Keep Seurat's own columns (orig.ident, nCount_*, nFeature_*) and add ours.
  own <- mc_so@meta.data[, setdiff(colnames(mc_so@meta.data), colnames(mc_meta)),
                         drop = FALSE]
  mc_so@meta.data <- cbind(own, mc_meta[colnames(mc_so), , drop = FALSE])
  mc_so
}


# ── Public functions ──────────────────────────────────────────────────────────

#' Aggregate cells into metacells of approximately \code{cell_size} cells
#'
#' Within each combination of \code{group_vars}, cells are assigned to
#' \code{floor(n / cell_size)} metacells, either randomly or by similarity in
#' an embedding (see \code{method}).  Groups with fewer than \code{cell_min}
#' cells are excluded.  With \code{method = "similarity"} the metacell
#' metadata gains a \code{compactness} column (mean cell-to-centroid distance
#' in the embedding).
#'
#' @param so_obj A Seurat object.
#' @param group_vars Character vector of metadata column names to group by,
#'   e.g. \code{c("cell_type", "sample")}.
#' @param cell_size Target cells per metacell.
#' @param cell_min Minimum cells required to form at least one metacell.
#'   Default: 70\% of \code{cell_size} (7 for size 10, 15 for size 25).
#' @param meta_vars Metadata columns to carry into the metacell object.
#'   Default: all columns not in \code{group_vars}.
#' @param assays Seurat assay(s) to aggregate (default \code{"SCT"}).
#' @param seed Integer random seed (default 42).
#' @param save_membership If non-NULL, path to save a CSV mapping each cell
#'   barcode to its metacell ID.
#' @param idents_col Name of the metadata column to set as the active
#'   \code{Idents} of the returned Seurat object (e.g. \code{"cell_type"} or
#'   a derived column like \code{"time_treat"}).  \code{NULL} (default) leaves
#'   Idents unchanged.
#' @param method \code{"random"} (default): cells are shuffled and dealt
#'   round-robin into metacells.  \code{"similarity"}: cells that are close in
#'   the embedding \code{reduction} are grouped together via balanced
#'   recursive PCA bisection.  Both give identical metacell counts and sizes
#'   (within one cell); only which cells share a metacell differs.  The
#'   similarity method is deterministic, so \code{seed} does not affect it.
#' @param reduction Dimensional reduction used when
#'   \code{method = "similarity"}, e.g. \code{"pca"}, \code{"harmony"}.
#' @param n_dims Use only the first \code{n_dims} dimensions of
#'   \code{reduction}.  \code{NULL} (default) uses all.
#' @param refine When \code{method = "similarity"}, polish the bisection with
#'   size-preserving swaps between neighbouring metacells.
#' @param n_init_pcs When \code{method = "similarity"}, try each split
#'   starting from PC1 ... PC\code{n_init_pcs} and keep the tightest cut.
#'   \code{1} (default) starts from PC1 only; \code{3} helps when a broad
#'   gradient dominates PC1 but distinct sub-states separate along a later
#'   PC.  Cost grows roughly linearly with \code{n_init_pcs}.
#' @param normalize If \code{TRUE} (default), log-normalise each assay of the
#'   result with \code{Seurat::NormalizeData()}.  Counts are always the raw
#'   sums; no gene is scaled (run \code{ScaleData()} on the genes you need).
#' @param from_log If \code{TRUE}, aggregate the log-normalised \code{data}
#'   layer instead of \code{counts}: values \code{log1p(count / total *
#'   scale_factor)} are converted back to raw counts with
#'   \code{count = expm1(x) * total / scale_factor} before summing, using each
#'   cell's total from \code{total_col}.  If no totals are found, each
#'   cell's smallest non-zero value is taken as a count of 1 (a message says
#'   so).  Not needed for SCT assays.  Default \code{FALSE}.
#' @param mean_normalized Only with \code{from_log = TRUE}.  Instead of
#'   recovering counts, store \code{log1p(mean(expm1(x)))} per metacell in the
#'   \code{data} layer (no \code{counts} layer): the average normalised
#'   profile, with every cell weighted equally.  Not for DESeq2 / edgeR.
#' @param total_col Metadata column with each cell's total count.
#'   \code{NULL} (default) uses \code{nCount_<assay>} if present.
#' @param scale_factor Target sum used in normalisation (Seurat default 1e4).
#'   \code{NULL} (default) infers it from the data.
#'
#' @details Cells are only ever pooled with cells that share every value of
#'   \code{group_vars}: put your sample column in \code{group_vars} (e.g.
#'   \code{c("cell_type", "sample")}) so that no metacell mixes samples.  The
#'   function stops if two groups would map to the same metacell ID.
#'
#' @return A named list with:
#'   \describe{
#'     \item{\code{obj}}{The aggregated Seurat object.}
#'     \item{\code{membership}}{A data frame mapping every input cell to its
#'       metacell ID (or "unassigned").}
#'   }
#'
#' @examples
#' \dontrun{
#' res <- make_metacells_by_size(
#'   so, group_vars = c("cell_type", "sample"),
#'   cell_size = 10, idents_col = "cell_type",
#'   save_membership = "metacell_membership.csv"
#' )
#' seurat_mc <- res$obj
#' }
#'
#' @export
make_metacells_by_size <- function(
    so_obj,
    group_vars,
    cell_size,
    cell_min    = NULL,
    meta_vars   = NULL,
    assays      = "SCT",
    seed        = 42L,
    save_membership = NULL,
    idents_col  = NULL,
    method      = c("random", "similarity"),
    reduction   = "pca",
    n_dims      = NULL,
    refine      = TRUE,
    n_init_pcs  = 1L,
    normalize   = TRUE,
    from_log    = FALSE,
    mean_normalized = FALSE,
    total_col   = NULL,
    scale_factor = NULL
) {
  .check_value_args(from_log, mean_normalized)
  method <- match.arg(method)
  if (n_init_pcs < 1L) stop("[scagg] n_init_pcs must be >= 1")
  if (is.null(cell_min)) cell_min <- .default_cell_min(cell_size)
  if (is.null(meta_vars)) {
    meta_vars <- setdiff(colnames(so_obj@meta.data), group_vars)
  }

  emb   <- .get_embedding(so_obj, method, reduction, n_dims)
  mdata <- so_obj@meta.data

  labels_fn <- function(idx) {
    n <- length(idx)
    if (n < cell_min) return(NULL)
    n_mc    <- if (n < cell_size) 1L else n %/% cell_size
    labels  <- .group_labels(idx, n_mc, emb, refine, n_init_pcs)
    key     <- paste(vapply(group_vars, function(v) as.character(mdata[idx[1], v]),
                            character(1)), collapse = "_")
    paste0("mc", labels, "_", key)
  }

  mdata <- .assign_metacells(mdata, group_vars, labels_fn, seed)
  so_obj@meta.data <- mdata

  extra <- if (is.null(emb)) NULL else .compactness_table(mdata, emb, group_vars, seed)
  .finish(so_obj, group_vars, meta_vars, assays, save_membership, mdata, idents_col,
          extra_meta = extra, normalize = normalize,
          values = list(from_log = from_log, mean_normalized = mean_normalized,
                        total_col = total_col, scale_factor = scale_factor))
}


#' Aggregate cells into exactly \code{n_metacells} metacells per group
#'
#' Within each combination of \code{group_vars}, cells are assigned to
#' exactly \code{n_metacells} groups, either randomly or by similarity in an
#' embedding (see \code{method}).
#'
#' @param so_obj A Seurat object.
#' @param group_vars Character vector of metadata column names to group by.
#' @param n_metacells Number of metacells to create per group.
#' @param min_cells Minimum cells required. Default: equal to \code{n_metacells}.
#' @param meta_vars Metadata columns to carry forward.
#' @param assays Seurat assay(s) to aggregate.
#' @param seed Random seed.
#' @param save_membership Path for membership CSV.
#' @param idents_col Metadata column to set as active Idents. \code{NULL}
#'   (default) leaves Idents unchanged.
#' @inheritParams make_metacells_by_size
#' @inherit make_metacells_by_size details
#'
#' @return Named list with \code{obj} (Seurat) and \code{membership} (data.frame).
#'
#' @examples
#' \dontrun{
#' res <- make_metacells_by_num(
#'   so, group_vars = c("cell_type", "sample"),
#'   n_metacells = 10, idents_col = "cell_type"
#' )
#' }
#'
#' @export
make_metacells_by_num <- function(
    so_obj,
    group_vars,
    n_metacells,
    min_cells   = n_metacells,
    meta_vars   = NULL,
    assays      = "SCT",
    seed        = 42L,
    save_membership = NULL,
    idents_col  = NULL,
    method      = c("random", "similarity"),
    reduction   = "pca",
    n_dims      = NULL,
    refine      = TRUE,
    n_init_pcs  = 1L,
    normalize   = TRUE,
    from_log    = FALSE,
    mean_normalized = FALSE,
    total_col   = NULL,
    scale_factor = NULL
) {
  .check_value_args(from_log, mean_normalized)
  method <- match.arg(method)
  if (n_init_pcs < 1L) stop("[scagg] n_init_pcs must be >= 1")
  if (is.null(meta_vars)) {
    meta_vars <- setdiff(colnames(so_obj@meta.data), group_vars)
  }

  emb   <- .get_embedding(so_obj, method, reduction, n_dims)
  mdata <- so_obj@meta.data

  labels_fn <- function(idx) {
    n <- length(idx)
    if (n < min_cells) return(NULL)
    actual_n <- min(n_metacells, n)
    labels   <- .group_labels(idx, actual_n, emb, refine, n_init_pcs)
    key      <- paste(vapply(group_vars, function(v) as.character(mdata[idx[1], v]),
                             character(1)), collapse = "_")
    paste0("mc", labels, "_", key)
  }

  mdata <- .assign_metacells(mdata, group_vars, labels_fn, seed)
  so_obj@meta.data <- mdata

  extra <- if (is.null(emb)) NULL else .compactness_table(mdata, emb, group_vars, seed)
  .finish(so_obj, group_vars, meta_vars, assays, save_membership, mdata, idents_col,
          extra_meta = extra, normalize = normalize,
          values = list(from_log = from_log, mean_normalized = mean_normalized,
                        total_col = total_col, scale_factor = scale_factor))
}


#' True pseudobulk aggregation: one sample per unique group combination
#'
#' Every cell is included; there is no random assignment.  Each unique
#' combination of \code{group_vars} values becomes one pseudobulk sample.
#'
#' @param so_obj A Seurat object.
#' @param group_vars Metadata column(s) to group by, e.g.
#'   \code{c("sample", "cell_type")}.
#' @param meta_vars Metadata columns to carry forward.
#' @param assays Seurat assay(s) to aggregate.
#' @param save_membership Path for membership CSV.
#' @param idents_col Metadata column to set as active Idents. \code{NULL}
#'   (default) leaves Idents unchanged.
#' @inheritParams make_metacells_by_size
#'
#' @return Named list with \code{obj} and \code{membership}.
#'
#' @examples
#' \dontrun{
#' res <- make_pseudobulk(
#'   so, group_vars = c("sample", "cell_type"),
#'   assays = c("SCT", "RNA"), idents_col = "sample"
#' )
#' }
#'
#' @export
make_pseudobulk <- function(
    so_obj,
    group_vars,
    meta_vars   = NULL,
    assays      = "SCT",
    save_membership = NULL,
    idents_col  = NULL,
    normalize   = TRUE,
    from_log    = FALSE,
    mean_normalized = FALSE,
    total_col   = NULL,
    scale_factor = NULL
) {
  .check_value_args(from_log, mean_normalized)
  if (is.null(meta_vars)) {
    meta_vars <- setdiff(colnames(so_obj@meta.data), group_vars)
  }

  mdata <- so_obj@meta.data
  mdata$metacell_id <- apply(
    mdata[, group_vars, drop = FALSE], 1,
    function(r) paste(r, collapse = "_")
  )
  .check_metacell_purity(mdata, group_vars)
  so_obj@meta.data <- mdata

  .finish(so_obj, group_vars, meta_vars, assays, save_membership, mdata, idents_col,
          normalize = normalize,
          values = list(from_log = from_log, mean_normalized = mean_normalized,
                        total_col = total_col, scale_factor = scale_factor))
}


# ── Shared finaliser ──────────────────────────────────────────────────────────

.finish <- function(so_obj, group_vars, meta_vars, assays,
                    save_membership, mdata, idents_col, extra_meta = NULL,
                    normalize = TRUE, values = list()) {
  mc_so <- .aggregate_and_meta(
    so_obj, group_vars, meta_vars, assays, extra_meta, normalize,
    from_log        = isTRUE(values$from_log),
    mean_normalized = isTRUE(values$mean_normalized),
    total_col       = values$total_col,
    scale_factor    = values$scale_factor)

  # Optionally set active Idents
  if (!is.null(idents_col)) {
    if (!idents_col %in% colnames(mc_so@meta.data)) {
      warning(sprintf(
        "[scagg] idents_col '%s' not found in metacell metadata; Idents not set.",
        idents_col
      ))
    } else {
      Seurat::Idents(mc_so) <- idents_col
      message("[scagg] Idents set to '", idents_col, "'")
    }
  }

  # Build membership table
  membership <- mdata
  membership$cell_barcode <- rownames(mdata)
  membership <- membership[, unique(c("cell_barcode", "metacell_id",
                                      group_vars, meta_vars)),
                           drop = FALSE]

  if (!is.null(save_membership)) {
    dir.create(dirname(save_membership), showWarnings = FALSE, recursive = TRUE)
    utils::write.csv(membership, file = save_membership, row.names = FALSE)
    message("[scagg] Membership saved → ", save_membership)
  }

  n_mc <- sum(membership$metacell_id != "unassigned" &
              !duplicated(membership$metacell_id))
  n_excl <- sum(membership$metacell_id == "unassigned")
  message(sprintf("[scagg] %d cells → %d metacells (%d excluded)",
                  nrow(membership), n_mc, n_excl))

  list(obj = mc_so, membership = membership)
}
