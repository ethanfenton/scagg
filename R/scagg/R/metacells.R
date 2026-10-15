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

# Aggregate a Seurat object and attach clean metadata.
.aggregate_and_meta <- function(so_obj, group_vars, meta_vars, assays,
                                extra_meta = NULL) {
  mdata <- so_obj@meta.data

  # Filter out unassigned cells
  keep <- mdata$metacell_id != "unassigned"
  so_sub <- subset(so_obj, cells = rownames(mdata)[keep])

  # Aggregate expression
  mc_so <- Seurat::AggregateExpression(
    so_sub,
    group.by  = "metacell_id",
    assays    = assays,
    return.seurat = TRUE
  )

  # Build per-metacell metadata from the original object
  cols_to_keep <- unique(c("metacell_id", group_vars, meta_vars))
  cols_present <- intersect(cols_to_keep, colnames(mdata))

  mc_mdata <- mdata[keep, cols_present, drop = FALSE]

  # For each metacell take the most common value of each column
  mc_mdata_agg <- do.call(rbind, lapply(
    split(mc_mdata, mc_mdata$metacell_id),
    function(df) {
      row <- as.data.frame(lapply(df, function(x) {
        if (is.numeric(x)) mean(x, na.rm = TRUE)
        else {
          tab <- sort(table(x), decreasing = TRUE)
          if (length(tab) == 0L) NA else names(tab)[1L]
        }
      }), stringsAsFactors = FALSE)
      row$n_cells <- nrow(df)
      row
    }
  ))

  # Join to AggregateExpression metadata
  mc_mdata_agg$metacell_id <- rownames(mc_mdata_agg)
  if (!is.null(extra_meta)) {
    mc_mdata_agg <- dplyr::left_join(mc_mdata_agg, extra_meta, by = "metacell_id")
  }

  # Seurat v5's AggregateExpression replaces "_" with "-" in group names, so
  # map its column names back to the original metacell IDs before joining.
  agg_names <- rownames(mc_so@meta.data)
  orig_ids  <- mc_mdata_agg$metacell_id
  id_map    <- c(stats::setNames(orig_ids, orig_ids),
                 stats::setNames(orig_ids, gsub("_", "-", orig_ids, fixed = TRUE)))
  id_map    <- id_map[!duplicated(names(id_map))]
  mc_so@meta.data$metacell_id <- unname(id_map[agg_names])
  if (anyNA(mc_so@meta.data$metacell_id)) {
    warning("[scagg] some aggregated metacells could not be matched to their ",
            "metadata; their metadata columns will be NA.")
  }

  # Merge (left_join keeps row order, so the Seurat cell names stay valid)
  merged <- dplyr::left_join(
    mc_so@meta.data,
    mc_mdata_agg,
    by = "metacell_id",
    suffix = c("", ".scagg")
  )
  rownames(merged) <- agg_names
  mc_so@meta.data  <- merged

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
    n_init_pcs  = 1L
) {
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
          extra_meta = extra)
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
    n_init_pcs  = 1L
) {
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
          extra_meta = extra)
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
    idents_col  = NULL
) {
  if (is.null(meta_vars)) {
    meta_vars <- setdiff(colnames(so_obj@meta.data), group_vars)
  }

  mdata <- so_obj@meta.data
  mdata$metacell_id <- apply(
    mdata[, group_vars, drop = FALSE], 1,
    function(r) paste(r, collapse = "_")
  )
  so_obj@meta.data <- mdata

  .finish(so_obj, group_vars, meta_vars, assays, save_membership, mdata, idents_col)
}


# ── Shared finaliser ──────────────────────────────────────────────────────────

.finish <- function(so_obj, group_vars, meta_vars, assays,
                    save_membership, mdata, idents_col, extra_meta = NULL) {
  mc_so <- .aggregate_and_meta(so_obj, group_vars, meta_vars, assays, extra_meta)

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
