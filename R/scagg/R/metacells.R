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

# Aggregate a Seurat object and attach clean metadata.
.aggregate_and_meta <- function(so_obj, group_vars, meta_vars, assays) {
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

  # AggregateExpression uses "-" instead of "_" in cell names; normalise
  mc_so@meta.data$metacell_id <- rownames(mc_so@meta.data)

  # Merge
  merged <- dplyr::left_join(
    mc_so@meta.data,
    mc_mdata_agg,
    by = "metacell_id"
  )
  rownames(merged) <- merged$metacell_id
  mc_so@meta.data  <- merged

  mc_so
}


# ── Public functions ──────────────────────────────────────────────────────────

#' Aggregate cells into metacells of approximately \code{cell_size} cells
#'
#' Within each combination of \code{group_vars}, cells are shuffled then
#' assigned round-robin to \code{floor(n / cell_size)} metacells.  Groups
#' with fewer than \code{cell_min} cells are excluded.
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
#'   cell_size = 10, save_membership = "metacell_membership.csv"
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
    save_membership = NULL
) {
  if (is.null(cell_min)) cell_min <- .default_cell_min(cell_size)
  if (is.null(meta_vars)) {
    meta_vars <- setdiff(colnames(so_obj@meta.data), group_vars)
  }

  mdata <- so_obj@meta.data

  labels_fn <- function(idx) {
    n <- length(idx)
    if (n < cell_min) return(NULL)
    n_mc    <- if (n < cell_size) 1L else n %/% cell_size
    labels  <- rep(seq_len(n_mc), length.out = n)
    labels  <- sample(labels)
    key     <- paste(vapply(group_vars, function(v) as.character(mdata[idx[1], v]),
                            character(1)), collapse = "_")
    paste0("mc", labels, "_", key)
  }

  mdata <- .assign_metacells(mdata, group_vars, labels_fn, seed)
  so_obj@meta.data <- mdata

  .finish(so_obj, group_vars, meta_vars, assays, save_membership, mdata)
}


#' Aggregate cells into exactly \code{n_metacells} metacells per group
#'
#' Within each combination of \code{group_vars}, cells are randomly assigned
#' to exactly \code{n_metacells} groups.
#'
#' @param so_obj A Seurat object.
#' @param group_vars Character vector of metadata column names to group by.
#' @param n_metacells Number of metacells to create per group.
#' @param min_cells Minimum cells required. Default: equal to \code{n_metacells}.
#' @param meta_vars Metadata columns to carry forward.
#' @param assays Seurat assay(s) to aggregate.
#' @param seed Random seed.
#' @param save_membership Path for membership CSV.
#'
#' @return Named list with \code{obj} (Seurat) and \code{membership} (data.frame).
#'
#' @examples
#' \dontrun{
#' res <- make_metacells_by_num(
#'   so, group_vars = c("cell_type", "sample"),
#'   n_metacells = 10
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
    save_membership = NULL
) {
  if (is.null(meta_vars)) {
    meta_vars <- setdiff(colnames(so_obj@meta.data), group_vars)
  }

  mdata <- so_obj@meta.data

  labels_fn <- function(idx) {
    n <- length(idx)
    if (n < min_cells) return(NULL)
    actual_n <- min(n_metacells, n)
    labels   <- rep(seq_len(actual_n), length.out = n)
    labels   <- sample(labels)
    key      <- paste(vapply(group_vars, function(v) as.character(mdata[idx[1], v]),
                             character(1)), collapse = "_")
    paste0("mc", labels, "_", key)
  }

  mdata <- .assign_metacells(mdata, group_vars, labels_fn, seed)
  so_obj@meta.data <- mdata

  .finish(so_obj, group_vars, meta_vars, assays, save_membership, mdata)
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
#'
#' @return Named list with \code{obj} and \code{membership}.
#'
#' @examples
#' \dontrun{
#' res <- make_pseudobulk(
#'   so, group_vars = c("sample", "cell_type"),
#'   assays = c("SCT", "RNA")
#' )
#' }
#'
#' @export
make_pseudobulk <- function(
    so_obj,
    group_vars,
    meta_vars   = NULL,
    assays      = "SCT",
    save_membership = NULL
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

  .finish(so_obj, group_vars, meta_vars, assays, save_membership, mdata)
}


# ── Shared finaliser ──────────────────────────────────────────────────────────

.finish <- function(so_obj, group_vars, meta_vars, assays, save_membership, mdata) {
  mc_so <- .aggregate_and_meta(so_obj, group_vars, meta_vars, assays)

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
