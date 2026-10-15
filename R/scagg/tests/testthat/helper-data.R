# 2 samples × 4 well-separated subclusters × 30 cells, with a "pca" reduction.
make_clustered_seurat <- function() {
  set.seed(1)
  n_sub <- 4L; per <- 30L; nd <- 10L
  centers <- matrix(stats::rnorm(n_sub * nd, 0, 10), n_sub)
  sub <- rep(rep(seq_len(n_sub), each = per), 2L)
  emb <- centers[sub, ] + matrix(stats::rnorm(length(sub) * nd), length(sub))
  counts <- matrix(stats::rpois(20L * length(sub), 5), 20L,
                   dimnames = list(paste0("g", 1:20), paste0("c", seq_along(sub))))
  so <- Seurat::CreateSeuratObject(counts)
  so$cell_type  <- "T"
  so$sample     <- rep(c("s1", "s2"), each = n_sub * per)
  so$subcluster <- as.character(sub)
  rownames(emb) <- colnames(so)
  colnames(emb) <- paste0("PC_", seq_len(nd))
  so[["pca"]] <- Seurat::CreateDimReducObject(embeddings = emb, key = "PC_",
                                              assay = "RNA")
  so
}

purity <- function(membership) {
  m <- membership[membership$metacell_id != "unassigned", ]
  mean(tapply(m$subcluster, m$metacell_id, function(s) max(table(s)) / length(s)))
}

quiet <- function(expr) suppressMessages(suppressWarnings(expr))
