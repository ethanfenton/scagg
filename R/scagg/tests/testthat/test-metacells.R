skip_if_not_installed("Seurat")

test_that("metacell metadata is populated (Seurat v5 '_' → '-' renaming)", {
  res <- quiet(make_metacells_by_size(make_clustered_seurat(), c("cell_type", "sample"),
                                      10, assays = "RNA"))
  md <- res$obj@meta.data
  expect_false(anyNA(md$n_cells))
  expect_false(anyNA(md$sample))
  expect_true(all(md$metacell_id %in% res$membership$metacell_id))
  expect_identical(colnames(res$obj), rownames(md))
})

test_that("similarity metacells recover subclusters; random ones do not", {
  so  <- make_clustered_seurat()
  sim <- quiet(make_metacells_by_size(so, c("cell_type", "sample"), 10,
                                      assays = "RNA", method = "similarity"))
  rnd <- quiet(make_metacells_by_size(so, c("cell_type", "sample"), 10, assays = "RNA"))
  expect_equal(purity(sim$membership), 1)
  expect_lt(purity(rnd$membership), 0.8)
  expect_true("compactness" %in% colnames(sim$obj@meta.data))
  expect_false("compactness" %in% colnames(rnd$obj@meta.data))
})

test_that("similarity and random give the same metacell sizes", {
  so  <- make_clustered_seurat()
  sim <- quiet(make_metacells_by_num(so, c("cell_type", "sample"), 7,
                                     assays = "RNA", method = "similarity"))
  rnd <- quiet(make_metacells_by_num(so, c("cell_type", "sample"), 7, assays = "RNA"))
  expect_equal(sort(unname(sim$obj$n_cells)), sort(unname(rnd$obj$n_cells)))
})

test_that("missing reduction gives a clear error", {
  expect_error(
    quiet(make_metacells_by_size(make_clustered_seurat(), "sample", 10, assays = "RNA",
                                 method = "similarity", reduction = "umap")),
    "not found"
  )
})

# Sum the raw counts per metacell from the membership table.
expected_sums <- function(so, membership) {
  m  <- membership[membership$metacell_id != "unassigned", ]
  cm <- as.matrix(SeuratObject::LayerData(so, assay = "RNA", layer = "counts"))
  t(rowsum(t(cm[, m$cell_barcode]), m$metacell_id))
}

test_that("summed counts are exact and nothing is scaled", {
  so  <- make_clustered_seurat()
  res <- quiet(make_metacells_by_size(so, c("cell_type", "sample"), 10,
                                      assays = "RNA", method = "similarity"))
  got <- as.matrix(SeuratObject::LayerData(res$obj, assay = "RNA", layer = "counts"))
  exp <- expected_sums(so, res$membership)
  expect_equal(got[, colnames(exp)], exp)
  expect_identical(colnames(res$obj), unname(res$obj$metacell_id))
  expect_false("scale.data" %in% SeuratObject::Layers(res$obj[["RNA"]]))
  expect_true("data" %in% SeuratObject::Layers(res$obj[["RNA"]]))

  raw <- quiet(make_metacells_by_size(so, "sample", 10, assays = "RNA",
                                      normalize = FALSE))
  expect_false("data" %in% SeuratObject::Layers(raw$obj[["RNA"]]))
})

test_that("split Seurat v5 counts layers are summed across layers", {
  so  <- make_clustered_seurat()
  ref <- quiet(make_pseudobulk(so, "sample", assays = "RNA"))
  so[["RNA"]] <- split(so[["RNA"]], f = so$subcluster)
  res <- quiet(make_pseudobulk(so, "sample", assays = "RNA"))
  expect_equal(
    as.matrix(SeuratObject::LayerData(res$obj, assay = "RNA", layer = "counts")),
    as.matrix(SeuratObject::LayerData(ref$obj, assay = "RNA", layer = "counts")))
})

test_that("similarity metacells never mix samples or cell types", {
  so <- make_clustered_seurat()
  # The same subclusters appear in both samples, so mixing would be tempting.
  so$cell_type <- ifelse(so$subcluster %in% c("1", "2"), "A", "B")
  for (fn in list(
    function() make_metacells_by_size(so, c("cell_type", "sample"), 10,
                                      assays = "RNA", method = "similarity"),
    function() make_metacells_by_num(so, c("cell_type", "sample"), 4,
                                     assays = "RNA", method = "similarity"))) {
    m <- quiet(fn())$membership
    m <- m[m$metacell_id != "unassigned", ]
    expect_true(all(tapply(m$sample,    m$metacell_id, function(x) length(unique(x))) == 1L))
    expect_true(all(tapply(m$cell_type, m$metacell_id, function(x) length(unique(x))) == 1L))
  }
  sample_only <- quiet(make_metacells_by_size(so, "sample", 10, assays = "RNA",
                                              method = "similarity"))$membership
  expect_true(all(tapply(sample_only$sample, sample_only$metacell_id,
                         function(x) length(unique(x))) == 1L))
})

test_that("group values whose joined IDs collide are rejected", {
  so <- make_clustered_seurat()
  so$cell_type <- rep(c("A_1", "A"), length.out = ncol(so))
  so$sample    <- rep(c("x", "1_x"), length.out = ncol(so))
  expect_error(quiet(make_metacells_by_size(so, c("cell_type", "sample"), 10,
                                            assays = "RNA")),
               "different")
  expect_error(quiet(make_pseudobulk(so, c("cell_type", "sample"), assays = "RNA")),
               "different")
})

test_that("from_log recovers the raw counts exactly", {
  so  <- Seurat::NormalizeData(make_clustered_seurat(), verbose = FALSE)
  ref <- quiet(make_metacells_by_size(so, "sample", 10, assays = "RNA"))
  ref_counts <- SeuratObject::LayerData(ref$obj, assay = "RNA", layer = "counts")
  got <- quiet(make_metacells_by_size(so, "sample", 10, assays = "RNA", from_log = TRUE))
  expect_equal(SeuratObject::LayerData(got$obj, assay = "RNA", layer = "counts"), ref_counts)

  # explicit scale factor; split v5 layers
  so2 <- so
  so2[["RNA"]] <- split(so2[["RNA"]], f = so2$subcluster)
  got2 <- quiet(make_metacells_by_size(so2, "sample", 10, assays = "RNA",
                                       from_log = TRUE, scale_factor = 1e4))
  expect_equal(SeuratObject::LayerData(got2$obj, assay = "RNA", layer = "counts"),
               ref_counts)
})

test_that("from_log falls back to the smallest-non-zero trick with a message", {
  so <- Seurat::NormalizeData(make_clustered_seurat(), verbose = FALSE)
  so$nCount_RNA <- NULL
  expect_message(
    res <- suppressWarnings(make_pseudobulk(so, "sample", assays = "RNA", from_log = TRUE)),
    "smallest-non-zero")
  ref <- quiet(make_pseudobulk(so, "sample", assays = "RNA"))
  # Poisson(5) counts: every cell has at least one count of 1
  expect_equal(SeuratObject::LayerData(res$obj, assay = "RNA", layer = "counts"),
               SeuratObject::LayerData(ref$obj, assay = "RNA", layer = "counts"))
  expect_error(quiet(make_pseudobulk(so, "sample", assays = "RNA", from_log = TRUE,
                                     total_col = "missing")), "not found")
})

test_that("mean_normalized averages in linear space and needs from_log", {
  so  <- Seurat::NormalizeData(make_clustered_seurat(), verbose = FALSE)
  res <- quiet(make_pseudobulk(so, "sample", assays = "RNA",
                               from_log = TRUE, mean_normalized = TRUE))
  x   <- as.matrix(SeuratObject::LayerData(so, assay = "RNA", layer = "data"))
  got <- as.matrix(SeuratObject::LayerData(res$obj, assay = "RNA", layer = "data"))
  for (s in c("s1", "s2")) {
    expect_equal(unname(got[, s]), unname(log1p(rowMeans(expm1(x[, so$sample == s])))))
  }
  expect_false("counts" %in% SeuratObject::Layers(res$obj[["RNA"]]))
  expect_error(make_pseudobulk(so, "sample", assays = "RNA", mean_normalized = TRUE),
               "from_log")
})
