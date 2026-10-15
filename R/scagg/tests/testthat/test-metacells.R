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
