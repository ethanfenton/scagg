test_that("partition sizes are balanced to within one cell", {
  for (case in list(c(53, 5), c(100, 10), c(7, 7), c(31, 1), c(1000, 37))) {
    n <- case[1]; k <- case[2]
    set.seed(n)
    emb <- matrix(stats::rnorm(n * 8), n)
    for (refine in c(FALSE, TRUE)) {
      sizes <- tabulate(scagg:::.balanced_partition(emb, k, refine = refine), k)
      expect_equal(sum(sizes), n)
      expect_true(all(sizes %in% c(n %/% k, ceiling(n / k))))
    }
  }
})

test_that("partition handles degenerate embeddings and bad k", {
  expect_equal(tabulate(scagg:::.balanced_partition(matrix(0, 20, 3), 4)), rep(5L, 4))
  expect_error(scagg:::.balanced_partition(matrix(0, 5, 2), 6))
})

test_that("refinement never increases within-group spread", {
  set.seed(3)
  emb <- matrix(stats::rnorm(500 * 10), 500)
  sse <- function(lab) sum(vapply(split(seq_len(nrow(emb)), lab), function(i) {
    b <- emb[i, , drop = FALSE]; sum(sweep(b, 2, colMeans(b))^2)
  }, numeric(1)))
  expect_lte(sse(scagg:::.balanced_partition(emb, 25, refine = TRUE)),
             sse(scagg:::.balanced_partition(emb, 25, refine = FALSE)) + 1e-9)
})

test_that("n_init_pcs finds a split along a later PC", {
  set.seed(0)
  n <- 400
  sub <- rep(0:1, each = n / 2)
  emb <- cbind(stats::runif(n, -5.48, 5.48),
               ifelse(sub == 1, 3, -3) + stats::rnorm(n, 0, 0.3),
               matrix(stats::rnorm(n * 8, 0, 0.3), n))
  pur <- function(lab) mean(tapply(sub, lab, function(s) max(table(s)) / length(s)))
  expect_equal(pur(scagg:::.balanced_partition(emb, 20, n_init_pcs = 3)), 1)
  expect_lt(pur(scagg:::.balanced_partition(emb, 20, n_init_pcs = 1)), 1)
  sizes <- tabulate(scagg:::.balanced_partition(emb, 20, n_init_pcs = 50))
  expect_true(all(sizes == 20))
})
