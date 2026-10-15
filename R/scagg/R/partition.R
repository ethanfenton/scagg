# ── Similarity-aware, size-balanced partitioning ─────────────────────────────
#
# Port of src/scagg/_partition.py; see that file for the algorithm notes.
#
# 1. Balanced recursive bisection: order cells by PC1 within the current
#    piece, polish the cut with size-constrained 2-means, and cut at
#    floor(n * k_left / k).  Every final group has floor(n/k) or ceil(n/k)
#    cells.  With n_init_pcs > 1 each split is tried from PC1..PCn and the
#    tightest cut wins.
# 2. Optional pairwise refinement: re-split the pooled cells of neighbouring
#    groups with both sizes fixed; never increases within-group spread.

# Split the rows of `emb` into `k` groups of near-equal size.
# Returns integer labels 1..k.
.balanced_partition <- function(emb, k, refine = TRUE, n_init_pcs = 1L,
                                max_iter = 10L, n_neighbors = 3L) {
  emb <- as.matrix(emb)
  n <- nrow(emb)
  if (k < 1L || k > n) stop(sprintf("k must be between 1 and n (%d), got %d", n, k))
  if (n_init_pcs < 1L) stop("n_init_pcs must be >= 1")

  labels <- integer(n)
  next_label <- 1L
  # Explicit stack (left branch processed first) → labels follow tree order.
  stack <- list(list(idx = seq_len(n), k = k))
  while (length(stack) > 0L) {
    top <- stack[[length(stack)]]
    stack[[length(stack)]] <- NULL
    idx <- top$idx; kk <- top$k
    if (kk == 1L) {
      labels[idx] <- next_label
      next_label <- next_label + 1L
      next
    }
    k_left <- kk %/% 2L
    n_left <- (length(idx) * k_left) %/% kk
    ord <- .balanced_split_order(emb[idx, , drop = FALSE], n_left, n_init_pcs)
    left  <- idx[ord[seq_len(n_left)]]
    right <- idx[ord[-seq_len(n_left)]]
    stack[[length(stack) + 1L]] <- list(idx = right, k = kk - k_left)
    stack[[length(stack) + 1L]] <- list(idx = left,  k = k_left)
  }

  if (refine && k > 1L) {
    labels <- .refine_pairs(emb, labels, k, max_iter, n_neighbors)
  }
  labels
}

# Mean Euclidean distance from each cell to its group centroid, per label.
.compactness <- function(emb, labels) {
  vapply(split(seq_len(nrow(emb)), labels), function(i) {
    block <- emb[i, , drop = FALSE]
    mean(sqrt(rowSums(sweep(block, 2L, colMeans(block))^2)))
  }, numeric(1))
}

# Row orderings of `emb` by score on each of its first `n_pcs` PCs.
.pc_orders <- function(emb, n_pcs) {
  if (nrow(emb) < 2L) return(list(seq_len(nrow(emb))))
  centred <- sweep(emb, 2L, colMeans(emb))
  vs <- svd(centred, nu = 0L, nv = min(n_pcs, dim(centred)))$v
  lapply(seq_len(ncol(vs)), function(j) {
    v <- vs[, j]
    v <- v * sign(v[which.max(abs(v))])   # fix SVD sign for determinism
    order(drop(centred %*% v))
  })
}

.sq_dist_to <- function(emb, centre) rowSums(sweep(emb, 2L, centre)^2)

.within_ss <- function(block) sum(sweep(block, 2L, colMeans(block))^2)

# Try each PC start, polish, keep the split with the lowest within-half SS
# (ties → earlier PC).
.balanced_split_order <- function(emb, n_left, n_init_pcs = 1L, max_iter = 10L) {
  orders <- .pc_orders(emb, n_init_pcs)
  if (n_left == 0L || n_left == nrow(emb)) return(orders[[1L]])
  best <- NULL; best_cost <- NULL
  for (start in orders) {
    ord  <- .polish_split(emb, start, n_left, max_iter)
    left <- ord[seq_len(n_left)]
    cost <- .within_ss(emb[left, , drop = FALSE]) + .within_ss(emb[-left, , drop = FALSE])
    if (is.null(best_cost) || cost < best_cost - 1e-9 * (best_cost + 1)) {
      best <- ord; best_cost <- cost
    }
  }
  best
}

# Size-constrained 2-means starting from ordering `ord`.
.polish_split <- function(emb, ord, n_left, max_iter) {
  left <- logical(length(ord))
  left[ord[seq_len(n_left)]] <- TRUE
  for (it in seq_len(max_iter)) {
    diff <- .sq_dist_to(emb, colMeans(emb[left, , drop = FALSE])) -
            .sq_dist_to(emb, colMeans(emb[!left, , drop = FALSE]))
    ord <- order(diff)
    new_left <- logical(length(ord))
    new_left[ord[seq_len(n_left)]] <- TRUE
    if (identical(new_left, left)) break
    left <- new_left
  }
  ord
}

.refine_pairs <- function(emb, labels, k, max_iter, n_neighbors) {
  m <- min(n_neighbors, k - 1L)
  members <- split(seq_along(labels), factor(labels, levels = seq_len(k)))
  for (it in seq_len(max_iter)) {
    centroids <- t(vapply(members, function(i) colMeans(emb[i, , drop = FALSE]),
                          numeric(ncol(emb))))
    if (ncol(emb) == 1L) centroids <- t(centroids)
    cdist <- as.matrix(stats::dist(centroids))^2
    diag(cdist) <- Inf
    nn <- t(apply(cdist, 1L, function(r) order(r)[seq_len(m)]))
    if (m == 1L) nn <- matrix(nn, ncol = 1L)
    pairs <- unique(t(apply(cbind(rep(seq_len(k), m), as.vector(nn)), 1L, sort)))
    pairs <- pairs[order(pairs[, 1L], pairs[, 2L]), , drop = FALSE]  # match Python

    improved <- FALSE
    for (p in seq_len(nrow(pairs))) {
      a <- pairs[p, 1L]; b <- pairs[p, 2L]
      idx <- c(members[[a]], members[[b]])
      n_a <- length(members[[a]])
      sub <- emb[idx, , drop = FALSE]
      in_a <- seq_along(idx) <= n_a
      da <- .sq_dist_to(sub, colMeans(sub[in_a, , drop = FALSE]))
      db <- .sq_dist_to(sub, colMeans(sub[!in_a, , drop = FALSE]))
      ord <- order(da - db)
      old_cost <- sum(da[in_a]) + sum(db[!in_a])
      new_cost <- sum(da[ord[seq_len(n_a)]]) + sum(db[ord[-seq_len(n_a)]])
      if (new_cost < old_cost - 1e-9 * (old_cost + 1)) {
        members[[a]] <- sort(idx[ord[seq_len(n_a)]])
        members[[b]] <- sort(idx[ord[-seq_len(n_a)]])
        improved <- TRUE
      }
    }
    if (!improved) break
  }
  out <- integer(length(labels))
  for (j in seq_len(k)) out[members[[j]]] <- j
  out
}
