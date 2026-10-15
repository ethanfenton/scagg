"""Similarity-aware, size-balanced partitioning of cells within a group.

The partition is built in two stages:

1. **Balanced recursive bisection.**  Cells are projected onto their first
   principal component, sorted, and cut at the rank that gives each side a
   share of cells proportional to the number of metacells it will hold.  The
   cut is then polished with size-constrained 2-means (re-sort cells by
   ``d²(left centroid) − d²(right centroid)`` and re-cut at the same rank).
   Each side is split again along *its own* first PC, until every piece is
   one metacell.  Cutting at ``floor(n * k_left / k)`` guarantees every final
   metacell has ``floor(n / k)`` or ``ceil(n / k)`` cells.

2. **Pairwise refinement (optional).**  For neighbouring metacells *a* and *b*
   the pooled cells are re-split, keeping both sizes fixed, by sorting on
   ``d²(cell, centroid_a) − d²(cell, centroid_b)``.  This is the optimal
   size-constrained two-way assignment for fixed centroids, so the total
   within-metacell sum of squares never increases.
"""

from __future__ import annotations

import numpy as np


def balanced_partition(
    emb: np.ndarray,
    k: int,
    *,
    refine: bool = True,
    n_init_pcs: int = 1,
    max_iter: int = 10,
    n_neighbors: int = 3,
) -> np.ndarray:
    """Split the rows of *emb* into *k* groups of near-equal size.

    *n_init_pcs* > 1 tries each split starting from PC1 … PC{n_init_pcs} and
    keeps the cut with the lowest within-half sum of squares.

    Returns integer labels ``0 .. k-1``; group sizes differ by at most one.
    """
    n = emb.shape[0]
    if not 1 <= k <= n:
        raise ValueError(f"k must be between 1 and n ({n}), got {k}")
    if n_init_pcs < 1:
        raise ValueError(f"n_init_pcs must be >= 1, got {n_init_pcs}")

    labels = np.empty(n, dtype=np.int64)
    next_label = 0

    # Explicit stack (left branch processed first) → labels follow tree order.
    stack = [(np.arange(n), k)]
    while stack:
        idx, kk = stack.pop()
        if kk == 1:
            labels[idx] = next_label
            next_label += 1
            continue
        k_left = kk // 2
        n_left = len(idx) * k_left // kk
        order = _balanced_split_order(emb[idx], n_left, n_init_pcs)
        stack.append((idx[order[n_left:]], kk - k_left))
        stack.append((idx[order[:n_left]], k_left))

    if refine and k > 1:
        labels = _refine_pairs(emb, labels, k, max_iter=max_iter, n_neighbors=n_neighbors)
    return labels


def compactness(emb: np.ndarray, labels: np.ndarray) -> dict:
    """Mean Euclidean distance from each cell to its group centroid, per label."""
    out = {}
    for lab in np.unique(labels):
        block = emb[labels == lab]
        out[lab] = float(np.linalg.norm(block - block.mean(axis=0), axis=1).mean())
    return out


# ── Internal ─────────────────────────────────────────────────────────────────

def _pc_orders(emb: np.ndarray, n_pcs: int) -> list[np.ndarray]:
    """Row orderings of *emb* by score on each of its first *n_pcs* PCs."""
    if emb.shape[0] < 2:
        return [np.arange(emb.shape[0])]
    centred = emb - emb.mean(axis=0)
    _, _, vt = np.linalg.svd(centred, full_matrices=False)
    orders = []
    for v in vt[:n_pcs]:
        v = v * np.sign(v[np.argmax(np.abs(v))])   # fix SVD sign for determinism
        orders.append(np.argsort(centred @ v, kind="stable"))
    return orders


def _balanced_split_order(
    emb: np.ndarray,
    n_left: int,
    n_init_pcs: int = 1,
    max_iter: int = 10,
) -> np.ndarray:
    """Order rows so that ``order[:n_left]`` / ``order[n_left:]`` is a good split.

    Starts from a PC ordering, then runs size-constrained 2-means: compute
    both halves' centroids and re-sort by ``d²(left) − d²(right)``.  This lets
    the cut follow gaps in the full embedding rather than only along one PC.
    With *n_init_pcs* > 1 this is repeated from each of the first PCs and the
    split with the lowest within-half sum of squares wins (ties → earlier PC).
    """
    orders = _pc_orders(emb, n_init_pcs)
    if n_left == 0 or n_left == emb.shape[0]:
        return orders[0]
    best_order, best_cost = None, None
    for start in orders:
        order, cost = _polish_split(emb, start, n_left, max_iter)
        if best_cost is None or cost < best_cost - 1e-9 * (best_cost + 1.0):
            best_order, best_cost = order, cost
    return best_order


def _polish_split(
    emb: np.ndarray,
    order: np.ndarray,
    n_left: int,
    max_iter: int,
) -> tuple[np.ndarray, float]:
    """Size-constrained 2-means from *order*; returns (order, within-half SS)."""
    left = np.zeros(len(order), dtype=bool)
    left[order[:n_left]] = True
    for _ in range(max_iter):
        ca, cb = emb[left].mean(axis=0), emb[~left].mean(axis=0)
        diff = ((emb - ca) ** 2).sum(axis=1) - ((emb - cb) ** 2).sum(axis=1)
        new_order = np.argsort(diff, kind="stable")
        new_left = np.zeros_like(left)
        new_left[new_order[:n_left]] = True
        order = new_order
        if np.array_equal(new_left, left):
            break
        left = new_left
    left = np.zeros(len(order), dtype=bool)
    left[order[:n_left]] = True
    cost = (((emb[left] - emb[left].mean(axis=0)) ** 2).sum()
            + ((emb[~left] - emb[~left].mean(axis=0)) ** 2).sum())
    return order, float(cost)


def _refine_pairs(
    emb: np.ndarray,
    labels: np.ndarray,
    k: int,
    *,
    max_iter: int,
    n_neighbors: int,
) -> np.ndarray:
    m = min(n_neighbors, k - 1)
    members = [np.flatnonzero(labels == j) for j in range(k)]
    for _ in range(max_iter):
        centroids = np.vstack([emb[idx].mean(axis=0) for idx in members])
        sq = (centroids ** 2).sum(axis=1)
        cdist = sq[:, None] + sq[None, :] - 2.0 * centroids @ centroids.T
        np.fill_diagonal(cdist, np.inf)
        nn = np.argpartition(cdist, m - 1, axis=1)[:, :m] if m < k - 1 else \
            np.argsort(cdist, axis=1)[:, :m]
        pairs = sorted({(min(a, b), max(a, b)) for a in range(k) for b in nn[a]})

        improved = False
        for a, b in pairs:
            idx = np.concatenate([members[a], members[b]])
            n_a = len(members[a])
            sub = emb[idx]
            ca, cb = sub[:n_a].mean(axis=0), sub[n_a:].mean(axis=0)
            da = ((sub - ca) ** 2).sum(axis=1)
            db = ((sub - cb) ** 2).sum(axis=1)

            order = np.argsort(da - db, kind="stable")
            old_cost = da[:n_a].sum() + db[n_a:].sum()
            new_cost = da[order[:n_a]].sum() + db[order[n_a:]].sum()
            if new_cost < old_cost - 1e-9 * (old_cost + 1.0):
                members[a] = np.sort(idx[order[:n_a]])
                members[b] = np.sort(idx[order[n_a:]])
                improved = True
        if not improved:
            break

    out = np.empty_like(labels)
    for j, idx in enumerate(members):
        out[idx] = j
    return out
