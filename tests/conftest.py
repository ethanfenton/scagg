"""Shared fixtures for scagg tests."""

import numpy as np
import pandas as pd
import pytest
import anndata as ad


@pytest.fixture
def small_adata():
    """200 cells, 50 genes, 2 cell types, 4 samples."""
    rng = np.random.default_rng(0)
    n_cells, n_genes = 200, 50
    X = rng.poisson(5, size=(n_cells, n_genes)).astype(float)
    obs = pd.DataFrame({
        "cell_type": np.tile(["A", "B"], n_cells // 2),
        "sample":    np.repeat(["s1", "s2", "s3", "s4"], n_cells // 4),
        "batch":     np.tile(["b1", "b2"], n_cells // 2),
        "n_counts":  rng.integers(500, 5000, n_cells).astype(float),
    })
    var = pd.DataFrame(index=[f"gene{i}" for i in range(n_genes)])
    return ad.AnnData(X=X, obs=obs, var=var)


@pytest.fixture
def clustered_adata():
    """2 samples × 4 well-separated subclusters × 30 cells, with an X_pca embedding."""
    rng = np.random.default_rng(1)
    n_sub, per_sub, n_dims = 4, 30, 10
    centers = rng.normal(0, 10, size=(n_sub, n_dims))
    subcluster = np.tile(np.repeat(np.arange(n_sub), per_sub), 2)
    n_cells = len(subcluster)
    emb = centers[subcluster] + rng.normal(0, 1, size=(n_cells, n_dims))
    # Shuffle so subclusters are not contiguous in the input order
    perm = rng.permutation(n_cells)
    obs = pd.DataFrame({
        "cell_type":  "T",
        "sample":     np.repeat(["s1", "s2"], n_sub * per_sub)[perm],
        "subcluster": subcluster[perm].astype(str),
    }, index=[f"cell{i}" for i in range(n_cells)])
    a = ad.AnnData(
        X=rng.poisson(5, size=(n_cells, 20)).astype(float),
        obs=obs,
        var=pd.DataFrame(index=[f"gene{i}" for i in range(20)]),
    )
    a.obsm["X_pca"] = emb[perm]
    return a
