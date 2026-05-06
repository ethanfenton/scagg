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
