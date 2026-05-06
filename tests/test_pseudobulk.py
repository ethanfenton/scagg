"""Tests for make_pseudobulk."""

import numpy as np
from scagg import make_pseudobulk


def test_pseudobulk_n_obs(small_adata):
    # 2 cell types × 4 samples = 8 pseudobulk samples
    result = make_pseudobulk(small_adata, group_vars=["cell_type", "sample"])
    assert result.n_obs == 8


def test_pseudobulk_sum_preserved(small_adata):
    # Total counts must be conserved
    result = make_pseudobulk(small_adata, group_vars=["sample"])
    assert np.isclose(result.X.sum(), small_adata.X.sum())


def test_pseudobulk_membership_shape(small_adata):
    result = make_pseudobulk(small_adata, group_vars=["sample"])
    m = result.uns["metacell_membership"]
    assert len(m) == small_adata.n_obs


def test_pseudobulk_meta_vars(small_adata):
    result = make_pseudobulk(
        small_adata,
        group_vars=["sample"],
        meta_vars=["cell_type", "batch"],
    )
    assert "cell_type" in result.obs.columns
    assert "batch" in result.obs.columns


def test_pseudobulk_save_membership(small_adata, tmp_path):
    path = tmp_path / "membership.csv"
    make_pseudobulk(
        small_adata,
        group_vars=["sample"],
        save_membership=str(path),
    )
    assert path.exists()
    import pandas as pd
    df = pd.read_csv(path, index_col=0)
    assert len(df) == small_adata.n_obs
