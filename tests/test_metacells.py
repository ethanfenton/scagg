"""Tests for metacell creation strategies."""

import numpy as np
import pytest
from scagg import make_metacells_by_size, make_metacells_by_num


# ── by_size ──────────────────────────────────────────────────────────────────

def test_by_size_n_obs_roughly_correct(small_adata):
    # 200 cells / 10 per metacell = ~20 metacells
    result = make_metacells_by_size(
        small_adata, group_vars=["cell_type", "sample"], cell_size=10
    )
    assert 10 <= result.n_obs <= 30


def test_by_size_counts_non_negative(small_adata):
    result = make_metacells_by_size(
        small_adata, group_vars=["sample"], cell_size=10
    )
    assert (result.X >= 0).all()


def test_by_size_reproducible(small_adata):
    r1 = make_metacells_by_size(
        small_adata, group_vars=["sample"], cell_size=10, seed=0
    )
    r2 = make_metacells_by_size(
        small_adata, group_vars=["sample"], cell_size=10, seed=0
    )
    assert np.array_equal(r1.X, r2.X)


def test_by_size_different_seeds_differ(small_adata):
    r1 = make_metacells_by_size(
        small_adata, group_vars=["sample"], cell_size=5, seed=0
    )
    r2 = make_metacells_by_size(
        small_adata, group_vars=["sample"], cell_size=5, seed=99
    )
    # Different seeds should produce different metacell contents (almost surely)
    assert not np.array_equal(r1.X, r2.X)


def test_by_size_membership(small_adata):
    result = make_metacells_by_size(
        small_adata, group_vars=["cell_type", "sample"], cell_size=10
    )
    m = result.uns["metacell_membership"]
    assert "metacell_id" in m.columns
    assert len(m) == small_adata.n_obs


def test_by_size_n_cells_in_metadata(small_adata):
    result = make_metacells_by_size(
        small_adata, group_vars=["cell_type", "sample"], cell_size=10
    )
    assert "n_cells" in result.obs.columns
    assert (result.obs["n_cells"] > 0).all()


# ── by_num ───────────────────────────────────────────────────────────────────

def test_by_num_exact_count(small_adata):
    n = 5
    result = make_metacells_by_num(
        small_adata, group_vars=["cell_type", "sample"], n_metacells=n
    )
    # 2 cell_types × 4 samples × 5 metacells = 40 (minus any excluded groups)
    assert result.n_obs <= 2 * 4 * n


def test_by_num_reproducible(small_adata):
    r1 = make_metacells_by_num(
        small_adata, group_vars=["sample"], n_metacells=4, seed=1
    )
    r2 = make_metacells_by_num(
        small_adata, group_vars=["sample"], n_metacells=4, seed=1
    )
    assert np.array_equal(r1.X, r2.X)


def test_by_num_min_cells_excludes(small_adata):
    # Setting min_cells very high should produce fewer metacells
    result_low  = make_metacells_by_num(
        small_adata, group_vars=["cell_type", "sample"],
        n_metacells=3, min_cells=1
    )
    result_high = make_metacells_by_num(
        small_adata, group_vars=["cell_type", "sample"],
        n_metacells=3, min_cells=999
    )
    assert result_high.n_obs <= result_low.n_obs


def test_by_num_uns_strategy(small_adata):
    result = make_metacells_by_num(
        small_adata, group_vars=["sample"], n_metacells=3
    )
    assert result.uns["scagg_strategy"].startswith("by_num")
