"""Aggregating log-normalised input (from_log / mean_normalized)."""

import numpy as np
import pandas as pd
import pytest
import anndata as ad
from scipy.sparse import csr_matrix

from scagg import make_metacells_by_size, make_pseudobulk


def _lognorm_adata(target_sum, sparse=True, with_totals=True):
    rng = np.random.default_rng(3)
    counts = rng.poisson(0.8, size=(120, 40)).astype(float)
    counts[:, 0] += 1                      # every cell has some counts
    totals = counts.sum(axis=1)
    if target_sum is None:
        target_sum = np.median(totals)     # scanpy normalize_total default
    lognorm = np.log1p(counts / totals[:, None] * target_sum)
    obs = pd.DataFrame({"sample": np.repeat(["s1", "s2", "s3"], 40)},
                       index=[f"c{i}" for i in range(120)])
    if with_totals:
        obs["total_counts"] = totals
    a = ad.AnnData(X=csr_matrix(lognorm) if sparse else lognorm, obs=obs)
    a.layers["counts"] = counts
    return a


@pytest.mark.parametrize("target_sum", [1e4, None])
@pytest.mark.parametrize("sparse", [True, False])
def test_from_log_recovers_exact_counts(target_sum, sparse):
    a = _lognorm_adata(target_sum, sparse)
    got = make_metacells_by_size(a, "sample", 10, from_log=True)
    ref = make_metacells_by_size(a, "sample", 10, layer="counts")
    assert np.array_equal(got.X, ref.X)
    assert got.uns["scagg_values"] == "recovered_counts"


def test_from_log_fallback_message(capsys):
    a = _lognorm_adata(1e4, with_totals=False)
    got = make_pseudobulk(a, "sample", from_log=True)
    assert "smallest-non-zero" in capsys.readouterr().out
    ref = make_pseudobulk(a, "sample", layer="counts")
    # exact whenever each cell's smallest non-zero count is 1
    assert np.array_equal(got.X, ref.X)


def test_from_log_explicit_scale_factor_and_total_col():
    a = _lognorm_adata(1e4)
    a.obs = a.obs.rename(columns={"total_counts": "my_total"})
    got = make_pseudobulk(a, "sample", from_log=True, total_col="my_total",
                          scale_factor=1e4)
    ref = make_pseudobulk(a, "sample", layer="counts")
    assert np.array_equal(got.X, ref.X)
    with pytest.raises(KeyError):
        make_pseudobulk(a, "sample", from_log=True, total_col="missing")


def test_from_log_warns_on_wrong_scale_factor():
    a = _lognorm_adata(1e4)
    with pytest.warns(UserWarning, match="whole numbers"):
        make_pseudobulk(a, "sample", from_log=True, scale_factor=3e3)


def test_mean_normalized():
    a = _lognorm_adata(1e4)
    got = make_pseudobulk(a, "sample", from_log=True, mean_normalized=True)
    x = a.X.toarray()
    for s in ["s1", "s2", "s3"]:
        exp = np.log1p(np.expm1(x[a.obs["sample"] == s]).mean(axis=0))
        assert np.allclose(got[s].X.ravel(), exp)
    assert got.uns["scagg_values"] == "mean_normalized"


def test_mean_normalized_requires_from_log():
    a = _lognorm_adata(1e4)
    with pytest.raises(ValueError, match="from_log"):
        make_pseudobulk(a, "sample", mean_normalized=True)


def test_summing_log_data_warns():
    a = _lognorm_adata(1e4)
    with pytest.warns(UserWarning, match="from_log"):
        make_pseudobulk(a, "sample")


def test_from_log_cells_without_a_count_of_one():
    rng = np.random.default_rng(5)
    counts = rng.poisson(6, size=(60, 15)).astype(float) * 2 + 2   # even, min ≥ 2
    counts[:, 0] += 1                                             # one odd count
    totals = counts.sum(axis=1)
    obs = pd.DataFrame({"sample": np.repeat(["s1", "s2"], 30)},
                       index=[f"c{i}" for i in range(60)])
    a = ad.AnnData(X=np.log1p(counts / totals[:, None] * 1e4), obs=obs)
    a.layers["counts"] = counts
    ref = make_pseudobulk(a, "sample", layer="counts")
    assert np.array_equal(make_pseudobulk(a, "sample", from_log=True).X, ref.X)
    a.obs["total_counts"] = totals
    assert np.array_equal(make_pseudobulk(a, "sample", from_log=True).X, ref.X)
