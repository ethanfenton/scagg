"""Tests for similarity-based (method="similarity") metacells."""

import numpy as np
import pytest
from scagg import make_metacells_by_size, make_metacells_by_num
from scagg._partition import balanced_partition


# ── balanced_partition ───────────────────────────────────────────────────────

@pytest.mark.parametrize("n,k", [(53, 5), (100, 10), (7, 7), (31, 1), (1000, 37)])
@pytest.mark.parametrize("refine", [False, True])
def test_partition_sizes_balanced(n, k, refine):
    emb = np.random.default_rng(n).normal(size=(n, 8))
    sizes = np.bincount(balanced_partition(emb, k, refine=refine), minlength=k)
    assert len(sizes) == k
    assert sizes.sum() == n
    assert set(sizes) <= {n // k, -(-n // k)}


def test_partition_deterministic():
    emb = np.random.default_rng(0).normal(size=(200, 5))
    assert np.array_equal(balanced_partition(emb, 9), balanced_partition(emb, 9))


def test_partition_rejects_bad_k():
    with pytest.raises(ValueError):
        balanced_partition(np.zeros((5, 2)), 6)


def test_partition_identical_cells():
    # Degenerate embedding (no variance) must still give balanced groups
    sizes = np.bincount(balanced_partition(np.zeros((20, 3)), 4))
    assert list(sizes) == [5, 5, 5, 5]


def test_refine_does_not_increase_spread():
    emb = np.random.default_rng(3).normal(size=(500, 10))

    def sse(lab):
        return sum(((emb[lab == j] - emb[lab == j].mean(0)) ** 2).sum()
                   for j in np.unique(lab))

    assert sse(balanced_partition(emb, 25, refine=True)) <= \
        sse(balanced_partition(emb, 25, refine=False)) + 1e-9


# ── public API ───────────────────────────────────────────────────────────────

def _purity(result):
    m = result.uns["metacell_membership"]
    m = m[m["metacell_id"] != "unassigned"]
    return m.groupby("metacell_id")["subcluster"].agg(
        lambda s: s.value_counts().iloc[0] / len(s)).mean()


def test_similarity_recovers_subclusters(clustered_adata):
    # 30 cells per subcluster, 10 per metacell → perfect purity is achievable
    sim = make_metacells_by_size(clustered_adata, ["cell_type", "sample"], 10,
                                 method="similarity")
    rnd = make_metacells_by_size(clustered_adata, ["cell_type", "sample"], 10)
    assert _purity(sim) == 1.0
    assert _purity(rnd) < 0.8


def test_similarity_same_sizes_as_random(clustered_adata):
    kw = dict(group_vars=["cell_type", "sample"], n_metacells=7)
    sim = make_metacells_by_num(clustered_adata, method="similarity", **kw)
    rnd = make_metacells_by_num(clustered_adata, **kw)
    assert sim.n_obs == rnd.n_obs
    assert sorted(sim.obs["n_cells"]) == sorted(rnd.obs["n_cells"])


def test_similarity_compactness_column(clustered_adata):
    sim = make_metacells_by_size(clustered_adata, ["sample"], 10, method="similarity")
    assert "compactness" in sim.obs.columns
    assert sim.obs["compactness"].notna().all()
    rnd = make_metacells_by_size(clustered_adata, ["sample"], 10)
    assert "compactness" not in rnd.obs.columns


def test_similarity_respects_strata(clustered_adata):
    sim = make_metacells_by_size(clustered_adata, ["sample"], 10, method="similarity")
    m = sim.uns["metacell_membership"]
    assert (m.groupby("metacell_id")["sample"].nunique() == 1).all()


def test_similarity_counts_preserved(clustered_adata):
    sim = make_metacells_by_num(clustered_adata, ["sample"], 6, method="similarity")
    assert np.isclose(sim.X.sum(), clustered_adata.X.sum())


def test_similarity_ignores_seed(clustered_adata):
    r1 = make_metacells_by_size(clustered_adata, ["sample"], 10, method="similarity", seed=0)
    r2 = make_metacells_by_size(clustered_adata, ["sample"], 10, method="similarity", seed=9)
    assert np.array_equal(r1.X, r2.X)


def test_similarity_n_dims(clustered_adata):
    r = make_metacells_by_size(clustered_adata, ["sample"], 10,
                               method="similarity", n_dims=2)
    assert r.uns["scagg_strategy"] == "by_size_10_similarity"


def test_missing_rep_raises(clustered_adata):
    with pytest.raises(KeyError, match="X_umap"):
        make_metacells_by_size(clustered_adata, ["sample"], 10,
                               method="similarity", use_rep="X_umap")


def test_bad_method_raises(clustered_adata):
    with pytest.raises(ValueError, match="method"):
        make_metacells_by_size(clustered_adata, ["sample"], 10, method="knn")


# ── n_init_pcs ───────────────────────────────────────────────────────────────

def _gradient_plus_subtypes(seed=0, n=400):
    """Broad gradient on x (var ≈ 10) dominates PC1; two subtypes split on y (±3)."""
    rng = np.random.default_rng(seed)
    sub = np.repeat([0, 1], n // 2)
    emb = np.column_stack([
        rng.uniform(-5.48, 5.48, n),
        np.where(sub, 3.0, -3.0) + rng.normal(0, 0.3, n),
        rng.normal(0, 0.3, (n, 8)),
    ])
    return emb, sub


def _subtype_purity(lab, sub):
    return np.mean([np.bincount(sub[lab == j], minlength=2).max() / (lab == j).sum()
                    for j in np.unique(lab)])


def test_n_init_pcs_finds_split_on_later_pc():
    emb, sub = _gradient_plus_subtypes()
    assert _subtype_purity(balanced_partition(emb, 20, n_init_pcs=3), sub) == 1.0
    assert _subtype_purity(balanced_partition(emb, 20, n_init_pcs=1), sub) < 1.0


def test_n_init_pcs_default_is_pc1_only():
    emb = np.random.default_rng(4).normal(size=(300, 10))
    assert np.array_equal(balanced_partition(emb, 17),
                          balanced_partition(emb, 17, n_init_pcs=1))


@pytest.mark.parametrize("n_init_pcs", [2, 3, 50])
def test_n_init_pcs_keeps_balance(n_init_pcs):
    emb = np.random.default_rng(5).normal(size=(211, 6))   # 50 > n_dims is allowed
    sizes = np.bincount(balanced_partition(emb, 13, n_init_pcs=n_init_pcs))
    assert set(sizes) <= {211 // 13, -(-211 // 13)}


def test_n_init_pcs_validated(clustered_adata):
    with pytest.raises(ValueError, match="n_init_pcs"):
        make_metacells_by_size(clustered_adata, ["sample"], 10,
                               method="similarity", n_init_pcs=0)
