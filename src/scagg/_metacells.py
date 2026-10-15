"""
AnnData-based metacell and pseudobulk creation.

All three public functions return an AnnData where:
  - obs     = one row per metacell / pseudobulk sample
  - X       = summed counts
  - uns['metacell_membership'] = DataFrame mapping every input cell to its group
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from pathlib import Path
from typing import Literal, Sequence

import anndata as ad

from ._aggregate import (
    _assign_by_num,
    _assign_by_size,
    _assign_pseudobulk,
    _aggregate_metadata,
    _default_cell_min,
    _membership_table,
    _sum_by_group,
)
from ._partition import compactness as _compactness

_METHODS = ("random", "similarity")


def make_pseudobulk(
    adata: ad.AnnData,
    group_vars: str | Sequence[str],
    *,
    meta_vars: Sequence[str] | None = None,
    layer: str | None = None,
    save_membership: str | Path | None = None,
) -> ad.AnnData:
    """Aggregate all cells that share the same values of *group_vars* into one
    pseudobulk sample.

    This is true pseudobulk aggregation (one sample per unique group-variable
    combination) rather than metacell creation.

    Parameters
    ----------
    adata:
        Input AnnData (cells × genes).
    group_vars:
        One or more obs column names to group by, e.g. ``["sample", "cell_type"]``.
    meta_vars:
        Additional obs columns to carry forward into the result (mode for
        categorical, mean for numeric). ``None`` = carry all obs columns.
    layer:
        Which layer to sum. ``None`` = use ``adata.X``.
    save_membership:
        If given, save the cell → pseudobulk mapping CSV to this path.

    Returns
    -------
    AnnData with one obs per unique group combination.
    ``adata.uns['metacell_membership']`` contains the full mapping table.
    """
    if isinstance(group_vars, str):
        group_vars = [group_vars]
    group_vars = list(group_vars)

    obs = adata.obs.copy()
    obs = _assign_pseudobulk(obs, group_vars, metacell_col="metacell_id")

    return _build_result(
        adata=adata,
        obs_with_mc=obs,
        group_vars=group_vars,
        meta_vars=meta_vars,
        layer=layer,
        save_membership=save_membership,
        strategy="pseudobulk",
    )


def make_metacells_by_size(
    adata: ad.AnnData,
    group_vars: str | Sequence[str],
    cell_size: int,
    *,
    cell_min: int | None = None,
    meta_vars: Sequence[str] | None = None,
    layer: str | None = None,
    seed: int = 42,
    save_membership: str | Path | None = None,
    method: Literal["random", "similarity"] = "random",
    use_rep: str = "X_pca",
    n_dims: int | None = None,
    refine: bool = True,
    n_init_pcs: int = 1,
) -> ad.AnnData:
    """Aggregate cells into metacells of approximately *cell_size* cells.

    Within each group (defined by *group_vars*) cells are assigned to
    ``floor(n / cell_size)`` metacells — randomly, or by similarity in an
    embedding (see *method*).  Groups with fewer than *cell_min* cells are
    excluded.

    Parameters
    ----------
    adata:
        Input AnnData.
    group_vars:
        Obs column(s) that define independent grouping strata. Metacells are
        created separately within each combination, e.g. ``["cell_type", "sample"]``.
    cell_size:
        Target number of cells per metacell.
    cell_min:
        Minimum cells required to form at least one metacell.
        Default: 70 % of *cell_size* (7 for cell_size=10, 15 for cell_size=25).
    meta_vars:
        Obs columns to carry into metacell metadata.
    layer:
        Layer to sum. ``None`` = ``adata.X``.
    seed:
        Random seed for reproducibility.
    save_membership:
        Path for the cell → metacell CSV.
    method:
        ``"random"`` (default): cells are shuffled and dealt round-robin into
        metacells.  ``"similarity"``: cells that are close in the embedding
        *use_rep* are grouped together via balanced recursive PCA bisection.
        Both methods give identical metacell counts and sizes (within ±1
        cell); only which cells share a metacell differs.  The similarity
        method is deterministic, so *seed* has no effect on it.
    use_rep:
        ``adata.obsm`` key of the embedding used when
        ``method="similarity"``, e.g. ``"X_pca"``, ``"X_harmony"``,
        ``"X_scVI"``.
    n_dims:
        Use only the first *n_dims* columns of *use_rep*.  ``None`` = all.
    refine:
        When ``method="similarity"``, polish the bisection by size-preserving
        swaps between neighbouring metacells (never increases
        within-metacell spread).
    n_init_pcs:
        When ``method="similarity"``, try each split starting from PC1 …
        PC{n_init_pcs} and keep the tightest cut.  ``1`` (default) starts from
        PC1 only; ``3`` helps when a broad gradient dominates PC1 but
        distinct sub-states separate along a later PC.  Cost grows roughly
        linearly with *n_init_pcs*.

    Returns
    -------
    AnnData with one obs per metacell.  With ``method="similarity"``,
    ``obs["compactness"]`` holds each metacell's mean cell-to-centroid
    distance in the embedding.
    """
    if isinstance(group_vars, str):
        group_vars = [group_vars]
    group_vars = list(group_vars)

    if cell_min is None:
        cell_min = _default_cell_min(cell_size)

    emb = _get_embedding(adata, method, use_rep, n_dims)
    if n_init_pcs < 1:
        raise ValueError(f"n_init_pcs must be >= 1, got {n_init_pcs}")
    rng = np.random.default_rng(seed)
    obs = adata.obs.copy()
    obs = _assign_by_size(obs, group_vars, cell_size, cell_min, rng, "metacell_id",
                          emb=emb, refine=refine, n_init_pcs=n_init_pcs)

    return _build_result(
        adata=adata,
        obs_with_mc=obs,
        group_vars=group_vars,
        meta_vars=meta_vars,
        layer=layer,
        save_membership=save_membership,
        strategy=_strategy_name(f"by_size_{cell_size}", method),
        emb=emb,
        seed=seed,
    )


def make_metacells_by_num(
    adata: ad.AnnData,
    group_vars: str | Sequence[str],
    n_metacells: int,
    *,
    min_cells: int | None = None,
    meta_vars: Sequence[str] | None = None,
    layer: str | None = None,
    seed: int = 42,
    save_membership: str | Path | None = None,
    method: Literal["random", "similarity"] = "random",
    use_rep: str = "X_pca",
    n_dims: int | None = None,
    refine: bool = True,
    n_init_pcs: int = 1,
) -> ad.AnnData:
    """Aggregate cells into exactly *n_metacells* metacells per group.

    Cells are assigned randomly or by similarity in an embedding (see
    *method*).

    Parameters
    ----------
    adata:
        Input AnnData.
    group_vars:
        Obs column(s) defining grouping strata.
    n_metacells:
        Number of metacells to create within each group.
    min_cells:
        Minimum cells needed to form *n_metacells* groups.
        Default: equal to *n_metacells*.
    meta_vars:
        Obs columns to carry into metacell metadata.
    layer:
        Layer to sum.
    seed:
        Random seed.
    save_membership:
        Path for the cell → metacell CSV.
    method:
        ``"random"`` (default): cells are shuffled and dealt round-robin into
        metacells.  ``"similarity"``: cells that are close in the embedding
        *use_rep* are grouped together via balanced recursive PCA bisection.
        Both methods give identical metacell counts and sizes (within ±1
        cell); only which cells share a metacell differs.  The similarity
        method is deterministic, so *seed* has no effect on it.
    use_rep:
        ``adata.obsm`` key of the embedding used when
        ``method="similarity"``, e.g. ``"X_pca"``, ``"X_harmony"``,
        ``"X_scVI"``.
    n_dims:
        Use only the first *n_dims* columns of *use_rep*.  ``None`` = all.
    refine:
        When ``method="similarity"``, polish the bisection by size-preserving
        swaps between neighbouring metacells (never increases
        within-metacell spread).
    n_init_pcs:
        When ``method="similarity"``, try each split starting from PC1 …
        PC{n_init_pcs} and keep the tightest cut.  ``1`` (default) starts from
        PC1 only; ``3`` helps when a broad gradient dominates PC1 but
        distinct sub-states separate along a later PC.  Cost grows roughly
        linearly with *n_init_pcs*.

    Returns
    -------
    AnnData with one obs per metacell.  With ``method="similarity"``,
    ``obs["compactness"]`` holds each metacell's mean cell-to-centroid
    distance in the embedding.
    """
    if isinstance(group_vars, str):
        group_vars = [group_vars]
    group_vars = list(group_vars)

    if min_cells is None:
        min_cells = n_metacells

    emb = _get_embedding(adata, method, use_rep, n_dims)
    if n_init_pcs < 1:
        raise ValueError(f"n_init_pcs must be >= 1, got {n_init_pcs}")
    rng = np.random.default_rng(seed)
    obs = adata.obs.copy()
    obs = _assign_by_num(obs, group_vars, n_metacells, min_cells, rng, "metacell_id",
                         emb=emb, refine=refine, n_init_pcs=n_init_pcs)

    return _build_result(
        adata=adata,
        obs_with_mc=obs,
        group_vars=group_vars,
        meta_vars=meta_vars,
        layer=layer,
        save_membership=save_membership,
        strategy=_strategy_name(f"by_num_{n_metacells}", method),
        emb=emb,
        seed=seed,
    )


# ── Internal helpers ─────────────────────────────────────────────────────────

def _strategy_name(base: str, method: str) -> str:
    return base if method == "random" else f"{base}_{method}"


def _get_embedding(
    adata: ad.AnnData,
    method: str,
    use_rep: str,
    n_dims: int | None,
) -> np.ndarray | None:
    """Validate *method* and return the embedding matrix (``None`` for random)."""
    if method not in _METHODS:
        raise ValueError(f"method must be one of {_METHODS}, got {method!r}")
    if method == "random":
        return None
    if use_rep not in adata.obsm:
        available = ", ".join(adata.obsm.keys()) or "none"
        raise KeyError(
            f"use_rep={use_rep!r} not found in adata.obsm (available: {available}). "
            "Compute one first, e.g. sc.pp.pca(adata)."
        )
    if not adata.obs_names.is_unique:
        raise ValueError("method='similarity' requires unique obs_names "
                         "(see adata.obs_names_make_unique()).")
    emb = np.asarray(adata.obsm[use_rep], dtype=np.float64)
    if n_dims is not None:
        if n_dims < 1:
            raise ValueError(f"n_dims must be >= 1, got {n_dims}")
        emb = emb[:, :n_dims]
    if not np.isfinite(emb).all():
        raise ValueError(f"adata.obsm[{use_rep!r}] contains NaN or infinite values.")
    return emb


def _add_compactness(
    adata: ad.AnnData,
    assigned: pd.DataFrame,
    mc_meta: pd.DataFrame,
    group_vars: list[str],
    emb: np.ndarray,
    seed: int,
) -> pd.DataFrame:
    """Add per-metacell compactness and log it against a size-matched random baseline."""
    pos = adata.obs_names.get_indexer(assigned.index)
    sub_emb = emb[pos]
    mc_ids = assigned["metacell_id"].to_numpy()
    mc_meta = mc_meta.copy()
    mc_meta["compactness"] = pd.Series(_compactness(sub_emb, mc_ids))

    # Baseline: shuffle metacell labels within each stratum (same sizes).
    rng = np.random.default_rng(seed)
    shuffled = mc_ids.copy()
    for _, rows in assigned.groupby(group_vars, sort=False).indices.items():
        shuffled[rows] = rng.permutation(mc_ids[rows])
    baseline = np.median(list(_compactness(sub_emb, shuffled).values()))
    observed = float(np.median(mc_meta["compactness"]))
    print(
        f"[scagg] median metacell compactness: {observed:.3g} "
        f"(random metacells of the same sizes: {baseline:.3g})",
        flush=True,
    )
    return mc_meta


# ── Internal builder ─────────────────────────────────────────────────────────

def _build_result(
    adata: ad.AnnData,
    obs_with_mc: pd.DataFrame,
    group_vars: list[str],
    meta_vars: Sequence[str] | None,
    layer: str | None,
    save_membership: str | Path | None,
    strategy: str,
    emb: np.ndarray | None = None,
    seed: int = 42,
) -> ad.AnnData:
    assigned = obs_with_mc[obs_with_mc["metacell_id"] != "unassigned"]
    n_unassigned = len(obs_with_mc) - len(assigned)
    if n_unassigned > 0:
        import warnings
        warnings.warn(
            f"{n_unassigned} cells excluded (below min_cells threshold).",
            stacklevel=3,
        )
    if len(assigned) == 0:
        import warnings
        warnings.warn("All cells excluded — no metacells formed. Returning empty AnnData.",
                      stacklevel=3)
        empty = ad.AnnData(obs=pd.DataFrame(), var=adata.var.copy())
        empty.uns["metacell_membership"] = _membership_table(obs_with_mc, "metacell_id")
        empty.uns["scagg_strategy"] = strategy
        empty.uns["scagg_group_vars"] = group_vars
        return empty

    membership = _membership_table(obs_with_mc, "metacell_id")

    X = adata[assigned.index].X if layer is None else adata[assigned.index].layers[layer]
    mc_ids = assigned["metacell_id"].values
    unique_ids, agg_X = _sum_by_group(X, mc_ids)

    if meta_vars is None:
        meta_vars_list = [c for c in adata.obs.columns if c not in group_vars]
    else:
        meta_vars_list = list(meta_vars)

    mc_meta = _aggregate_metadata(assigned, "metacell_id", meta_vars_list, group_vars)
    mc_meta = mc_meta.loc[unique_ids]   # ensure same order as agg_X rows
    if emb is not None:
        mc_meta = _add_compactness(adata, assigned, mc_meta, group_vars, emb, seed)

    new_adata = ad.AnnData(
        X=agg_X,
        obs=mc_meta,
        var=adata.var.copy(),
    )
    new_adata.uns["metacell_membership"] = membership
    new_adata.uns["scagg_strategy"] = strategy
    new_adata.uns["scagg_group_vars"] = group_vars

    if save_membership is not None:
        path = Path(save_membership)
        path.parent.mkdir(parents=True, exist_ok=True)
        membership.to_csv(path)
        print(f"[scagg] Membership saved → {path}", flush=True)

    print(
        f"[scagg] {strategy}: {len(adata)} cells → {len(new_adata)} metacells "
        f"({n_unassigned} excluded)",
        flush=True,
    )
    return new_adata
