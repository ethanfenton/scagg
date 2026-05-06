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
from typing import Sequence

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
) -> ad.AnnData:
    """Randomly aggregate cells into metacells of approximately *cell_size* cells.

    Within each group (defined by *group_vars*) cells are shuffled and then
    assigned round-robin to ``floor(n / cell_size)`` metacells.  Groups with
    fewer than *cell_min* cells are excluded.

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

    Returns
    -------
    AnnData with one obs per metacell.
    """
    if isinstance(group_vars, str):
        group_vars = [group_vars]
    group_vars = list(group_vars)

    if cell_min is None:
        cell_min = _default_cell_min(cell_size)

    rng = np.random.default_rng(seed)
    obs = adata.obs.copy()
    obs = _assign_by_size(obs, group_vars, cell_size, cell_min, rng, "metacell_id")

    return _build_result(
        adata=adata,
        obs_with_mc=obs,
        group_vars=group_vars,
        meta_vars=meta_vars,
        layer=layer,
        save_membership=save_membership,
        strategy=f"by_size_{cell_size}",
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
) -> ad.AnnData:
    """Randomly aggregate cells into exactly *n_metacells* metacells per group.

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

    Returns
    -------
    AnnData with one obs per metacell.
    """
    if isinstance(group_vars, str):
        group_vars = [group_vars]
    group_vars = list(group_vars)

    if min_cells is None:
        min_cells = n_metacells

    rng = np.random.default_rng(seed)
    obs = adata.obs.copy()
    obs = _assign_by_num(obs, group_vars, n_metacells, min_cells, rng, "metacell_id")

    return _build_result(
        adata=adata,
        obs_with_mc=obs,
        group_vars=group_vars,
        meta_vars=meta_vars,
        layer=layer,
        save_membership=save_membership,
        strategy=f"by_num_{n_metacells}",
    )


# ── Internal builder ─────────────────────────────────────────────────────────

def _build_result(
    adata: ad.AnnData,
    obs_with_mc: pd.DataFrame,
    group_vars: list[str],
    meta_vars: Sequence[str] | None,
    layer: str | None,
    save_membership: str | Path | None,
    strategy: str,
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
