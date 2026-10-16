"""Core aggregation helpers shared by all three strategies."""

from __future__ import annotations

import numpy as np
import pandas as pd
from math import floor
from typing import Sequence

try:
    from scipy.sparse import issparse
except ImportError:
    def issparse(x):  # type: ignore[misc]
        return False

from ._partition import balanced_partition


# ── Cell-to-metacell assignment ───────────────────────────────────────────────

def _assign_by_size(
    obs: pd.DataFrame,
    group_vars: Sequence[str],
    cell_size: int,
    cell_min: int,
    rng: np.random.Generator,
    metacell_col: str,
    emb: np.ndarray | None = None,
    refine: bool = True,
    n_init_pcs: int = 1,
) -> pd.DataFrame:
    obs = obs.copy()
    obs[metacell_col] = "unassigned"

    for keys, grp in obs.groupby(list(group_vars), sort=False):
        n = len(grp)
        if n < cell_min:
            continue
        n_mc = 1 if n < cell_size else n // cell_size
        labels = _group_labels(obs, grp, n_mc, rng, emb, refine, n_init_pcs)
        key_str = _key_str(keys)
        obs.loc[grp.index, metacell_col] = [f"mc{l}_{key_str}" for l in labels]

    return obs


def _assign_by_num(
    obs: pd.DataFrame,
    group_vars: Sequence[str],
    n_metacells: int,
    min_cells: int,
    rng: np.random.Generator,
    metacell_col: str,
    emb: np.ndarray | None = None,
    refine: bool = True,
    n_init_pcs: int = 1,
) -> pd.DataFrame:
    obs = obs.copy()
    obs[metacell_col] = "unassigned"

    for keys, grp in obs.groupby(list(group_vars), sort=False):
        n = len(grp)
        if n < min_cells:
            continue
        actual_n = min(n_metacells, n)
        labels = _group_labels(obs, grp, actual_n, rng, emb, refine, n_init_pcs)
        key_str = _key_str(keys)
        obs.loc[grp.index, metacell_col] = [f"mc{l}_{key_str}" for l in labels]

    return obs


def _group_labels(obs, grp, n_mc, rng, emb, refine, n_init_pcs) -> np.ndarray:
    """Metacell labels 1..n_mc for the cells of one group.

    ``emb is None`` → random balanced assignment; otherwise similarity-based
    balanced partition of the group's rows of *emb*.
    """
    n = len(grp)
    if emb is None:
        labels = np.tile(np.arange(1, n_mc + 1), -(-n // n_mc))[:n]
        return rng.permutation(labels)
    pos = obs.index.get_indexer(grp.index)
    return balanced_partition(emb[pos], n_mc, refine=refine, n_init_pcs=n_init_pcs) + 1


def _assign_pseudobulk(
    obs: pd.DataFrame,
    group_vars: Sequence[str],
    metacell_col: str,
) -> pd.DataFrame:
    obs = obs.copy()
    obs[metacell_col] = obs[list(group_vars)].astype(str).agg("_".join, axis=1)
    return obs


# ── Expression aggregation ────────────────────────────────────────────────────

def _sum_by_group(X, group_labels: np.ndarray):
    """Return (groups, matrix) where matrix[i] = sum over the rows of group i.

    Rows labelled ``None`` are skipped.  The sum is one sparse indicator
    product, so time and memory scale with the input, not groups × cells.
    """
    from scipy.sparse import csr_matrix

    group_labels = np.asarray(group_labels, dtype=object)
    keep = np.flatnonzero(pd.notna(group_labels))
    codes, order = pd.factorize(group_labels[keep])   # first-appearance order
    dtype = X.dtype if np.issubdtype(X.dtype, np.floating) else np.int64
    ind = csr_matrix(
        (np.ones(len(keep), dtype=dtype), (codes, keep)),
        shape=(len(order), X.shape[0]),
    )
    out = ind @ X
    out = out.toarray() if issparse(out) else np.asarray(out)
    return order, out


# ── Metadata aggregation ──────────────────────────────────────────────────────

def _aggregate_metadata(
    obs: pd.DataFrame,
    metacell_col: str,
    meta_vars: list[str],
    group_vars: list[str],
) -> pd.DataFrame:
    """One row per metacell; mode for categorical, mean for numeric."""
    keep = [c for c in meta_vars + group_vars if c in obs.columns and c != metacell_col]
    records = []
    for mc_id, grp in obs.groupby(metacell_col, sort=False):
        row: dict = {"metacell_id": mc_id, "n_cells": len(grp)}
        for col in keep:
            s = grp[col]
            if pd.api.types.is_numeric_dtype(s):
                row[col] = s.mean()
            else:
                mode = s.mode()
                row[col] = mode.iloc[0] if len(mode) > 0 else None
        records.append(row)
    return pd.DataFrame(records).set_index("metacell_id")


# ── Membership table ─────────────────────────────────────────────────────────

def _membership_table(obs: pd.DataFrame, metacell_col: str) -> pd.DataFrame:
    """Map every cell barcode to its metacell."""
    cols = [metacell_col] + [c for c in obs.columns if c != metacell_col]
    t = obs[cols].copy()
    t.index.name = "cell_barcode"
    return t


# ── Safety check ─────────────────────────────────────────────────────────────

def _check_metacell_purity(obs: pd.DataFrame, group_vars: Sequence[str],
                           metacell_col: str = "metacell_id") -> None:
    """Raise unless every metacell holds cells from one *group_vars* combination.

    Different groups can map to the same ID when values contain ``"_"``
    (``("A_1", "x")`` and ``("A", "1_x")`` both give ``"A_1_x"``), which would
    otherwise sum cells from different samples together.
    """
    assigned = obs[obs[metacell_col] != "unassigned"]
    if assigned.empty:
        return
    n_groups = (assigned[list(group_vars)].astype(str)
                .groupby(assigned[metacell_col].to_numpy())
                .nunique().max(axis=1))
    bad = n_groups.index[n_groups > 1]
    if len(bad) > 0:
        raise ValueError(
            f"{len(bad)} metacell ID(s) would combine cells from different "
            f"{' × '.join(group_vars)} groups (e.g. {bad[0]!r}), because group "
            f"values joined with '_' collide. Rename the values in "
            f"{list(group_vars)} so that no joined combination is ambiguous."
        )


# ── Utility ──────────────────────────────────────────────────────────────────

def _key_str(keys) -> str:
    if isinstance(keys, (str, int, float)):
        return str(keys)
    return "_".join(str(k) for k in keys)


def _default_cell_min(cell_size: int) -> int:
    if cell_size == 10:
        return 7
    if cell_size == 25:
        return 15
    return floor(cell_size * 0.7)
