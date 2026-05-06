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


# ── Cell-to-metacell assignment ───────────────────────────────────────────────

def _assign_by_size(
    obs: pd.DataFrame,
    group_vars: Sequence[str],
    cell_size: int,
    cell_min: int,
    rng: np.random.Generator,
    metacell_col: str,
) -> pd.DataFrame:
    obs = obs.copy()
    obs[metacell_col] = "unassigned"

    for keys, grp in obs.groupby(list(group_vars), sort=False):
        n = len(grp)
        if n < cell_min:
            continue
        n_mc = 1 if n < cell_size else n // cell_size
        labels = np.tile(np.arange(1, n_mc + 1), -(-n // n_mc))[:n]
        labels = rng.permutation(labels)
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
) -> pd.DataFrame:
    obs = obs.copy()
    obs[metacell_col] = "unassigned"

    for keys, grp in obs.groupby(list(group_vars), sort=False):
        n = len(grp)
        if n < min_cells:
            continue
        actual_n = min(n_metacells, n)
        labels = np.tile(np.arange(1, actual_n + 1), -(-n // actual_n))[:n]
        labels = rng.permutation(labels)
        key_str = _key_str(keys)
        obs.loc[grp.index, metacell_col] = [f"mc{l}_{key_str}" for l in labels]

    return obs


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
    """Return (groups, matrix) where matrix[i] = sum over group i."""
    unique_groups = pd.unique(group_labels)  # preserves order
    rows = []
    for g in unique_groups:
        mask = group_labels == g
        block = X[mask]
        if issparse(block):
            rows.append(np.asarray(block.sum(axis=0)).flatten())
        else:
            rows.append(np.asarray(block).sum(axis=0))
    return unique_groups, np.vstack(rows)


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
