"""Undo log-normalisation, x = log1p(count / total * scale_factor), so values
can be combined in linear space: exact raw counts, or per-cell normalised
values for averaging.
"""

from __future__ import annotations

import warnings

import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix, diags, issparse

# obs columns that commonly hold per-cell total counts, tried in order.
_TOTAL_COLS = ("total_counts", "n_counts", "nCount_RNA")


def check_value_args(from_log: bool, mean_normalized: bool) -> None:
    if mean_normalized and not from_log:
        raise ValueError("mean_normalized=True requires from_log=True "
                         "(it averages log-normalised values in linear space).")


def expm1(X):
    """``expm1`` of a (sparse) matrix, as float64, keeping sparsity."""
    if issparse(X):
        Y = csr_matrix(X, dtype=np.float64, copy=True)
        np.expm1(Y.data, out=Y.data)
        return Y
    return np.expm1(np.asarray(X, dtype=np.float64))


def _row_min_positive(Y) -> np.ndarray:
    """Smallest positive value in each row (NaN for rows without one)."""
    if issparse(Y):
        data = np.where(Y.data > 0, Y.data, np.inf)
        out = np.full(Y.shape[0], np.inf)
        nonempty = np.diff(Y.indptr) > 0
        if nonempty.any():
            out[nonempty] = np.minimum.reduceat(data, Y.indptr[:-1][nonempty])
    else:
        out = np.where(Y > 0, Y, np.inf).min(axis=1)
    out[~np.isfinite(out)] = np.nan
    return out


def _row_max(Y, values) -> np.ndarray:
    """Per-row max of *values* laid out like Y's non-zeros (0 for empty rows)."""
    out = np.zeros(Y.shape[0])
    nonempty = np.diff(Y.indptr) > 0
    if nonempty.any():
        out[nonempty] = np.maximum.reduceat(values, Y.indptr[:-1][nonempty])
    return out


def count_unit(Y, max_multiple: int = 10, tol: float = 0.05) -> np.ndarray:
    """Per-row value of a single count in linear-space data *Y* (cells × genes).

    Every non-zero value is a whole multiple of it, so the unit is
    min / c for the smallest c in 1..max_multiple that makes all of the row's
    values whole multiples.  c = 1 is the plain smallest-non-zero trick; larger
    c handles cells with no count of exactly 1.  NaN for rows without counts.
    """
    Y = csr_matrix(Y)
    min_pos = _row_min_positive(Y)
    unit = min_pos.copy()
    todo = np.isfinite(min_pos)
    row = np.repeat(np.arange(Y.shape[0]), np.diff(Y.indptr))
    for c in range(1, max_multiple + 1):
        if not todo.any():
            break
        u = min_pos / c
        ratio = Y.data / u[row]
        dev = np.where(todo[row], np.abs(ratio - np.rint(ratio)), 0.0)
        ok = todo & (_row_max(Y, dev) < tol)
        unit[ok] = u[ok]
        todo &= ~ok
    unit[todo] = min_pos[todo]          # unresolved: plain trick (checked later)
    return unit


def _get_totals(obs: pd.DataFrame, total_col: str | None) -> np.ndarray | None:
    if total_col is not None:
        if total_col not in obs.columns:
            raise KeyError(f"total_col={total_col!r} not found in obs.")
        return pd.to_numeric(obs[total_col], errors="coerce").to_numpy(dtype=np.float64)
    for col in _TOTAL_COLS:
        if col in obs.columns:
            print(f"[scagg] from_log: using per-cell totals from obs[{col!r}]", flush=True)
            return pd.to_numeric(obs[col], errors="coerce").to_numpy(dtype=np.float64)
    return None


def recover_counts(X, obs: pd.DataFrame, total_col: str | None,
                   scale_factor: float | None):
    """Raw counts from log-normalised *X* (cells × genes).

    count = expm1(x) * total / scale_factor, using each cell's total count.
    Without totals, falls back to assuming the smallest non-zero value in each
    cell is a count of 1.  *scale_factor* ``None`` infers it from the data
    (Seurat uses 1e4; scanpy's ``normalize_total`` defaults to the median
    total).
    """
    Y = expm1(X)
    unit = count_unit(Y)                  # value of one count in each cell
    totals = _get_totals(obs, total_col)

    if totals is None:
        print("[scagg] from_log: no per-cell totals found (looked for "
              f"{', '.join(_TOTAL_COLS)}; set total_col). Falling back to the "
              "smallest-non-zero trick: each cell's values are scaled so its "
              "smallest non-zero value is the smallest whole count consistent "
              "with the rest (usually 1).", flush=True)
        size = unit
    else:
        valid = np.isfinite(totals) & (totals > 0)
        if scale_factor is None:
            est = (unit * totals)[valid & np.isfinite(unit)]
            if est.size == 0:
                raise ValueError("Cannot infer scale_factor: no cell has both a "
                                 "valid total and a non-zero value.")
            scale_factor = float(np.median(est))
            print(f"[scagg] from_log: inferred scale_factor = {scale_factor:.6g}",
                  flush=True)
        size = np.full(len(totals), np.nan)
        size[valid] = scale_factor / totals[valid]
        if not valid.all():
            print(f"[scagg] from_log: {int((~valid).sum())} cell(s) have a missing "
                  "or zero total; falling back to the smallest-non-zero trick "
                  "for them.", flush=True)
            size[~valid] = unit[~valid]

    size = np.where(np.isfinite(size) & (size > 0), size, 1.0)   # all-zero cells
    C = diags(1.0 / size) @ Y if issparse(Y) else Y / size[:, None]
    if issparse(C):
        C = csr_matrix(C)
    vals = C.data if issparse(C) else C
    dev = float(np.abs(vals - np.rint(vals)).max()) if vals.size else 0.0
    if dev > 0.1:
        warnings.warn(
            f"Recovered counts are up to {dev:.3g} away from whole numbers, so the "
            "input may not be log1p(count / total * scale_factor) (e.g. it was "
            "scaled, batch-corrected, or used a different total). Check "
            "total_col and scale_factor.",
            stacklevel=4,
        )
    np.rint(vals, out=vals)
    return C


def warn_if_not_counts(X) -> None:
    """Warn when the matrix about to be summed does not look like counts."""
    vals = X.data if issparse(X) else np.asarray(X).ravel()
    vals = vals[:1_000_000]
    if vals.size and (np.any(vals < 0) or np.any(vals != np.rint(vals))):
        warnings.warn(
            "The matrix being summed has non-integer or negative values. Sums "
            "are only meaningful for counts; if this is log-normalised data, "
            "pass from_log=True (or layer='counts').",
            stacklevel=4,
        )
