"""
scagg — Single-cell Aggregation

Create pseudobulk samples or metacells from AnnData / Seurat objects.
"""

from importlib.metadata import version, PackageNotFoundError

try:
    __version__ = version("scagg")
except PackageNotFoundError:
    __version__ = "0.0.0"

from ._metacells import (
    make_pseudobulk,
    make_metacells_by_size,
    make_metacells_by_num,
)

__all__ = [
    "make_pseudobulk",
    "make_metacells_by_size",
    "make_metacells_by_num",
    "__version__",
]
