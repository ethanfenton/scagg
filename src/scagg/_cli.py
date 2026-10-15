"""Command-line interface for scagg."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from . import __version__


def _build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="scagg",
        description=(
            "scagg: pseudobulk and metacell aggregation for single-cell data.\n\n"
            "Sub-commands:\n"
            "  pseudobulk          Aggregate all cells per group (true pseudobulk)\n"
            "  metacells-by-size   Make metacells of ~N cells each\n"
            "  metacells-by-num    Make exactly N metacells per group\n"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--version", action="version", version=f"scagg {__version__}")
    sub = p.add_subparsers(dest="command")

    _add_common = lambda sp: (
        sp.add_argument("--input", "-i", required=True, metavar="H5AD",
                        help="Input AnnData (.h5ad) file"),
        sp.add_argument("--output", "-o", required=True, metavar="H5AD",
                        help="Output AnnData (.h5ad) file"),
        sp.add_argument("--group-vars", "-g", nargs="+", required=True, metavar="COL",
                        help="obs column(s) to group by (e.g. cell_type sample)"),
        sp.add_argument("--meta-vars", nargs="*", metavar="COL",
                        help="obs columns to carry into metacell metadata "
                             "(default: all obs columns)"),
        sp.add_argument("--layer", metavar="LAYER",
                        help="AnnData layer to aggregate (default: X)"),
        sp.add_argument("--save-membership", metavar="CSV",
                        help="Path to save cell → metacell membership CSV"),
    )

    def _add_method(sp):
        sp.add_argument("--method", choices=["random", "similarity"], default="random",
                        help="How cells are grouped within each stratum: randomly, or "
                             "by similarity in an embedding [default: random]")
        sp.add_argument("--use-rep", default="X_pca", metavar="KEY",
                        help="obsm embedding for --method similarity [default: X_pca]")
        sp.add_argument("--n-dims", type=int, default=None, metavar="N",
                        help="Use the first N embedding dimensions (default: all)")
        sp.add_argument("--no-refine", action="store_true",
                        help="Skip size-preserving refinement for --method similarity")
        sp.add_argument("--n-init-pcs", type=int, default=1, metavar="N",
                        help="For --method similarity, try each split from PC1..PCN "
                             "and keep the tightest [default: 1]")

    # pseudobulk
    sp_pb = sub.add_parser("pseudobulk", help="True pseudobulk aggregation")
    _add_common(sp_pb)

    # metacells-by-size
    sp_sz = sub.add_parser("metacells-by-size",
                            help="Create metacells of approximately N cells")
    _add_common(sp_sz)
    sp_sz.add_argument("--cell-size", "-c", type=int, required=True, metavar="N",
                       help="Target cells per metacell")
    sp_sz.add_argument("--cell-min", type=int, default=None, metavar="N",
                       help="Minimum cells to form a metacell group "
                            "(default: 70%% of --cell-size)")
    sp_sz.add_argument("--seed", type=int, default=42, metavar="N",
                       help="Random seed [default: 42]")
    _add_method(sp_sz)

    # metacells-by-num
    sp_num = sub.add_parser("metacells-by-num",
                             help="Create exactly N metacells per group")
    _add_common(sp_num)
    sp_num.add_argument("--n-metacells", "-n", type=int, required=True, metavar="N",
                        help="Number of metacells per group")
    sp_num.add_argument("--min-cells", type=int, default=None, metavar="N",
                        help="Minimum cells needed (default: --n-metacells)")
    sp_num.add_argument("--seed", type=int, default=42, metavar="N",
                        help="Random seed [default: 42]")
    _add_method(sp_num)

    return p


def main(argv: list[str] | None = None) -> None:
    parser = _build_parser()
    args = parser.parse_args(argv)

    if args.command is None:
        parser.print_help()
        sys.exit(0)

    try:
        import anndata as ad
    except ImportError:
        print("ERROR: anndata is required. Install with: pip install anndata",
              file=sys.stderr)
        sys.exit(1)

    from ._metacells import make_pseudobulk, make_metacells_by_size, make_metacells_by_num

    print(f"[scagg] Reading {args.input} ...", flush=True)
    adata = ad.read_h5ad(args.input)

    kwargs = dict(
        group_vars=args.group_vars,
        meta_vars=args.meta_vars,
        layer=args.layer,
        save_membership=args.save_membership,
    )

    if args.command == "pseudobulk":
        result = make_pseudobulk(adata, **kwargs)
    else:
        mc_kwargs = dict(
            seed=args.seed, method=args.method, use_rep=args.use_rep,
            n_dims=args.n_dims, refine=not args.no_refine,
            n_init_pcs=args.n_init_pcs, **kwargs,
        )
        if args.command == "metacells-by-size":
            result = make_metacells_by_size(
                adata, cell_size=args.cell_size, cell_min=args.cell_min, **mc_kwargs,
            )
        else:
            result = make_metacells_by_num(
                adata, n_metacells=args.n_metacells, min_cells=args.min_cells,
                **mc_kwargs,
            )

    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    result.write_h5ad(str(out))
    print(f"[scagg] Saved → {out}", flush=True)


if __name__ == "__main__":
    main()
