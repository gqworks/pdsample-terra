"""Download helpers for brain-expression resources.

Nothing in the default pipeline requires these files. A missing table leaves
LoF weights at 1.

    python scripts/fetch_expression_resources.py print-urls
    python scripts/fetch_expression_resources.py gtex --input downloads/gtex.gct.gz --output resources/expression/brain_expression.tsv
    python scripts/fetch_expression_resources.py pext --input pext.tsv --output resources/expression/brain_expression.tsv
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from pdvar.expression_resources import build_expression_table, gtex_brain_median, print_urls, read_pext


def main() -> None:
    parser = argparse.ArgumentParser(description="Build a brain expression TSV for pdvar isoform weights.")
    sub = parser.add_subparsers(dest="command", required=True)
    sub.add_parser("print-urls", help="Print GTEx and gnomAD pext locations. Does not download.")
    gtex = sub.add_parser("gtex", help="Convert a GTEx v8 gene-median GCT into the expression TSV.")
    gtex.add_argument("--input", required=True)
    gtex.add_argument("--output", required=True)
    pext = sub.add_parser("pext", help="Convert a user-exported pext TSV (chrom, start, end, gene, brain_pext).")
    pext.add_argument("--input", required=True)
    pext.add_argument("--output", required=True)
    both = sub.add_parser("build", help="Concatenate a GTEx table and a pext table.")
    both.add_argument("--gtex")
    both.add_argument("--pext")
    both.add_argument("--output", required=True)
    args = parser.parse_args()
    if args.command == "print-urls":
        print(print_urls())
        return
    if args.command == "gtex":
        table = gtex_brain_median(args.input)
    elif args.command == "pext":
        table = read_pext(args.input)
    else:
        gtex_table = gtex_brain_median(args.gtex) if args.gtex else None
        pext_table = read_pext(args.pext) if args.pext else None
        table = build_expression_table(gtex_table, pext_table)
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(output, sep="\t", index=False)
    print(f"wrote {len(table):,} rows -> {output}")


if __name__ == "__main__":
    main()
