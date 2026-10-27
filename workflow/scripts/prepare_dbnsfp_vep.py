#!/usr/bin/env python3
"""Build a VEP-ready bgzip+tabix dbNSFP file from academic per-chromosome gz tables.

Academic dbNSFP 5.4a ships as ``dbNSFP5.4a_variant.chr{1-22,X,Y,M}.gz`` (plain gzip,
not tabix). The Ensembl VEP dbNSFP plugin needs one bgzip+tabix file. This script
concatenates the chrom tables (optionally restricted to a panel BED) and writes
``out_prefix.gz`` + ``out_prefix.gz.tbi``.
"""

from __future__ import annotations

import argparse
import gzip
import re
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

_CHR_RE = re.compile(r"(?:^|[._])chr([0-9]+|[XYM]|MT)(?:\.|$)", re.IGNORECASE)


def _norm_chrom(value: str) -> str:
    text = value.strip()
    if text.lower().startswith("chr"):
        text = text[3:]
    if text.upper() == "MT":
        return "M"
    return text.upper() if text.upper() in {"X", "Y", "M"} else text


def _chrom_from_path(path: Path) -> str:
    match = _CHR_RE.search(path.name)
    if not match:
        raise SystemExit(f"Could not parse chromosome from dbNSFP filename: {path.name}")
    return _norm_chrom(match.group(1))


def _load_panel_intervals(bed: Path | None) -> dict[str, list[tuple[int, int]]] | None:
    """Return 1-based inclusive intervals keyed by dbNSFP chrom, or None for genome-wide."""
    if bed is None:
        return None
    intervals: dict[str, list[tuple[int, int]]] = defaultdict(list)
    with bed.open(encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            parts = line.split("\t")
            if len(parts) < 3:
                continue
            chrom = _norm_chrom(parts[0])
            start0 = int(parts[1])
            end0 = int(parts[2])
            # BED half-open 0-based -> inclusive 1-based (dbNSFP pos(1-based)).
            start = start0 + 1
            end = end0
            if end >= start:
                intervals[chrom].append((start, end))
    for chrom, ranges in intervals.items():
        ranges.sort()
    return dict(intervals)


def _in_intervals(pos: int, ranges: list[tuple[int, int]]) -> bool:
    for start, end in ranges:
        if pos < start:
            return False
        if pos <= end:
            return True
    return False


def _open_text(path: Path):
    if path.name.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return path.open(encoding="utf-8", errors="replace")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out-prefix", required=True, help="Output prefix (writes .gz and .gz.tbi)")
    parser.add_argument(
        "--panel-bed",
        default="",
        help="Optional BED (chr-prefixed OK). When set, keep only panel positions.",
    )
    parser.add_argument(
        "chrom_files",
        nargs="+",
        help="Per-chromosome dbNSFP5.4a_variant.chr*.gz files",
    )
    args = parser.parse_args()

    chrom_paths = [Path(p) for p in args.chrom_files]
    for path in chrom_paths:
        if not path.is_file():
            raise SystemExit(f"Missing dbNSFP chrom file: {path}")

    intervals = _load_panel_intervals(Path(args.panel_bed) if args.panel_bed else None)
    out_prefix = Path(args.out_prefix)
    out_tsv = out_prefix.with_suffix(".tmp.tsv")
    out_gz = Path(str(out_prefix) + ".gz") if not str(out_prefix).endswith(".gz") else out_prefix

    # Stable chrom order: 1..22, X, Y, M, then any extras.
    rank = {str(i): i for i in range(1, 23)}
    rank.update({"X": 23, "Y": 24, "M": 25})
    chrom_paths = sorted(chrom_paths, key=lambda p: (rank.get(_chrom_from_path(p), 99), p.name))

    header: list[str] | None = None
    kept = 0
    with out_tsv.open("w", encoding="utf-8", newline="\n") as out:
        for path in chrom_paths:
            chrom = _chrom_from_path(path)
            if intervals is not None and chrom not in intervals:
                print(f"skip {path.name}: no panel intervals on chr{chrom}", file=sys.stderr)
                continue
            ranges = intervals.get(chrom, []) if intervals is not None else None
            with _open_text(path) as handle:
                first = handle.readline()
                if not first:
                    continue
                cols = first.lstrip("#").rstrip("\n").split("\t")
                if header is None:
                    header = cols
                    out.write("#" + "\t".join(header) + "\n")
                elif cols != header:
                    raise SystemExit(f"Header mismatch in {path.name}")
                try:
                    pos_idx = header.index("pos(1-based)")
                except ValueError as exc:
                    raise SystemExit("dbNSFP header missing pos(1-based)") from exc
                for line in handle:
                    if not line.strip():
                        continue
                    parts = line.rstrip("\n").split("\t")
                    if len(parts) <= pos_idx:
                        continue
                    if ranges is not None:
                        try:
                            pos = int(parts[pos_idx])
                        except ValueError:
                            continue
                        if not _in_intervals(pos, ranges):
                            continue
                    out.write(line if line.endswith("\n") else line + "\n")
                    kept += 1

    if header is None:
        raise SystemExit("No dbNSFP rows written (empty inputs?)")

    # Sort for tabix (chrom then pos). Panel extracts are small; full builds need disk for sort.
    sorted_tsv = Path(str(out_prefix) + ".sorted.tsv")
    with out_tsv.open(encoding="utf-8") as raw, sorted_tsv.open("w", encoding="utf-8", newline="\n") as sorted_out:
        sorted_out.write(raw.readline())
        body = subprocess.run(
            ["sort", "-k1,1", "-k2,2n"],
            check=True,
            input=raw.read(),
            capture_output=True,
            text=True,
        )
        sorted_out.write(body.stdout)
    out_tsv.unlink(missing_ok=True)
    with out_gz.open("wb") as gz_handle:
        subprocess.run(["bgzip", "-c", str(sorted_tsv)], check=True, stdout=gz_handle)
    sorted_tsv.unlink(missing_ok=True)
    subprocess.run(["tabix", "-f", "-s", "1", "-b", "2", "-e", "2", str(out_gz)], check=True)
    print(f"wrote {out_gz} ({kept} data rows)", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
