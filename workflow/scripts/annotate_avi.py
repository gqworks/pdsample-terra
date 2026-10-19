#!/usr/bin/env python3
"""Look up AlphaGenome AVI scores for SNVs in a VCF from the precomputed table.

The table is the AlphaGenome "AVI SNV scores" download: a bgzipped, Tabix
indexed TSV of every GRCh38 SNV with ``CHROM POS REF ALT raw_score PHRED``.
Lookups are local. No API key and no network call are needed, so the task runs
inside restricted Terra / AnVIL workspaces.

No ``--table`` writes a header-only TSV and exits 0 so the rest of PdSample
still ranks variants.

Phred is the genome-wide tail rank published with the Atlas:
``-10 * log10(1 - CDF)``. Phred 20 is the top 1% of SNVs and Phred 30 is the
top 0.1%. It is not a pathogenicity call. The table does not carry the CDF, so
``avi_quantile`` is derived from the published Phred as
``1 - 10 ** (-Phred / 10)``. It is only as precise as the rounded Phred.
"""

from __future__ import annotations

import argparse
import csv
import math
import sys
from collections import defaultdict
from pathlib import Path

HEADER = ("chrom", "pos", "ref", "alt", "avi_raw", "avi_quantile", "avi_phred", "avi_url")
_ACGT = set("ACGT")
_PORTAL = "https://deepmind.google.com/science/alphagenome/atlas"

# Column names accepted in the table header, compared case-insensitively
# after stripping a leading '#'.
_COLUMN_ALIASES = {
    "chrom": ("chrom", "chr", "chromosome", "contig"),
    "pos": ("pos", "position"),
    "ref": ("ref",),
    "alt": ("alt",),
    "raw": ("raw_score", "avi_raw", "avi_score", "avi", "score", "raw"),
    "phred": ("phred", "avi_phred", "phred_score"),
}
# Order used when the table has no header line.
_DEFAULT_COLUMNS = ("chrom", "pos", "ref", "alt", "raw", "phred")
_MISSING = {"", ".", "na", "nan", "none", "null"}

# Positions closer than this are fetched in one Tabix query. The table has up
# to three rows per position, so a small gap keeps wasted parsing low while
# cutting the number of seeks for dense WGS call sets.
DEFAULT_MERGE_GAP = 250
# Fail when this share of looked-up positions has a table row with a different
# REF. That points at a build or table mix-up, not at missing scores.
MAX_REF_MISMATCH_FRACTION = 0.10
MIN_POSITIONS_FOR_MISMATCH_CHECK = 50


def phred_to_quantile(phred: float) -> float:
    """Genome-wide CDF quantile implied by a published Phred."""
    value = float(phred)
    if not math.isfinite(value) or value < 0.0:
        raise ValueError(f"AVI Phred must be finite and >= 0, got {phred}")
    return round(1.0 - 10.0 ** (-value / 10.0), 10)


def atlas_chrom(chrom: str) -> str:
    """chr-prefixed GRCh38 name. ``MT`` becomes ``chrM``."""
    text = str(chrom).strip()
    rest = text[3:] if text.lower().startswith("chr") else text
    if rest.upper() in {"M", "MT"}:
        return "chrM"
    return f"chr{rest}"


def _bare_chrom(chrom: str) -> str:
    name = atlas_chrom(chrom)[3:]
    return "MT" if name == "M" else name


def atlas_url(chrom: str, pos: int, ref: str, alt: str) -> str:
    query = f"{atlas_chrom(chrom)}:{pos}:{ref}>{alt}"
    return f"{_PORTAL}?m=variant&q={query}"


def _is_snv(ref: str, alt: str) -> bool:
    return len(ref) == 1 and len(alt) == 1 and ref in _ACGT and alt in _ACGT and ref != alt


def iter_snvs(path: Path):
    """PASS SNVs as ``(chrom, pos, ref, alt)``. Multi-allelic rows are split."""
    if str(path).endswith(".gz"):
        import gzip

        handle_cm = gzip.open(path, "rt", encoding="utf-8", errors="replace")
    else:
        handle_cm = path.open("r", encoding="utf-8", errors="replace")
    with handle_cm as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 8:
                continue
            chrom, pos_text, _rsid, ref, alts, _qual, filt = parts[:7]
            if filt not in {"PASS", "."}:
                continue
            try:
                pos = int(pos_text)
            except ValueError:
                continue
            ref = ref.upper()
            for alt in alts.split(","):
                alt = alt.upper()
                if _is_snv(ref, alt):
                    yield chrom, pos, ref, alt


def column_indices(header_lines) -> dict[str, int]:
    """Map the needed fields to column numbers from the table header."""
    names = None
    for line in header_lines:
        text = line.decode() if isinstance(line, bytes) else str(line)
        if text.startswith("#"):
            names = text.lstrip("#").rstrip("\n").split("\t")
    if names is None:
        print(
            "AVI table has no header line. Assuming " + " ".join(_DEFAULT_COLUMNS) + ".",
            file=sys.stderr,
        )
        return {field: index for index, field in enumerate(_DEFAULT_COLUMNS)}
    lowered = [name.strip().lower() for name in names]
    found: dict[str, int] = {}
    for field, aliases in _COLUMN_ALIASES.items():
        for alias in aliases:
            if alias in lowered:
                found[field] = lowered.index(alias)
                break
    missing = [field for field in _COLUMN_ALIASES if field not in found]
    if missing:
        raise ValueError(f"AVI table header {names} has no column for {missing}")
    return found


def _number(text: str) -> float | None:
    if text.strip().lower() in _MISSING:
        return None
    value = float(text)
    return value if math.isfinite(value) else None


def _windows(positions: list[int], gap: int):
    """Merge sorted 1-based positions into inclusive ``(start, end)`` windows."""
    start = end = positions[0]
    for pos in positions[1:]:
        if pos - end <= gap:
            end = pos
        else:
            yield start, end
            start = end = pos
    yield start, end


def _table_contig(chrom: str, contigs: set[str]) -> str | None:
    for name in (atlas_chrom(chrom), _bare_chrom(chrom), str(chrom).strip()):
        if name in contigs:
            return name
    return None


def lookup_snvs(vcf: Path, table: Path, index: Path | None = None, merge_gap: int = DEFAULT_MERGE_GAP):
    """AVI rows for the PASS SNVs in ``vcf``, in VCF order. Misses are dropped."""
    import pysam

    snvs = []
    seen = set()
    for chrom, pos, ref, alt in iter_snvs(vcf):
        key = (atlas_chrom(chrom), pos, ref, alt)
        if key in seen:
            continue
        seen.add(key)
        snvs.append((chrom, pos, ref, alt))

    tbx = pysam.TabixFile(str(table), index=str(index) if index else None)
    try:
        cols = column_indices(tbx.header)
        width = max(cols.values()) + 1
        contigs = set(tbx.contigs)

        # table contig -> {pos: {alt: key}} for the wanted SNVs.
        wanted: dict[str, dict[int, dict[str, tuple]]] = defaultdict(lambda: defaultdict(dict))
        ref_at: dict[tuple[str, int], str] = {}
        absent_contigs = set()
        for chrom, pos, ref, alt in snvs:
            contig = _table_contig(chrom, contigs)
            if contig is None:
                absent_contigs.add(chrom)
                continue
            wanted[contig][pos][alt] = (atlas_chrom(chrom), pos, ref, alt)
            ref_at[(contig, pos)] = ref

        found: dict[tuple, tuple[float | None, float | None]] = {}
        checked = set()
        mismatched = set()
        for contig, by_pos in wanted.items():
            for start, end in _windows(sorted(by_pos), merge_gap):
                # Tabix regions are 0-based half-open; VCF POS is 1-based.
                for line in tbx.fetch(contig, start - 1, end):
                    parts = line.split("\t")
                    if len(parts) < width:
                        continue
                    try:
                        pos = int(parts[cols["pos"]])
                    except ValueError:
                        continue
                    alts = by_pos.get(pos)
                    if alts is None:
                        continue
                    ref = parts[cols["ref"]].upper()
                    checked.add((contig, pos))
                    if ref != ref_at[(contig, pos)]:
                        mismatched.add((contig, pos))
                        continue
                    key = alts.get(parts[cols["alt"]].upper())
                    if key is None:
                        continue
                    raw = _number(parts[cols["raw"]])
                    phred = _number(parts[cols["phred"]])
                    if raw is None and phred is None:
                        continue
                    found[key] = (raw, phred)
    finally:
        tbx.close()

    if absent_contigs:
        print(f"AVI table has no contig for {sorted(absent_contigs)}", file=sys.stderr)
    if mismatched:
        share = len(mismatched) / max(len(checked), 1)
        print(
            f"AVI table REF differs from the VCF at {len(mismatched)} of {len(checked)} positions",
            file=sys.stderr,
        )
        if len(checked) >= MIN_POSITIONS_FOR_MISMATCH_CHECK and share > MAX_REF_MISMATCH_FRACTION:
            raise ValueError(
                "AVI table REF disagrees with the VCF at "
                f"{share:.1%} of positions. Check that both are GRCh38."
            )

    rows = []
    for chrom, pos, ref, alt in snvs:
        key = (atlas_chrom(chrom), pos, ref, alt)
        if key not in found:
            continue
        raw, phred = found[key]
        rows.append(
            {
                "chrom": key[0],
                "pos": pos,
                "ref": ref,
                "alt": alt,
                "avi_raw": "" if raw is None else raw,
                "avi_quantile": "" if phred is None else phred_to_quantile(phred),
                "avi_phred": "" if phred is None else phred,
                "avi_url": atlas_url(chrom, pos, ref, alt),
            }
        )
    print(f"AVI table scored {len(rows)} of {len(snvs)} PASS SNVs", file=sys.stderr)
    return rows


def write_table(path: Path, rows: list[dict[str, object]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=HEADER, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Annotate SNVs with AlphaGenome AVI scores from the precomputed Tabix table."
    )
    parser.add_argument("--vcf", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--table", type=Path, help="bgzipped AVI SNV table (.tsv.gz)")
    parser.add_argument("--index", type=Path, help="Tabix index. Default: <table>.tbi or .csi")
    parser.add_argument("--merge-gap", type=int, default=DEFAULT_MERGE_GAP)
    args = parser.parse_args(argv)
    if args.table is None or not str(args.table).strip():
        print("No AVI table given. AVI scores were not looked up.", file=sys.stderr)
        write_table(args.out, [])
        return 0
    for label, path in (("VCF", args.vcf), ("AVI table", args.table), ("AVI index", args.index)):
        if path is not None and not path.exists():
            print(f"{label} not found: {path}", file=sys.stderr)
            return 1
    rows = lookup_snvs(args.vcf, args.table, args.index, args.merge_gap)
    write_table(args.out, rows)
    print(f"Wrote {len(rows)} AVI scores to {args.out}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
