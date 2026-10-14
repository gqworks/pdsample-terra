#!/usr/bin/env python3
"""Look up AlphaGenome Atlas AVI scores for SNVs in a VCF.

Reads ``ALPHAGENOME_API_KEY`` from the environment. The key is not printed.
An empty key writes a header-only TSV and exits 0 so the rest of PdSample
still ranks variants. Atlas (alphagenome 0.10.0) serves GRCh38 SNVs only.

Phred is the genome-wide tail rank published with the Atlas:
``-10 * log10(1 - CDF)``. Phred 20 is the top 1% of SNVs and Phred 30 is the
top 0.1%. It is not a pathogenicity call.
"""

from __future__ import annotations

import argparse
import csv
import math
import os
import sys
from pathlib import Path

HEADER = ("chrom", "pos", "ref", "alt", "avi_raw", "avi_quantile", "avi_phred", "avi_url")
_ACGT = set("ACGT")
_PORTAL = "https://deepmind.google.com/science/alphagenome/atlas"


def avi_phred(cdf: float) -> float:
    """Phred-scaled genome-wide tail from an Atlas CDF quantile."""
    value = float(cdf)
    if not math.isfinite(value) or value < 0.0 or value > 1.0:
        raise ValueError(f"AVI CDF quantile must be in [0, 1], got {cdf}")
    tail = max(1.0 - value, 1e-10)
    # Six decimals keeps Phred 20 and 30 exact when the CDF is 0.99 or 0.999.
    return round(-10.0 * math.log10(tail), 6)


def atlas_chrom(chrom: str) -> str:
    """chr-prefixed GRCh38 name. ``MT`` becomes ``chrM``."""
    text = str(chrom).strip()
    rest = text[3:] if text.lower().startswith("chr") else text
    if rest.upper() in {"M", "MT"}:
        return "chrM"
    return f"chr{rest}"


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


def _fatal_lookup(exc: BaseException) -> bool:
    if isinstance(exc, (PermissionError, TimeoutError)):
        return True
    text = str(exc).lower()
    return "api key" in text or "unauthenticated" in text or "permission" in text


def _row_from_scores(chrom: str, pos: int, ref: str, alt: str, scores) -> dict[str, object] | None:
    adata = scores.get("AVI_SCORE") if scores else None
    if adata is None or getattr(adata, "n_obs", 0) < 1:
        return None
    returned = adata.obs["variant"].iloc[0]
    expected = f"{atlas_chrom(chrom)}:{pos}:{ref}>{alt}"
    if str(returned) != expected:
        print(f"AVI returned {returned}, not {expected}", file=sys.stderr)
        return None
    raw = float(adata.X[0, 0])
    layers = getattr(adata, "layers", None) or {}
    if "quantiles" not in layers:
        print(f"AVI response for {expected} has no quantile layer", file=sys.stderr)
        return None
    cdf = float(layers["quantiles"][0, 0])
    return {
        "chrom": atlas_chrom(chrom),
        "pos": pos,
        "ref": ref,
        "alt": alt,
        "avi_raw": raw,
        "avi_quantile": cdf,
        "avi_phred": avi_phred(cdf),
        "avi_url": atlas_url(chrom, pos, ref, alt),
    }


def lookup_snvs(vcf: Path, api_key: str):
    from alphagenome.atlas import atlas
    from alphagenome.data import genome

    client = atlas.create(api_key, timeout=30)
    rows = []
    seen = set()
    for chrom, pos, ref, alt in iter_snvs(vcf):
        key = (atlas_chrom(chrom), pos, ref, alt)
        if key in seen:
            continue
        seen.add(key)
        variant = genome.Variant.from_str(f"{key[0]}:{pos}:{ref}>{alt}")
        try:
            scores = client.query_variant(variant, requested_scorers=["AVI_SCORE"])
        except (ValueError, PermissionError, TimeoutError) as exc:
            if _fatal_lookup(exc):
                raise
            print(f"AVI miss for {variant}: {exc}", file=sys.stderr)
            continue
        row = _row_from_scores(chrom, pos, ref, alt, scores)
        if row is not None:
            rows.append(row)
    return rows


def write_table(path: Path, rows: list[dict[str, object]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=HEADER, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Annotate SNVs with AlphaGenome Atlas AVI scores.")
    parser.add_argument("--vcf", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args(argv)
    key = os.environ.get("ALPHAGENOME_API_KEY", "").strip()
    if not key:
        print("ALPHAGENOME_API_KEY is empty. AVI scores were not queried.", file=sys.stderr)
        write_table(args.out, [])
        return 0
    if not args.vcf.exists():
        print(f"VCF not found: {args.vcf}", file=sys.stderr)
        return 1
    rows = lookup_snvs(args.vcf, key)
    write_table(args.out, rows)
    print(f"Wrote {len(rows)} AVI scores to {args.out}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
