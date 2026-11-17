#!/usr/bin/env python3
"""Rebuild MonoPDAus panel gene/exon BEDs from a GENCODE/Ensembl GTF.

Gene spans are the union of GENCODE-primary (or transcript_type=protein_coding
with tag containing 'GENCODE primary' / Ensembl canonical) exons for each panel
ENSG. That avoids non-primary transcripts that pull a panel gene into a
neighbour (e.g. OPA3-203 into GPR4).

Exon rows are exons belonging to those transcripts only (never every exon that
merely overlaps the gene interval), with optional padding.

If --gtf is omitted, downloads GENCODE v44 annotation to --gtf-cache.

Outputs (BED 0-based half-open, chr-prefixed Broad GRCh38 contigs):
  workflow/assets/monopdaus-genes.grch38.bed
      chrom start end symbol gene_id   (gene body + calling_flank)
  resources/panel/monopdaus-genes.grch38.bed
      chrom start end gene_id symbol   (gene body + sv_flank)
  resources/panel/monopdaus-exons.grch38.bed
      chrom start end gene_id symbol   (exons + exon_pad)
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import re
import sys
import time
import urllib.request
from collections import defaultdict
from pathlib import Path

GENCODE_V44_URL = (
    "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/"
    "release_44/gencode.v44.annotation.gtf.gz"
)

# Panel BED symbols use HGNC/clinical names (GBA1); inheritance CSV may still say GBA.
DEFAULT_ALIASES = {"GBA": "GBA1"}

ATTR_RE = re.compile(r'(\w+)\s+"([^"]*)"')


def _parse_attrs(blob: str) -> dict[str, str]:
    """Parse GTF attributes; repeated keys (e.g. tag) are comma-joined."""
    out: dict[str, str] = {}
    for key, value in ATTR_RE.findall(blob):
        if key in out:
            out[key] = f"{out[key]},{value}"
        else:
            out[key] = value
    return out


def _style_chrom(seq_region: str) -> str:
    if seq_region in {"MT", "M", "chrM", "chrMT"}:
        return "chrM"
    if seq_region.startswith("chr"):
        return seq_region
    return f"chr{seq_region}"


def load_panel_genes(path: Path, aliases: dict[str, str] | None = None) -> dict[str, str]:
    """Return gene_id -> panel symbol."""
    aliases = {**(DEFAULT_ALIASES), **(aliases or {})}
    out: dict[str, str] = {}
    with path.open(encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            gid = (row.get("Gene_ID") or row.get("gene_id") or "").strip()
            sym = (row.get("Symbol") or row.get("gene") or "").strip()
            if gid and sym:
                out[gid] = aliases.get(sym, sym)
    if not out:
        raise SystemExit(f"no genes in {path}")
    return out


def ensure_gtf(path: Path, url: str) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists() and path.stat().st_size > 0:
        return path
    print(f"downloading {url} -> {path}", flush=True)
    urllib.request.urlretrieve(url, path)
    return path


def _open_text(path: Path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open(encoding="utf-8")


def _is_primary(attrs: dict[str, str]) -> bool:
    tags = {t.strip() for t in (attrs.get("tag") or "").split(",") if t.strip()}
    return bool(tags & {"GENCODE primary", "gencode_primary", "Ensembl_canonical", "MANE_Select"})


def extract_from_gtf(
    gtf: Path,
    panel: dict[str, str],
) -> tuple[dict[str, tuple[str, int, int]], dict[str, set[tuple[int, int]]]]:
    """Return gene_id -> (chrom, start1, end1) and gene_id -> exon set (1-based).

    Prefer exons on GENCODE-primary / Ensembl_canonical protein-coding transcripts.
    Fall back to all protein-coding transcripts, then any exon for the gene.
    """
    # gene_id -> transcript_id -> list of (chrom, start1, end1, primary, protein_coding)
    tx_exons: dict[str, dict[str, list[tuple[str, int, int]]]] = defaultdict(lambda: defaultdict(list))
    tx_flags: dict[str, dict[str, tuple[bool, bool]]] = defaultdict(dict)

    wanted = set(panel)
    with _open_text(gtf) as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            chrom, _src, feature, start_s, end_s, _score, _strand, _frame, attrs_s = parts
            if feature not in {"exon", "transcript"}:
                continue
            attrs = _parse_attrs(attrs_s)
            gene_id = (attrs.get("gene_id") or "").split(".")[0]
            if gene_id not in wanted:
                continue
            transcript_id = (attrs.get("transcript_id") or "").split(".")[0]
            if not transcript_id:
                continue
            biotype = attrs.get("transcript_type") or attrs.get("transcript_biotype") or ""
            primary = _is_primary(attrs)
            protein = biotype == "protein_coding"
            if feature == "transcript":
                tx_flags[gene_id][transcript_id] = (primary, protein)
                continue
            start1, end1 = int(start_s), int(end_s)
            tx_exons[gene_id][transcript_id].append((_style_chrom(chrom), start1, end1))
            # Some GTFs omit transcript rows in basic files; record flags from exon attrs too.
            if transcript_id not in tx_flags[gene_id]:
                tx_flags[gene_id][transcript_id] = (primary, protein)
            else:
                prev_p, prev_c = tx_flags[gene_id][transcript_id]
                tx_flags[gene_id][transcript_id] = (prev_p or primary, prev_c or protein)

    gene_spans: dict[str, tuple[str, int, int]] = {}
    gene_exon_sets: dict[str, set[tuple[int, int]]] = {}

    for gene_id in sorted(wanted):
        transcripts = tx_exons.get(gene_id) or {}
        if not transcripts:
            raise SystemExit(f"{gene_id} ({panel[gene_id]}) not found in {gtf}")

        def pick(pred) -> list[str]:
            return [tid for tid, exons in transcripts.items() if pred(tid) and exons]

        primary_pc = pick(lambda tid: tx_flags[gene_id].get(tid, (False, False)) == (True, True))
        any_primary = pick(lambda tid: tx_flags[gene_id].get(tid, (False, False))[0])
        any_pc = pick(lambda tid: tx_flags[gene_id].get(tid, (False, False))[1])
        chosen = primary_pc or any_primary or any_pc or list(transcripts)

        exons: set[tuple[str, int, int]] = set()
        for tid in chosen:
            exons.update(transcripts[tid])
        chrom = next(iter(exons))[0]
        start1 = min(s for _, s, _ in exons)
        end1 = max(e for _, _, e in exons)
        gene_spans[gene_id] = (chrom, start1, end1)
        gene_exon_sets[gene_id] = {(s, e) for _, s, e in exons}

    return gene_spans, gene_exon_sets


def build(args: argparse.Namespace) -> None:
    root = Path(__file__).resolve().parents[2]
    panel = load_panel_genes(Path(args.inheritance_csv))
    gtf_path = Path(args.gtf) if args.gtf else Path(args.gtf_cache)
    if not args.gtf:
        ensure_gtf(gtf_path, args.gtf_url)
    elif not gtf_path.exists():
        raise SystemExit(f"GTF not found: {gtf_path}")

    print(f"parsing {gtf_path} for {len(panel)} genes...", flush=True)
    gene_spans, gene_exons = extract_from_gtf(gtf_path, panel)

    calling_flank = int(args.calling_flank)
    sv_flank = int(args.sv_flank)
    exon_pad = int(args.exon_pad)

    asset_rows: list[tuple[str, int, int, str, str]] = []
    gene_rows: list[tuple[str, int, int, str, str]] = []
    exon_rows: list[tuple[str, int, int, str, str]] = []
    provenance_genes: list[dict] = []

    for gene_id, symbol in sorted(panel.items(), key=lambda kv: kv[1]):
        chrom, body_start, body_end = gene_spans[gene_id]
        body0, body1 = body_start - 1, body_end
        asset_rows.append((chrom, max(0, body0 - calling_flank), body1 + calling_flank, symbol, gene_id))
        gene_rows.append((chrom, max(0, body0 - sv_flank), body1 + sv_flank, gene_id, symbol))
        for start1, end1 in sorted(gene_exons[gene_id]):
            exon_rows.append(
                (
                    chrom,
                    max(0, (start1 - 1) - exon_pad),
                    end1 + exon_pad,
                    gene_id,
                    symbol,
                )
            )
        provenance_genes.append(
            {
                "gene_id": gene_id,
                "symbol": symbol,
                "chrom": chrom,
                "panel_body_1based": [body_start, body_end],
                "n_exons": len(gene_exons[gene_id]),
            }
        )
        print(f"{symbol} {chrom}:{body_start}-{body_end} exons={len(gene_exons[gene_id])}", flush=True)

    def _sort_key(row: tuple) -> tuple:
        chrom = row[0]
        order = chrom.replace("chr", "")
        chrom_rank = (0, int(order)) if order.isdigit() else (1, order)
        return chrom_rank, int(row[1]), int(row[2]), str(row[3]), str(row[4])

    asset_rows.sort(key=_sort_key)
    gene_rows.sort(key=_sort_key)
    exon_rows.sort(key=_sort_key)

    assets_path = Path(args.assets_bed)
    genes_path = Path(args.genes_bed)
    exons_path = Path(args.exons_bed)
    prov_path = Path(args.provenance)

    for path, rows in (
        (assets_path, asset_rows),
        (genes_path, gene_rows),
        (exons_path, exon_rows),
    ):
        path.parent.mkdir(parents=True, exist_ok=True)
        with path.open("w", encoding="utf-8", newline="\n") as fh:
            for row in rows:
                fh.write("\t".join(str(x) for x in row) + "\n")

    provenance = {
        "source": str(gtf_path.name),
        "source_url": args.gtf_url if not args.gtf else None,
        "assembly": "GRCh38",
        "n_genes": len(panel),
        "calling_flank_bp": calling_flank,
        "sv_flank_bp": sv_flank,
        "exon_pad_bp": exon_pad,
        "coordinate_system": "BED 0-based half-open; GTF start/end are 1-based inclusive",
        "contig_naming": "GTF names normalized to chr*; MT -> chrM (Broad GRCh38)",
        "gene_span_policy": (
            "union of exons on GENCODE-primary/Ensembl_canonical protein-coding "
            "transcripts; fallback primary, then protein-coding, then all"
        ),
        "exon_policy": (
            "unique exons from the same transcript set as gene_span_policy; "
            "never label foreign overlapping genes' exons"
        ),
        "content": "gene body plus flank; exons padded; overlapping genes not merged",
        "queried_on": time.strftime("%Y-%m-%d"),
        "genes": provenance_genes,
    }
    prov_path.parent.mkdir(parents=True, exist_ok=True)
    prov_path.write_text(json.dumps(provenance, indent=2) + "\n", encoding="utf-8")

    print(
        f"wrote {assets_path} ({len(asset_rows)}), {genes_path} ({len(gene_rows)}), "
        f"{exons_path} ({len(exon_rows)}), {prov_path}",
        flush=True,
    )


def main(argv: list[str] | None = None) -> None:
    root = Path(__file__).resolve().parents[2]
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument(
        "--inheritance-csv",
        default=str(root / "resources/panel/monopdaus-gene-inheritance.csv"),
    )
    p.add_argument("--gtf", default="", help="Local GENCODE/Ensembl GTF (.gtf or .gtf.gz)")
    p.add_argument(
        "--gtf-cache",
        default=str(root / "resources/gencode/gencode.v44.annotation.gtf.gz"),
    )
    p.add_argument("--gtf-url", default=GENCODE_V44_URL)
    p.add_argument(
        "--assets-bed",
        default=str(root / "workflow/assets/monopdaus-genes.grch38.bed"),
    )
    p.add_argument(
        "--genes-bed",
        default=str(root / "resources/panel/monopdaus-genes.grch38.bed"),
    )
    p.add_argument(
        "--exons-bed",
        default=str(root / "resources/panel/monopdaus-exons.grch38.bed"),
    )
    p.add_argument(
        "--provenance",
        default=str(root / "workflow/assets/panel_bed_provenance.json"),
    )
    p.add_argument("--calling-flank", type=int, default=5000)
    # Keep SV gene-overlap flank aligned with calling so far-neighbour events
    # (e.g. GPR4 under an inflated OPA3 annotation) are not labelled as panel hits.
    p.add_argument("--sv-flank", type=int, default=5000)
    p.add_argument("--exon-pad", type=int, default=100)
    build(p.parse_args(argv))


if __name__ == "__main__":
    main()
