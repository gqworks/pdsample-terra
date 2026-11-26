#!/usr/bin/env python3
"""Paralog-aware GBA1 exon 9–11 hotspot force-call.

Reads mapping to GBA1 or GBAP1 are reassigned using Gauchian discriminating
SNPs, then only gene-assigned fragments are piled at pathogenic hotspots
(L444P / A495P / D448H). Emits a VCF (including hom-ref QC rows) and JSON
summary. Requires pysam.
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class Discriminator:
    chrom: str
    gene_pos: int
    gene_allele: str
    pseudo_pos: int
    pseudo_allele: str
    label: str


@dataclass(frozen=True)
class Hotspot:
    token: str
    cdna: str
    protein: str
    chrom: str
    pos: int
    ref: str
    alt: str
    pseudo_pos: int
    pd_class: str


def bare_chrom(chrom: str) -> str:
    text = str(chrom or "").strip()
    return text[3:] if text.lower().startswith("chr") else text


def load_discriminators(path: Path) -> list[Discriminator]:
    rows: list[Discriminator] = []
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(
            (line for line in handle if line.strip() and not line.lstrip().startswith("#")),
            delimiter="\t",
        )
        for row in reader:
            rows.append(
                Discriminator(
                    chrom=str(row["chrom"]).strip(),
                    gene_pos=int(row["gene_pos"]),
                    gene_allele=str(row["gene_allele"]).strip().upper(),
                    pseudo_pos=int(row["pseudo_pos"]),
                    pseudo_allele=str(row["pseudo_allele"]).strip().upper(),
                    label=str(row.get("label") or "").strip(),
                )
            )
    if not rows:
        raise ValueError(f"{path} has no discriminator rows")
    return rows


def load_hotspots(path: Path) -> list[Hotspot]:
    rows: list[Hotspot] = []
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(
            (line for line in handle if line.strip() and not line.lstrip().startswith("#")),
            delimiter="\t",
        )
        for row in reader:
            rows.append(
                Hotspot(
                    token=str(row["token"]).strip(),
                    cdna=str(row["cdna"]).strip(),
                    protein=str(row.get("protein") or "").strip(),
                    chrom=str(row["chrom"]).strip(),
                    pos=int(row["pos"]),
                    ref=str(row["ref"]).strip().upper(),
                    alt=str(row["alt"]).strip().upper(),
                    pseudo_pos=int(row["pseudo_pos"]),
                    pd_class=str(row.get("pd_class") or "").strip().lower(),
                )
            )
    if not rows:
        raise ValueError(f"{path} has no hotspot rows")
    return rows


def score_assignment(
    observed: dict[tuple[str, int], str],
    discriminators: list[Discriminator],
    *,
    locus: str,
) -> tuple[int, int]:
    """Return (gene_score, pseudo_score) from bases seen on one fragment."""
    gene = 0
    pseudo = 0
    for disc in discriminators:
        if locus == "gene":
            key = (bare_chrom(disc.chrom), disc.gene_pos)
            gene_base, pseudo_base = disc.gene_allele, disc.pseudo_allele
        else:
            key = (bare_chrom(disc.chrom), disc.pseudo_pos)
            gene_base, pseudo_base = disc.gene_allele, disc.pseudo_allele
        base = observed.get(key)
        if not base:
            continue
        base = base.upper()
        if base == gene_base:
            gene += 1
        elif base == pseudo_base:
            pseudo += 1
    return gene, pseudo


def assign_locus(gene_score: int, pseudo_score: int) -> str:
    if gene_score == 0 and pseudo_score == 0:
        return "unassigned"
    if gene_score > pseudo_score:
        return "gene"
    if pseudo_score > gene_score:
        return "pseudo"
    return "unassigned"


def genotype_from_counts(
    counts: Counter[str],
    *,
    ref: str,
    alt: str,
    min_alt_reads: int,
    min_alt_fraction: float,
    min_gene_dp: int,
) -> tuple[str, int, int, int]:
    """Return GT, DP, AD_ref, AD_alt from gene-assigned base counts."""
    ad_ref = int(counts.get(ref, 0))
    ad_alt = int(counts.get(alt, 0))
    dp = sum(counts.values())
    if dp < min_gene_dp:
        return "./.", dp, ad_ref, ad_alt
    if ad_alt >= min_alt_reads and (ad_alt / dp) >= min_alt_fraction:
        if ad_ref == 0 and ad_alt >= min_alt_reads:
            return "1/1", dp, ad_ref, ad_alt
        return "0/1", dp, ad_ref, ad_alt
    return "0/0", dp, ad_ref, ad_alt


def _read_base_at(read, chrom: str, pos: int) -> str | None:
    """Return the read base aligned to 1-based reference pos, if any."""
    try:
        pairs = read.get_aligned_pairs(matches_only=True)
    except Exception:
        return None
    for query_pos, ref_pos in pairs:
        if ref_pos is None or query_pos is None:
            continue
        if int(ref_pos) + 1 == pos and bare_chrom(read.reference_name or "") == bare_chrom(chrom):
            seq = read.query_sequence
            if not seq or query_pos >= len(seq):
                return None
            return seq[query_pos].upper()
    return None


def _fragment_windows(discriminators: list[Discriminator], pad: int) -> tuple[tuple[str, int, int], tuple[str, int, int]]:
    chrom = discriminators[0].chrom
    gene_start = min(d.gene_pos for d in discriminators) - pad
    gene_end = max(d.gene_pos for d in discriminators) + pad
    pseudo_start = min(d.pseudo_pos for d in discriminators) - pad
    pseudo_end = max(d.pseudo_pos for d in discriminators) + pad
    return (chrom, max(1, gene_start), gene_end), (chrom, max(1, pseudo_start), pseudo_end)


def forcecall(
    bam_path: Path,
    *,
    discriminators: list[Discriminator],
    hotspots: list[Hotspot],
    sample_id: str,
    min_alt_reads: int = 3,
    min_alt_fraction: float = 0.15,
    min_gene_dp: int = 8,
    mapq_floor: int = 0,
    pad: int = 150,
) -> tuple[list[dict], dict]:
    """Pile gene-assigned fragments at hotspots. Returns (site dicts, summary)."""
    import pysam

    gene_win, pseudo_win = _fragment_windows(discriminators, pad)
    gene_positions = {d.gene_pos for d in discriminators}
    pseudo_positions = {d.pseudo_pos for d in discriminators}
    hotspot_gene_pos = {h.pos for h in hotspots}
    hotspot_pseudo_pos = {h.pseudo_pos for h in hotspots}

    # query_name -> list of (locus_hint, observed bases, alignment locus)
    by_qname: dict[str, list[tuple[str, dict[tuple[str, int], str]]]] = defaultdict(list)

    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        for label, (chrom, start, end) in (("gene", gene_win), ("pseudo", pseudo_win)):
            # pysam fetch is 0-based half-open
            try:
                iterator = bam.fetch(chrom, start - 1, end)
            except ValueError:
                # try without chr prefix / with chr
                alt = bare_chrom(chrom) if chrom.lower().startswith("chr") else f"chr{chrom}"
                iterator = bam.fetch(alt, start - 1, end)
            for read in iterator:
                if read.is_unmapped or read.is_duplicate or read.is_qcfail:
                    continue
                if read.mapping_quality < mapq_floor:
                    continue
                observed: dict[tuple[str, int], str] = {}
                positions = gene_positions | hotspot_gene_pos if label == "gene" else pseudo_positions | hotspot_pseudo_pos
                for pos in positions:
                    base = _read_base_at(read, chrom, pos)
                    if base:
                        observed[(bare_chrom(chrom), pos)] = base
                if not observed:
                    continue
                qname = read.query_name or f"anon:{read.reference_start}"
                by_qname[qname].append((label, observed))

    # Merge observations per fragment; score assignment.
    gene_assigned_obs: list[dict[tuple[str, int], str]] = []
    assignment_counts = Counter()
    for _qname, parts in by_qname.items():
        merged: dict[tuple[str, int], str] = {}
        locus_hints = []
        for hint, obs in parts:
            locus_hints.append(hint)
            merged.update(obs)
        # Prefer scoring in the coordinate system with more discriminator hits.
        g_gene, g_pseudo = score_assignment(merged, discriminators, locus="gene")
        p_gene, p_pseudo = score_assignment(merged, discriminators, locus="pseudo")
        # Use the stronger of the two views (gene coords vs pseudo coords).
        if (g_gene + g_pseudo) >= (p_gene + p_pseudo):
            gene_s, pseudo_s = g_gene, g_pseudo
        else:
            gene_s, pseudo_s = p_gene, p_pseudo
        locus = assign_locus(gene_s, pseudo_s)
        assignment_counts[locus] += 1
        if locus == "gene":
            gene_assigned_obs.append(merged)

    sites: list[dict] = []
    for hotspot in hotspots:
        counts: Counter[str] = Counter()
        chrom_key = bare_chrom(hotspot.chrom)
        for obs in gene_assigned_obs:
            base = obs.get((chrom_key, hotspot.pos))
            if not base:
                # Gene-assigned fragment may only cover the GBAP1 image of the site.
                base = obs.get((chrom_key, hotspot.pseudo_pos))
            if base:
                counts[base] += 1
        gt, dp, ad_ref, ad_alt = genotype_from_counts(
            counts,
            ref=hotspot.ref,
            alt=hotspot.alt,
            min_alt_reads=min_alt_reads,
            min_alt_fraction=min_alt_fraction,
            min_gene_dp=min_gene_dp,
        )
        sites.append(
            {
                "chrom": hotspot.chrom if hotspot.chrom.lower().startswith("chr") else f"chr{hotspot.chrom}",
                "pos": hotspot.pos,
                "ref": hotspot.ref,
                "alt": hotspot.alt,
                "token": hotspot.token,
                "cdna": hotspot.cdna,
                "protein": hotspot.protein,
                "pd_class": hotspot.pd_class,
                "gt": gt,
                "dp": dp,
                "ad_ref": ad_ref,
                "ad_alt": ad_alt,
                "base_counts": dict(counts),
            }
        )

    summary = {
        "sample_id": sample_id,
        "fragments_scored": int(sum(assignment_counts.values())),
        "assignment_counts": dict(assignment_counts),
        "gene_assigned_fragments": int(assignment_counts.get("gene", 0)),
        "sites": [
            {
                "token": s["token"],
                "cdna": s["cdna"],
                "gt": s["gt"],
                "dp": s["dp"],
                "ad_ref": s["ad_ref"],
                "ad_alt": s["ad_alt"],
            }
            for s in sites
        ],
    }
    return sites, summary


def write_vcf(path: Path, sample_id: str, sites: list[dict]) -> None:
    lines = [
        "##fileformat=VCFv4.2",
        "##source=gba1_paralog_forcecall",
        '##INFO=<ID=PDVAR_GBA1_PARALOG,Number=0,Type=Flag,Description="Gene-assigned GBA1 hotspot from paralog-aware force-call">',
        '##INFO=<ID=TOKEN,Number=1,Type=String,Description="Gauchian-style protein token">',
        '##INFO=<ID=HGVSC,Number=1,Type=String,Description="HGVS c. notation">',
        '##INFO=<ID=HGVSP,Number=1,Type=String,Description="HGVS p. notation">',
        '##INFO=<ID=GBA1_PD_CLASS,Number=1,Type=String,Description="GBA1-PD Browser class">',
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
        '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Gene-assigned depth">',
        '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Gene-assigned ref,alt depths">',
        f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample_id}",
    ]
    for site in sorted(sites, key=lambda row: (bare_chrom(row["chrom"]), int(row["pos"]))):
        info = (
            f"PDVAR_GBA1_PARALOG;TOKEN={site['token']};HGVSC={site['cdna']};"
            f"HGVSP={site['protein']};GBA1_PD_CLASS={site['pd_class'] or '.'}"
        )
        filt = "PASS" if site["gt"] not in {"./.", "."} else "LowGeneDP"
        lines.append(
            "\t".join(
                [
                    str(site["chrom"]),
                    str(site["pos"]),
                    site["token"] or ".",
                    site["ref"],
                    site["alt"],
                    ".",
                    filt,
                    info,
                    "GT:DP:AD",
                    f"{site['gt']}:{site['dp']}:{site['ad_ref']},{site['ad_alt']}",
                ]
            )
        )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bam", required=True, type=Path)
    parser.add_argument("--discriminators", required=True, type=Path)
    parser.add_argument("--hotspots", required=True, type=Path)
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--out-vcf", required=True, type=Path)
    parser.add_argument("--out-json", required=True, type=Path)
    parser.add_argument("--min-alt-reads", type=int, default=3)
    parser.add_argument("--min-alt-fraction", type=float, default=0.15)
    parser.add_argument("--min-gene-dp", type=int, default=8)
    parser.add_argument("--mapq-floor", type=int, default=0)
    args = parser.parse_args(argv)

    discriminators = load_discriminators(args.discriminators)
    hotspots = load_hotspots(args.hotspots)
    sites, summary = forcecall(
        args.bam,
        discriminators=discriminators,
        hotspots=hotspots,
        sample_id=args.sample_id,
        min_alt_reads=args.min_alt_reads,
        min_alt_fraction=args.min_alt_fraction,
        min_gene_dp=args.min_gene_dp,
        mapq_floor=args.mapq_floor,
    )
    write_vcf(args.out_vcf, args.sample_id, sites)
    args.out_json.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    carriers = [s["token"] for s in sites if s["gt"] not in {"0/0", "./.", "."}]
    print(
        f"gba1_paralog_forcecall: {len(sites)} site(s); "
        f"gene_assigned_fragments={summary['gene_assigned_fragments']}; "
        f"carriers={','.join(carriers) or 'none'}",
        file=sys.stderr,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
