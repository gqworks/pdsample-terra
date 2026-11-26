#!/usr/bin/env python3
"""Paralog-aware GBA1 exon 9–11 hotspot force-call.

Reads overlapping GBA1 or GBAP1 exon 9–11 windows are scored against the
GRCh38 gene template and a *patched* GBAP1 template (GRCh38 GBAP1 sites that
erroneously carry GBA-like bases are restored to the major GBAP1 alleles).
Gene-better fragments are piled at pathogenic hotspots (L444P / A495P / D448H).
Emits a VCF (including hom-ref QC rows) and JSON summary. Requires pysam.
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path


_COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


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


@dataclass(frozen=True)
class Patch:
    chrom: str
    pos: int
    grch38_allele: str
    patched_allele: str
    note: str


def bare_chrom(chrom: str) -> str:
    text = str(chrom or "").strip()
    return text[3:] if text.lower().startswith("chr") else text


def revcomp(seq: str) -> str:
    return seq.translate(_COMPLEMENT)[::-1]


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


def load_patches(path: Path) -> list[Patch]:
    rows: list[Patch] = []
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(
            (line for line in handle if line.strip() and not line.lstrip().startswith("#")),
            delimiter="\t",
        )
        for row in reader:
            rows.append(
                Patch(
                    chrom=str(row["chrom"]).strip(),
                    pos=int(row["pos"]),
                    grch38_allele=str(row["grch38_allele"]).strip().upper(),
                    patched_allele=str(row["patched_allele"]).strip().upper(),
                    note=str(row.get("note") or "").strip(),
                )
            )
    if not rows:
        raise ValueError(f"{path} has no patch rows")
    return rows


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


def _fetch_sequence(fasta, chrom: str, start: int, end: int) -> tuple[str, str, int, int]:
    """Return (seq, contigname, start, end) for 1-based inclusive coords."""
    names = list(fasta.references)
    candidates = [chrom]
    bare = bare_chrom(chrom)
    if chrom.lower().startswith("chr"):
        candidates.append(bare)
    else:
        candidates.append(f"chr{chrom}")
    for name in candidates:
        if name in fasta.references or name in names:
            # pysam fetch is 0-based half-open
            seq = fasta.fetch(name, start - 1, end).upper()
            return seq, name, start, end
    raise KeyError(f"chromosome {chrom!r} not in reference (tried {candidates})")


def apply_patches(
    seq: str,
    *,
    window_start: int,
    patches: list[Patch],
    chrom: str,
) -> tuple[str, list[dict]]:
    """Apply GBAP1 GRCh38→major patches inside a window sequence."""
    chars = list(seq)
    applied: list[dict] = []
    chrom_key = bare_chrom(chrom)
    for patch in patches:
        if bare_chrom(patch.chrom) != chrom_key:
            continue
        if patch.pos < window_start or patch.pos > window_start + len(chars) - 1:
            continue
        idx = patch.pos - window_start
        observed = chars[idx]
        if observed not in {patch.grch38_allele, patch.patched_allele}:
            raise ValueError(
                f"GBAP1 patch at {patch.chrom}:{patch.pos}: expected "
                f"{patch.grch38_allele} or {patch.patched_allele}, found {observed}"
            )
        chars[idx] = patch.patched_allele
        applied.append(
            {
                "pos": patch.pos,
                "from": observed,
                "to": patch.patched_allele,
                "note": patch.note,
            }
        )
    return "".join(chars), applied


def best_ungapped_matches(query: str, template: str) -> int:
    """Maximum ACGT matches for an ungapped placement of query on template."""
    q = (query or "").upper()
    t = (template or "").upper()
    if not q or not t:
        return 0
    if len(q) > len(t):
        # Rare (very long read); slide template over query instead.
        return best_ungapped_matches(t, q)
    best = 0
    qlen = len(q)
    for i in range(len(t) - qlen + 1):
        window = t[i : i + qlen]
        matches = 0
        for a, b in zip(q, window):
            if a in "ACGT" and b in "ACGT" and a == b:
                matches += 1
        if matches > best:
            best = matches
    return best


def score_read_vs_templates(query: str, gene_seq: str, pseudo_seq: str) -> tuple[int, int]:
    """Return (gene_matches, pseudo_matches) using the best strand for each."""
    if not query:
        return 0, 0
    rc = revcomp(query)
    gene = max(best_ungapped_matches(query, gene_seq), best_ungapped_matches(rc, gene_seq))
    pseudo = max(best_ungapped_matches(query, pseudo_seq), best_ungapped_matches(rc, pseudo_seq))
    return gene, pseudo


def forcecall(
    bam_path: Path,
    *,
    reference: Path,
    discriminators: list[Discriminator],
    hotspots: list[Hotspot],
    patches: list[Patch],
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

    with pysam.FastaFile(str(reference)) as fasta:
        gene_seq, gene_contig, gene_start, gene_end = _fetch_sequence(
            fasta, gene_win[0], gene_win[1], gene_win[2]
        )
        pseudo_raw, pseudo_contig, pseudo_start, pseudo_end = _fetch_sequence(
            fasta, pseudo_win[0], pseudo_win[1], pseudo_win[2]
        )
    pseudo_seq, patches_applied = apply_patches(
        pseudo_raw,
        window_start=pseudo_start,
        patches=patches,
        chrom=pseudo_win[0],
    )

    # query_name -> (best observed bases, best read sequence seen)
    by_qname_obs: dict[str, dict[tuple[str, int], str]] = defaultdict(dict)
    by_qname_seq: dict[str, str] = {}

    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        for label, (chrom, start, end) in (("gene", gene_win), ("pseudo", pseudo_win)):
            try:
                iterator = bam.fetch(chrom, start - 1, end)
            except ValueError:
                alt = bare_chrom(chrom) if chrom.lower().startswith("chr") else f"chr{chrom}"
                iterator = bam.fetch(alt, start - 1, end)
            for read in iterator:
                if read.is_unmapped or read.is_duplicate or read.is_qcfail:
                    continue
                if read.mapping_quality < mapq_floor:
                    continue
                seq = read.query_sequence
                if not seq:
                    continue
                observed: dict[tuple[str, int], str] = {}
                positions = (
                    gene_positions | hotspot_gene_pos if label == "gene" else pseudo_positions | hotspot_pseudo_pos
                )
                for pos in positions:
                    base = _read_base_at(read, chrom, pos)
                    if base:
                        observed[(bare_chrom(chrom), pos)] = base
                qname = read.query_name or f"anon:{read.reference_start}"
                if observed:
                    by_qname_obs[qname].update(observed)
                # Keep the longest sequence for the fragment (mates may differ).
                prev = by_qname_seq.get(qname)
                if prev is None or len(seq) > len(prev):
                    by_qname_seq[qname] = seq

    gene_assigned_obs: list[dict[tuple[str, int], str]] = []
    assignment_counts = Counter()
    for qname, seq in by_qname_seq.items():
        gene_s, pseudo_s = score_read_vs_templates(seq, gene_seq, pseudo_seq)
        locus = assign_locus(gene_s, pseudo_s)
        assignment_counts[locus] += 1
        if locus == "gene":
            gene_assigned_obs.append(by_qname_obs.get(qname, {}))

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
        "assignment_mode": "patched_grch38_template",
        "fragments_scored": int(sum(assignment_counts.values())),
        "assignment_counts": dict(assignment_counts),
        "gene_assigned_fragments": int(assignment_counts.get("gene", 0)),
        "templates": {
            "gene": {"contig": gene_contig, "start": gene_start, "end": gene_end, "length": len(gene_seq)},
            "pseudo_patched": {
                "contig": pseudo_contig,
                "start": pseudo_start,
                "end": pseudo_end,
                "length": len(pseudo_seq),
                "patches_applied": patches_applied,
            },
        },
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
    parser.add_argument("--reference", required=True, type=Path, help="GRCh38 FASTA (indexed)")
    parser.add_argument("--discriminators", required=True, type=Path)
    parser.add_argument("--hotspots", required=True, type=Path)
    parser.add_argument("--patches", required=True, type=Path, help="GBAP1 GRCh38→major allele patches")
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
    patches = load_patches(args.patches)
    sites, summary = forcecall(
        args.bam,
        reference=args.reference,
        discriminators=discriminators,
        hotspots=hotspots,
        patches=patches,
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
        f"mode=patched_grch38_template; "
        f"gene_assigned_fragments={summary['gene_assigned_fragments']}; "
        f"carriers={','.join(carriers) or 'none'}",
        file=sys.stderr,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
