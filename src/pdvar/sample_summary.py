"""Per-sample summary genotypes and scores (not candidate shortlist rows).

Fields:

- APOE ε2/ε3/ε4 from rs429358 + rs7412 (GRCh38). Not on the 133-gene Mendelian
  panel; called from ``summary_vcf``.
- m.2158T>C (SHLP2 p.Lys4Arg, rs41349444) from the mt VCF.
- Haplogroup from Haplogrep (Phylotree rCRS).
- ``panel_prs``: cited PGS000902 weights on panel-gene SNPs.
- Mean MitImpact MLC across mt SNVs, excluding Haplogrep FoundPolys
  (haplogroup-defining / haplotagging variants present in the sample).
"""

from __future__ import annotations

import csv
import re
from pathlib import Path

from .prs import (
    CITATION as PRS_CITATION,
    default_weights_path,
    index_genotypes,
    load_prs_weights,
    score_prs,
)
from .vcf_stream import STREAM_MAX_BYTES, chrom_bare, iter_vcf

# NCBI dbSNP / Manolio APOE ε2/ε3/ε4 definition (rs429358 p.Cys130Arg, rs7412 p.Arg176Cys).
APOE_RS429358 = ("19", 44908684, "T", "C", "rs429358")  # C = ε4
APOE_RS7412 = ("19", 44908822, "C", "T", "rs7412")  # T = ε2
APOE_CITATION = "rs429358 (p.Cys130Arg) + rs7412 (p.Arg176Cys); NCBI dbSNP / standard APOE ε2/ε3/ε4"

# Hudson G et al. Neurology 2013;80:2042-2048. Kim SJ et al. Mol Psychiatry 2024 (SHLP2 K4R).
MT_2158 = ("M", 2158, "T", "C", "rs41349444")
MT_2158_CITATION = (
    "m.2158T>C (SHLP2 p.Lys4Arg, rs41349444); Hudson G et al. Neurology 2013;80:2042-2048; "
    "Kim SJ et al. Mol Psychiatry 2024. doi:10.1038/s41380-023-02344-0"
)

HAPLOGREP_CITATION = (
    "Haplogrep 3 (Schönherr et al. Nucleic Acids Res 2023) with Phylotree rCRS; "
    "Weissensteiner et al. Nucleic Acids Res 2016"
)
MLC_CITATION = (
    "Lake NJ et al. Nature 2024;625:593-602. doi:10.1038/s41586-024-08048-x "
    "(MLC via MitImpact); haplotagging exclusion = Haplogrep FoundPolys"
)

_POLY_POS = re.compile(r"(?<![0-9])(\d{1,5})(?![0-9])")
_BLANK = {"", ".", "na", "nan", "none", "<na>"}


def _fmt(record: dict, sample_id: str) -> dict:
    samples = record.get("samples") or {}
    if sample_id in samples:
        return samples[sample_id]
    if len(samples) == 1:
        return next(iter(samples.values()))
    return {}


def _gt_bases(gt: str, ref: str, alt: str) -> list[str] | None:
    if not gt or gt in {".", "./.", ".|."}:
        return None
    alleles = gt.replace("|", "/").split("/")
    alts = [item.strip().upper() for item in str(alt).split(",") if item.strip()]
    try:
        bases = []
        for allele in alleles:
            token = str(allele).strip()
            if token in {"", "."}:
                return None
            if token == "0":
                bases.append(str(ref).upper())
            else:
                bases.append(alts[int(token) - 1])
        return bases
    except (ValueError, IndexError):
        return None


def _lookup(index: dict, chrom: str, pos: int) -> dict | None:
    return index.get((chrom_bare(chrom), int(pos)))


def _count_allele(bases: list[str] | None, allele: str) -> int | None:
    if bases is None:
        return None
    want = str(allele).upper()
    return sum(base == want for base in bases)


def apoe_genotype(index: dict[tuple[str, int], dict], sample_id: str) -> dict:
    """Unphased ε2/ε3/ε4 from rs429358 + rs7412."""
    out = {
        "apoe_genotype": None,
        "apoe_rs429358_gt": None,
        "apoe_rs7412_gt": None,
        "apoe_status": "not_assessed",
        "apoe_citation": APOE_CITATION,
    }
    rec_e4 = _lookup(index, APOE_RS429358[0], APOE_RS429358[1])
    rec_e2 = _lookup(index, APOE_RS7412[0], APOE_RS7412[1])
    if rec_e4 is None and rec_e2 is None:
        return out

    def _site(rec, chrom, pos, ref, alt):
        if rec is None:
            # All-sites summary VCF should emit 0/0; absence is no coverage.
            return None, None
        fmt = _fmt(rec, sample_id)
        bases = _gt_bases(str(fmt.get("GT") or ""), str(rec["ref"]), str(rec["alt"]))
        gt = "/".join(sorted(bases)) if bases else None
        return bases, gt

    bases_e4, gt_e4 = _site(rec_e4, *APOE_RS429358[:4])
    bases_e2, gt_e2 = _site(rec_e2, *APOE_RS7412[:4])
    out["apoe_rs429358_gt"] = gt_e4
    out["apoe_rs7412_gt"] = gt_e2
    n_e4 = _count_allele(bases_e4, "C")  # rs429358 C
    n_e2 = _count_allele(bases_e2, "T")  # rs7412 T
    if n_e4 is None or n_e2 is None:
        out["apoe_status"] = "partial"
        return out
    key = (int(n_e4), int(n_e2))
    mapping = {
        (0, 0): "e3/e3",
        (1, 0): "e3/e4",
        (2, 0): "e4/e4",
        (0, 1): "e2/e3",
        (0, 2): "e2/e2",
        (1, 1): "e2/e4_or_e1/e3",
        (1, 2): "e1/e2",
        (2, 1): "e1/e4",
        (2, 2): "e1/e1",
    }
    out["apoe_genotype"] = mapping.get(key)
    out["apoe_status"] = "called" if out["apoe_genotype"] else "not_assessed"
    return out


def mt_2158_genotype(mtdna_vcf: str | Path | None, sample_id: str, mtdna_called: bool) -> dict:
    """m.2158T>C from the MitochondriaPipeline VCF. Absent ALT after a successful call is T."""
    out = {
        "mt_2158_genotype": None,
        "mt_2158_heteroplasmy": None,
        "mt_2158_status": "not_assessed",
        "mt_2158_citation": MT_2158_CITATION,
    }
    if not mtdna_vcf:
        return out
    path = Path(mtdna_vcf)
    if not path.exists():
        return out
    chrom, pos, ref, alt, _rsid = MT_2158
    hit = None
    for record in iter_vcf(path, allow_full_scan=True):
        if chrom_bare(record["chrom"]) != chrom or int(record["pos"]) != pos:
            continue
        if str(record["ref"]).upper() != ref:
            continue
        if alt.upper() in {item.strip().upper() for item in str(record["alt"]).split(",")}:
            hit = record
            break
    if hit is None:
        if mtdna_called:
            out["mt_2158_genotype"] = "T"
            out["mt_2158_heteroplasmy"] = 0.0
            out["mt_2158_status"] = "called"
        return out
    fmt = _fmt(hit, sample_id)
    af_raw = fmt.get("AF") or (hit.get("info") or {}).get("AF")
    try:
        af = float(str(af_raw).split(",")[0]) if af_raw not in (None, "", ".") else None
    except (TypeError, ValueError):
        af = None
    if af is None:
        bases = _gt_bases(str(fmt.get("GT") or ""), str(hit["ref"]), str(hit["alt"]))
        if bases and any(base == alt for base in bases):
            af = 1.0 if all(base == alt for base in bases) else 0.5
        else:
            af = 0.0
    out["mt_2158_heteroplasmy"] = af
    if af >= 0.95:
        out["mt_2158_genotype"] = "C"
    elif af >= 0.01:
        out["mt_2158_genotype"] = "T/C"
    else:
        out["mt_2158_genotype"] = "T"
    out["mt_2158_status"] = "called"
    return out


def _norm_header(name: str) -> str:
    return re.sub(r"[^a-z0-9]", "", str(name).strip().lower())


def parse_haplogrep(path: str | Path | None) -> dict:
    """Read Haplogrep classify TSV (extend-report). FoundPolys = haplotagging sites."""
    out = {
        "haplogroup": None,
        "haplogroup_quality": None,
        "haplogroup_status": "not_assessed",
        "haplogroup_citation": HAPLOGREP_CITATION,
        "haplotag_positions": [],
        "haplotag_source": None,
    }
    if not path:
        return out
    file_path = Path(path)
    if not file_path.exists():
        return out
    text = file_path.read_text(encoding="utf-8", errors="replace")
    lines = [line for line in text.splitlines() if line.strip() and not line.startswith("##")]
    if not lines:
        return out
    reader = csv.DictReader(lines, delimiter="\t")
    rows = list(reader)
    if not rows:
        return out
    row = rows[0]
    lookup = {_norm_header(key): key for key in row}
    hg_key = next((lookup[name] for name in ("haplogroup", "hg", "hgmajor") if name in lookup), None)
    q_key = next((lookup[name] for name in ("quality", "rankquality", "score") if name in lookup), None)
    found_key = next(
        (lookup[name] for name in ("foundpolys", "found_polys", "foundpolymorphisms", "found") if name in lookup),
        None,
    )
    haplogroup = str(row.get(hg_key) or "").strip().strip('"') if hg_key else ""
    if haplogroup.lower() in _BLANK:
        haplogroup = ""
    quality = None
    if q_key:
        raw = str(row.get(q_key) or "").strip().strip('"').rstrip("%")
        try:
            quality = float(raw)
        except (TypeError, ValueError):
            quality = None
    positions: list[int] = []
    if found_key:
        blob = str(row.get(found_key) or "")
        for match in _POLY_POS.finditer(blob):
            pos = int(match.group(1))
            if 1 <= pos <= 16569:
                positions.append(pos)
    out["haplogroup"] = haplogroup or None
    out["haplogroup_quality"] = quality
    out["haplogroup_status"] = "called" if haplogroup else "not_assessed"
    out["haplotag_positions"] = sorted(set(positions))
    out["haplotag_source"] = "haplogrep_found_polys" if found_key else None
    return out


def mean_mlc(
    mtdna_vcf: str | Path | None,
    haplotag_positions: list[int] | None,
    haplotag_applied: bool,
) -> dict:
    """Mean MitImpact MLC on mt SNVs, excluding Haplogrep FoundPolys when available."""
    out = {
        "mean_mlc": None,
        "mlc_n_variants": 0,
        "mlc_n_excluded_haplotag": 0,
        "mlc_n_missing_score": 0,
        "mlc_haplotag_exclusion": "not_applied",
        "mlc_status": "not_assessed",
        "mlc_citation": MLC_CITATION,
    }
    if not mtdna_vcf:
        return out
    path = Path(mtdna_vcf)
    if not path.exists():
        return out
    exclude = set(haplotag_positions or [])
    applied = bool(haplotag_applied and exclude)
    scores: list[float] = []
    n_excluded = 0
    n_missing = 0
    n_snv = 0
    for record in iter_vcf(path, allow_full_scan=True):
        if chrom_bare(record["chrom"]) not in {"M", "MT"}:
            continue
        ref = str(record["ref"]).upper()
        alt = str(record["alt"]).split(",")[0].upper()
        if len(ref) != 1 or len(alt) != 1:
            continue
        n_snv += 1
        pos = int(record["pos"])
        if applied and pos in exclude:
            n_excluded += 1
            continue
        raw = (record.get("info") or {}).get("MI_MLC_score")
        if raw in (None, "", "."):
            n_missing += 1
            continue
        try:
            scores.append(float(raw))
        except (TypeError, ValueError):
            n_missing += 1
    out["mlc_n_variants"] = n_snv
    out["mlc_n_excluded_haplotag"] = n_excluded
    out["mlc_n_missing_score"] = n_missing
    if applied:
        out["mlc_haplotag_exclusion"] = "haplogrep_found_polys"
    elif haplotag_applied:
        out["mlc_haplotag_exclusion"] = "empty_found_polys"
    else:
        out["mlc_haplotag_exclusion"] = "not_applied"
    if scores:
        out["mean_mlc"] = round(sum(scores) / len(scores), 6)
        out["mlc_status"] = "called"
    elif n_snv:
        out["mlc_status"] = "not_assessed"
    return out


def _nuclear_genotypes(sample: dict) -> dict[tuple[str, int], dict]:
    """Index APOE/PRS sites from the all-sites summary VCF.

    Do not scan ``vep_vcf`` / ``small_vcf`` when the summary VCF exists: those
    callsets are variant-only and, on ``wgs_small_variants``, too large to read
    end to end. Fall back to a panel-sized VCF only when it is under
    ``STREAM_MAX_BYTES``.
    """
    genotypes = index_genotypes(sample.get("summary_vcf"))
    if genotypes:
        return genotypes
    for key in ("vep_vcf", "small_vcf"):
        path = sample.get(key)
        if not path:
            continue
        file_path = Path(str(path))
        if file_path.exists() and file_path.stat().st_size <= STREAM_MAX_BYTES:
            genotypes.update(index_genotypes(file_path))
    return genotypes


def build_sample_summary(sample: dict, cfg: dict | None = None) -> dict:
    """Compute the five sample-level fields from Assess inputs."""
    sample_id = str(sample.get("sample_id") or "")
    genotypes = _nuclear_genotypes(sample)
    apoe = apoe_genotype(genotypes, sample_id)

    mtdna_path = sample.get("mtdna_vcf")
    mtdna_called = str(sample.get("mtdna_annotation") or "").lower() in {"called", "prior"}
    # A present mt VCF still counts as called for the 2158 reference-T fallback.
    if mtdna_path and Path(str(mtdna_path)).exists():
        mtdna_called = True
    mt2158 = mt_2158_genotype(mtdna_path, sample_id, mtdna_called)

    haplo = parse_haplogrep(sample.get("haplogrep_tsv"))
    mlc = mean_mlc(
        mtdna_path,
        haplo.get("haplotag_positions") or [],
        haplotag_applied=haplo.get("haplogroup_status") == "called" and haplo.get("haplotag_source") == "haplogrep_found_polys",
    )

    weights_path = sample.get("prs_weights") or default_weights_path(cfg)
    weights = load_prs_weights(weights_path, panel_only=True)
    prs = score_prs(weights, genotypes, sample_id)

    summary = {
        "sample_id": sample_id,
        **apoe,
        **mt2158,
        "haplogroup": haplo.get("haplogroup"),
        "haplogroup_quality": haplo.get("haplogroup_quality"),
        "haplogroup_status": haplo.get("haplogroup_status"),
        "haplogroup_citation": haplo.get("haplogroup_citation"),
        "haplotag_n": len(haplo.get("haplotag_positions") or []),
        **prs,
        **mlc,
    }
    if prs.get("prs_citation") is None:
        summary["prs_citation"] = PRS_CITATION
    return summary


SUMMARY_TSV_COLUMNS = [
    "sample_id",
    "apoe_genotype",
    "apoe_rs429358_gt",
    "apoe_rs7412_gt",
    "apoe_status",
    "mt_2158_genotype",
    "mt_2158_heteroplasmy",
    "mt_2158_status",
    "haplogroup",
    "haplogroup_quality",
    "haplogroup_status",
    "panel_prs",
    "prs_sites_used",
    "prs_sites_missing",
    "prs_sites_total",
    "prs_status",
    "mean_mlc",
    "mlc_n_variants",
    "mlc_n_excluded_haplotag",
    "mlc_haplotag_exclusion",
    "mlc_status",
]


def summary_to_row(summary: dict) -> dict:
    return {key: summary.get(key) for key in SUMMARY_TSV_COLUMNS}
