"""Ingest one sample's ExpansionHunter, Gauchian, and mtDNA outputs.

No STR pathogenicity cutoff is applied unless config sets both
``pathogenic_min_repeats`` and a citation.
"""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

from .normalize import derive_zygosity
from .vcf_stream import chrom_bare, iter_vcf, read_samples

MODULES = {
    "expansionhunter": {
        "not_assessed": "STR not assessed",
        "carrier": "STR carrier",
        "absent": "no STR",
    },
    "gauchian": {
        "not_assessed": "GBA1 not assessed",
        "carrier": "GBA1 carrier",
        "absent": "no GBA1 variant",
    },
    "mtdna": {
        "not_assessed": "mtDNA not assessed",
        "carrier": "mtDNA carrier",
        "absent": "no mtDNA variant",
    },
}

VARIANT_COLUMNS = [
    "sample_id", "family_id", "chrom", "pos", "end", "ref", "alt", "variant_id",
    "gene", "consequence", "genotype", "zygosity", "gq", "allele_balance", "depth",
    "cohort_af", "gnomad_af", "gnomad_popmax_af", "cadd", "revel", "spliceai",
    "clinvar", "clinvar_class", "variant_class", "hit_source", "heteroplasmy",
    "repeat_count", "str_locus", "str_expanded", "str_citation", "copy_number",
    "dosage_mechanism", "dosage_event", "scored_dosage",
]


def empty_variant_table() -> pd.DataFrame:
    return pd.DataFrame(columns=VARIANT_COLUMNS)


def _blank_row(**values: object) -> dict:
    row = {column: pd.NA for column in VARIANT_COLUMNS}
    row.update(values)
    return row


def _float_or_none(value: object) -> float | None:
    if value is None or value is True or value in {"", "."}:
        return None
    try:
        return float(str(value).split(",")[0])
    except ValueError:
        return None


def _rep_counts(value: object) -> list[float]:
    if value is None or value is True or value in {"", "."}:
        return []
    counts = []
    for part in str(value).replace("|", "/").split("/"):
        try:
            counts.append(float(part))
        except ValueError:
            continue
    return counts


def _carrier_gt(gt: str) -> bool:
    alleles = [part for part in gt.replace("|", "/").split("/") if part not in {"", "."}]
    return any(part not in {"0"} for part in alleles)


def _expanded(locus: str, counts: list[float], cfg: dict) -> tuple[object, str]:
    meta = (cfg.get("str_loci") or {}).get(locus) or {}
    cutoff = meta.get("pathogenic_min_repeats")
    citation = meta.get("citation")
    if cutoff is None or not citation or not counts:
        return pd.NA, ""
    return bool(max(counts) >= float(cutoff)), str(citation)


def _gene_for_locus(locus: str, cfg: dict) -> str:
    meta = (cfg.get("str_loci") or {}).get(locus) or {}
    return str(meta.get("gene") or locus or "STR")


def _assessment_rows(
    module: str,
    universe: list[str],
    carriers: set[str],
    seen: set[str],
    configured: bool,
) -> list[dict]:
    meta = MODULES[module]
    rows = []
    for sample in universe:
        if not configured:
            status, detail = "not_assessed", meta["not_assessed"]
        elif sample not in seen:
            status, detail = "not_assessed", meta["not_assessed"]
        elif sample in carriers:
            status, detail = "carrier", meta["carrier"]
        else:
            status, detail = "no_variant", meta["absent"]
        rows.append({"sample_id": sample, "module": module, "callset_status": detail, "state": status})
    return rows


def _prefix(samples: dict[str, dict], prefix: str | None) -> dict[str, dict]:
    if not prefix:
        return samples
    return {sample: fmt for sample, fmt in samples.items() if str(sample).startswith(prefix)}


def parse_expansionhunter(path: str | Path, cfg: dict | None = None, sample_prefix: str | None = None) -> tuple[pd.DataFrame, set[str], set[str]]:
    """Return variant rows, header samples, and carrier samples."""
    cfg = cfg or {}
    path = Path(path)
    header = [sample for sample in read_samples(path) if not sample_prefix or sample.startswith(sample_prefix)]
    rows: list[dict] = []
    carriers: set[str] = set()
    for record in iter_vcf(path):
        info = record["info"]
        locus = str(info.get("REPID") or info.get("VARID") or record["id"] or "STR")
        if locus in {".", ""}:
            locus = "STR"
        gene = _gene_for_locus(locus, cfg)
        ref_repeats = _float_or_none(info.get("REF"))
        for sample, fmt in _prefix(record["samples"], sample_prefix).items():
            gt = str(fmt.get("GT") or "./.")
            counts = _rep_counts(fmt.get("REPCN") or info.get("REPCN"))
            if not _carrier_gt(gt):
                continue
            carriers.add(str(sample))
            expanded, citation = _expanded(locus, counts, cfg)
            repeat_count = max(counts) if counts else pd.NA
            alt = record["alt"].split(",")[0]
            chrom = chrom_bare(record["chrom"])
            rows.append(
                _blank_row(
                    sample_id=str(sample),
                    family_id=str(sample),
                    chrom=chrom,
                    pos=record["pos"],
                    end=_float_or_none(info.get("END")) or record["pos"],
                    ref=record["ref"],
                    alt=alt,
                    variant_id=f"{chrom}-{record['pos']}-{record['ref']}-{alt}-{locus}",
                    gene=gene,
                    consequence="short_tandem_repeat_variant",
                    genotype=gt,
                    gq=_float_or_none(fmt.get("GQ")),
                    depth=_float_or_none(fmt.get("DP")),
                    clinvar="",
                    clinvar_class="none",
                    variant_class="STR",
                    hit_source="expansionhunter",
                    repeat_count=repeat_count if repeat_count is not pd.NA else (ref_repeats if ref_repeats is not None else pd.NA),
                    str_locus=locus,
                    str_expanded=expanded,
                    str_citation=citation,
                    scored_dosage=False,
                )
            )
    table = pd.DataFrame(rows) if rows else empty_variant_table()
    if not table.empty:
        table = derive_zygosity(table)
    return table, set(header), carriers


def _column(df: pd.DataFrame, names: list[str]) -> pd.Series | None:
    lookup = {str(column).lower(): column for column in df.columns}
    for name in names:
        if name.lower() in lookup:
            return df[lookup[name.lower()]]
    return None


def _as_bool(value: object) -> bool:
    text = str(value).strip().lower()
    return text in {"1", "true", "yes", "y", "pathogenic", "p", "p/lp", "likely_pathogenic", "likely pathogenic"}


def read_gauchian_table(path: str | Path) -> pd.DataFrame:
    path = Path(path)
    if path.suffix.lower() == ".json":
        with path.open(encoding="utf-8") as handle:
            payload = json.load(handle)
        if isinstance(payload, dict):
            if isinstance(payload.get("samples"), list):
                payload = payload["samples"]
            else:
                payload = [payload]
        if not isinstance(payload, list):
            raise ValueError(f"{path} is not a Gauchian JSON object or list")
        return pd.DataFrame(payload)
    from .io import read_table

    return read_table(path)


def parse_gauchian(path: str | Path, sample_prefix: str | None = None, default_sample: str | None = None) -> tuple[pd.DataFrame, set[str], set[str]]:
    """Map explicit GBA1 deletions and CN 0 to transcript ablation.

    A recombinant allele stays ``structural_variant`` unless the file itself
    marks it pathogenic. CN below the diploid baseline is treated as a
    deletion, which is the same mechanistic rule used for SV copy loss.
    """
    raw = read_gauchian_table(path)
    if raw.empty:
        return empty_variant_table(), set(), set()
    samples = _column(raw, ["sample_id", "sample", "Sample"])
    if samples is None:
        if default_sample is None:
            raise ValueError(f"{path} has no sample_id column")
        samples = pd.Series([default_sample] * len(raw))
    copy_number = _column(raw, ["gba_copy_number", "copy_number", "cn", "gba_cn"])
    recombinant = _column(raw, ["recombinant_allele", "recombinant", "call", "variant"])
    deletion = _column(raw, ["deletion", "deletion_breakpoint", "is_deletion"])
    pathogenic = _column(raw, ["pathogenic", "classification", "clinical_significance"])
    rows = []
    seen: set[str] = set()
    carriers: set[str] = set()
    for index in range(len(raw)):
        sample = str(samples.iloc[index])
        if sample_prefix and not sample.startswith(sample_prefix):
            continue
        seen.add(sample)
        cn = _float_or_none(copy_number.iloc[index]) if copy_number is not None else None
        rec = "" if recombinant is None else str(recombinant.iloc[index] or "")
        dele = "" if deletion is None else str(deletion.iloc[index] or "")
        path_flag = False if pathogenic is None else _as_bool(pathogenic.iloc[index])
        rec_text = rec.strip()
        if rec_text.lower() in {"", "nan", "none", ".", "false", "0"}:
            rec_text = ""
        deletion_text = dele.strip()
        explicit_deletion = bool(deletion_text) and deletion_text.lower() not in {"nan", "none", ".", "false", "0"}
        copy_loss = cn is not None and cn < 2
        if "del" in rec_text.lower():
            explicit_deletion = True
        if not explicit_deletion and not copy_loss and not rec_text and not path_flag:
            continue
        carriers.add(sample)
        if explicit_deletion or (cn is not None and cn == 0) or (copy_loss and not rec_text):
            consequence = "transcript_ablation"
            mechanism = "lof"
            event = "deletion"
            scored = True
            variant_class = "SV"
            alt = "<DEL>"
        else:
            consequence = "structural_variant"
            mechanism = "recombinant"
            event = "recombinant"
            scored = False
            variant_class = "SV"
            alt = rec_text or "<RECOMBINANT>"
        zyg = "hom" if cn == 0 else "het"
        genotype = "1/1" if zyg == "hom" else "0/1"
        clinvar_class = "P/LP" if path_flag else "none"
        rows.append(
            _blank_row(
                sample_id=sample,
                family_id=sample,
                chrom="1",
                pos=155235000,
                end=155235000,
                ref="N",
                alt=alt,
                variant_id=f"1-155235000-N-{alt}-{sample}",
                gene="GBA1",
                consequence=consequence,
                genotype=genotype,
                zygosity=zyg,
                clinvar="" if not path_flag else "Pathogenic",
                clinvar_class=clinvar_class,
                variant_class=variant_class,
                hit_source="gauchian",
                copy_number=cn if cn is not None else pd.NA,
                dosage_mechanism=mechanism,
                dosage_event=event,
                scored_dosage=scored,
            )
        )
    table = pd.DataFrame(rows) if rows else empty_variant_table()
    if not table.empty:
        table = derive_zygosity(table)
    return table, seen, carriers


def parse_mtdna(path: str | Path, sample_prefix: str | None = None) -> tuple[pd.DataFrame, set[str], set[str]]:
    path = Path(path)
    header = [sample for sample in read_samples(path) if not sample_prefix or sample.startswith(sample_prefix)]
    rows: list[dict] = []
    carriers: set[str] = set()
    for record in iter_vcf(path):
        info = record["info"]
        gene = str(info.get("GENE") or info.get("SYMBOL") or info.get("gene") or "MT")
        for sample, fmt in _prefix(record["samples"], sample_prefix).items():
            gt = str(fmt.get("GT") or "./.")
            heteroplasmy = _float_or_none(fmt.get("AF") or info.get("AF"))
            if not _carrier_gt(gt) and not (heteroplasmy is not None and heteroplasmy > 0):
                continue
            carriers.add(str(sample))
            alt = record["alt"].split(",")[0]
            chrom = chrom_bare(record["chrom"])
            rows.append(
                _blank_row(
                    sample_id=str(sample),
                    family_id=str(sample),
                    chrom=chrom,
                    pos=record["pos"],
                    end=record["pos"],
                    ref=record["ref"],
                    alt=alt,
                    variant_id=f"{chrom}-{record['pos']}-{record['ref']}-{alt}",
                    gene=gene,
                    consequence="mitochondrial_variant",
                    genotype=gt if _carrier_gt(gt) else "0/1",
                    gq=_float_or_none(fmt.get("GQ")),
                    depth=_float_or_none(fmt.get("DP")),
                    clinvar="",
                    clinvar_class="none",
                    variant_class="MT",
                    hit_source="mtdna",
                    heteroplasmy=heteroplasmy if heteroplasmy is not None else pd.NA,
                    scored_dosage=False,
                )
            )
    table = pd.DataFrame(rows) if rows else empty_variant_table()
    if not table.empty:
        table = derive_zygosity(table)
    return table, set(header), carriers


def _concat(frames: list[pd.DataFrame]) -> pd.DataFrame:
    frames = [frame for frame in frames if frame is not None and not frame.empty]
    if not frames:
        return empty_variant_table()
    return pd.concat(frames, ignore_index=True)



def _restrict(table: pd.DataFrame, sample_id: str, seen: set[str]) -> tuple[pd.DataFrame, bool, bool]:
    """Keep rows for one sample. A single foreign sample id is relabelled."""
    sample_id = str(sample_id)
    seen = {str(item) for item in seen}
    if table is None or table.empty:
        return empty_variant_table() if table is None or table.empty else table, sample_id in seen, False
    ids = set(table["sample_id"].astype(str))
    if sample_id in ids:
        kept = table[table["sample_id"].astype(str) == sample_id].copy()
        return kept, True, not kept.empty
    # A one-sample file whose header name differs from sample_id is this sample.
    # A multi-sample file that merely has one carrier row is not relabelled.
    if len(seen) == 1:
        kept = table.copy()
        kept["sample_id"] = sample_id
        if "family_id" in kept.columns:
            kept["family_id"] = sample_id
        return kept, True, not kept.empty
    return table.iloc[0:0].copy(), sample_id in seen, False


def ingest_sample_callsets(
    sample_id: str,
    *,
    expansionhunter: Path | None = None,
    gauchian: Path | None = None,
    mtdna: Path | None = None,
    cfg: dict | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Ingest the three per-sample callsets. A missing file stays not assessed."""
    cfg = cfg or {}
    sample_id = str(sample_id)
    frames: list[pd.DataFrame] = []
    assessments: list[dict] = []

    def add(module: str, path: Path | None, parser) -> None:
        if path is None or not Path(path).exists():
            assessments.extend(_assessment_rows(module, [sample_id], set(), set(), False))
            return
        parsed = parser(Path(path))
        table, seen, carriers = parsed
        table, present, carrier = _restrict(table, sample_id, seen)
        frames.append(table)
        assessments.extend(
            _assessment_rows(
                module,
                [sample_id],
                {sample_id} if carrier else set(),
                {sample_id} if present else set(),
                True,
            )
        )

    add("expansionhunter", expansionhunter, lambda path: parse_expansionhunter(path, cfg))
    add("gauchian", gauchian, lambda path: parse_gauchian(path, default_sample=sample_id))
    add("mtdna", mtdna, lambda path: parse_mtdna(path))
    return _concat(frames), pd.DataFrame(assessments)
