"""Ingest one sample's ExpansionHunter, Gauchian, and mtDNA outputs.

No STR pathogenicity cutoff is applied unless config sets both
``pathogenic_min_repeats`` and a citation.
"""

from __future__ import annotations

import json
import re
from pathlib import Path

import pandas as pd

from .normalize import derive_zygosity
from .vcf_stream import chrom_bare, iter_vcf, read_samples

# Gauchian 1.0.2 ``other_variants`` entries look like ``c.681T>G(p.Asn227Lys)``.
_OTHER_VARIANT = re.compile(r"(c\.[^(,\s]+)(?:\((p\.[^)]+)\))?")
_GBA1_LOCUS_POS = 155235000

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
    "cohort_af", "gnomad_af", "gnomad_popmax_af", "cadd", "revel",
    "sift", "sift_pred", "polyphen", "polyphen_pred", "eigen",
    "muttaster", "muttaster_pred", "fathmm", "fathmm_pred",
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


def _empty_gauchian_text(value: object) -> bool:
    text = str(value or "").strip().lower()
    return text in {"", "nan", "none", ".", "false", "0"}


def _normalize_gauchian_payload(payload: object, path: Path) -> list[dict]:
    """Accept list rows, a flat sample object, or native ``{sample_id: metrics}`` JSON."""
    if isinstance(payload, list):
        rows = []
        for item in payload:
            if isinstance(item, dict):
                rows.append(item)
            else:
                raise ValueError(f"{path} list entries must be objects")
        return rows
    if not isinstance(payload, dict):
        raise ValueError(f"{path} is not a Gauchian JSON object or list")
    if isinstance(payload.get("samples"), list):
        return [item for item in payload["samples"] if isinstance(item, dict)]
    # Native Gauchian 1.0.2: {"SAMPLE": {"deletion_CN": 2, "other_variants": "...", ...}}
    flat_keys = {"sample_id", "sample", "gba_copy_number", "deletion_cn", "copy_number", "other_variants"}
    lowered = {str(key).lower() for key in payload}
    if payload and all(isinstance(value, dict) for value in payload.values()) and not (lowered & flat_keys):
        rows = []
        for sample_id, metrics in payload.items():
            if not isinstance(metrics, dict):
                continue
            row = dict(metrics)
            row.setdefault("sample_id", sample_id)
            rows.append(row)
        return rows
    return [payload]


def read_gauchian_table(path: str | Path) -> pd.DataFrame:
    path = Path(path)
    if path.suffix.lower() == ".json":
        with path.open(encoding="utf-8") as handle:
            payload = json.load(handle)
        rows = _normalize_gauchian_payload(payload, path)
        return pd.DataFrame(rows) if rows else pd.DataFrame()
    from .io import read_table

    return read_table(path)


def gba_c_dot(value: object) -> str:
    """Normalize HGVS c. notation for GBA1 dedup (``ENST...:c.681T>G`` / ``c.681T>G(p.*)``)."""
    text = str(value or "").strip()
    if "c." not in text:
        return ""
    body = text.split("c.", 1)[1]
    body = body.split("(", 1)[0].split(" ", 1)[0].strip()
    return f"c.{body}" if body else ""


def parse_other_variants(text: object) -> list[tuple[str, str]]:
    """Split Gauchian ``other_variants`` into ``(hgvsc, hgvsp)`` pairs."""
    if _empty_gauchian_text(text):
        return []
    out: list[tuple[str, str]] = []
    for chunk in str(text).split(","):
        chunk = chunk.strip()
        if not chunk:
            continue
        match = _OTHER_VARIANT.fullmatch(chunk) or _OTHER_VARIANT.search(chunk)
        if not match:
            continue
        hgvsc = match.group(1).strip()
        hgvsp = (match.group(2) or "").strip()
        if hgvsc:
            out.append((hgvsc, hgvsp))
    return out


def _consequence_from_protein(hgvsp: str) -> str:
    text = (hgvsp or "").lower()
    if "ter" in text or "*" in text:
        return "stop_gained"
    if "fs" in text:
        return "frameshift_variant"
    if text.startswith("p."):
        return "missense_variant"
    return "sequence_variant"


def _variant_raw_depth_ab(raw_counts: object, hgvsc: str, hgvsp: str) -> tuple[float | None, float | None]:
    """Gauchian ``variant_raw_count`` values are ``alt,ref`` supporting counts."""
    if not isinstance(raw_counts, dict):
        return None, None
    keys = [f"{hgvsc}({hgvsp})" if hgvsp else hgvsc, hgvsc]
    for key in keys:
        value = raw_counts.get(key)
        if value is None:
            continue
        parts = str(value).split(",")
        if len(parts) < 2:
            continue
        try:
            alt = float(parts[0])
            ref = float(parts[1])
        except ValueError:
            continue
        depth = alt + ref
        if depth <= 0:
            continue
        return depth, alt / depth
    return None, None


def drop_redundant_gauchian_snvs(df: pd.DataFrame) -> pd.DataFrame:
    """Prefer VEP/small-variant GBA1 alleles when Gauchian reports the same c. HGVS."""
    if df.empty or "hit_source" not in df.columns or "gene" not in df.columns:
        return df
    gene = df["gene"].astype(str).str.upper()
    hit = df["hit_source"].fillna("").astype(str)
    event = df["dosage_event"].fillna("").astype(str) if "dosage_event" in df.columns else pd.Series("", index=df.index)
    hgvsc = df["hgvsc"] if "hgvsc" in df.columns else pd.Series("", index=df.index)
    small = gene.eq("GBA1") & hit.eq("small_variant")
    gauchian = gene.eq("GBA1") & hit.eq("gauchian") & event.eq("other_variant")
    if not small.any() or not gauchian.any():
        return df
    seen = {gba_c_dot(value) for value in hgvsc[small] if gba_c_dot(value)}
    if not seen:
        return df
    drop = gauchian & hgvsc.map(lambda value: gba_c_dot(value) in seen)
    return df.loc[~drop].reset_index(drop=True)


def parse_gauchian(path: str | Path, sample_prefix: str | None = None, default_sample: str | None = None) -> tuple[pd.DataFrame, set[str], set[str]]:
    """Map Gauchian CN loss, recombinants, and catalog ``other_variants`` into GBA1 rows.

    Copy number 0 / explicit deletion becomes homozygous ``transcript_ablation``.
    A recombinant allele stays ``structural_variant`` unless the file itself marks
    it pathogenic. Pathogenic SNVs in ``other_variants`` (Gauchian 1.0.2) become
    missense/null rows with ``hit_source=gauchian`` so they are not missed when
    HaplotypeCaller fails in the GBAP1 homology region.
    """
    raw = read_gauchian_table(path)
    if raw.empty:
        return empty_variant_table(), set(), set()
    samples = _column(raw, ["sample_id", "sample", "Sample"])
    if samples is None:
        if default_sample is None:
            raise ValueError(f"{path} has no sample_id column")
        samples = pd.Series([default_sample] * len(raw))
    copy_number = _column(raw, ["gba_copy_number", "copy_number", "cn", "gba_cn", "deletion_CN", "deletion_cn"])
    recombinant = _column(
        raw,
        [
            "recombinant_allele",
            "recombinant",
            "call",
            "variant",
            "GBAP1_like_variant_exon9_11",
            "gbap1_like_variant_exon9_11",
        ],
    )
    deletion = _column(raw, ["deletion", "deletion_breakpoint", "is_deletion"])
    pathogenic = _column(raw, ["pathogenic", "classification", "clinical_significance"])
    other_variants = _column(raw, ["other_variants", "Other_variants"])
    raw_counts_col = _column(raw, ["variant_raw_count", "variant_raw_counts"])
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
        if _empty_gauchian_text(rec_text):
            rec_text = ""
        deletion_text = dele.strip()
        explicit_deletion = bool(deletion_text) and not _empty_gauchian_text(deletion_text)
        copy_loss = cn is not None and cn < 2
        if "del" in rec_text.lower():
            explicit_deletion = True
        other_text = "" if other_variants is None else other_variants.iloc[index]
        parsed_other = parse_other_variants(other_text)
        raw_counts = raw_counts_col.iloc[index] if raw_counts_col is not None else None
        structural = explicit_deletion or copy_loss or bool(rec_text) or path_flag
        if not structural and not parsed_other:
            continue
        carriers.add(sample)
        if structural:
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
                    pos=_GBA1_LOCUS_POS,
                    end=_GBA1_LOCUS_POS,
                    ref="N",
                    alt=alt,
                    variant_id=f"1-{_GBA1_LOCUS_POS}-N-{alt}-{sample}",
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
        for offset, (hgvsc, hgvsp) in enumerate(parsed_other):
            depth, allele_balance = _variant_raw_depth_ab(raw_counts, hgvsc, hgvsp)
            consequence = _consequence_from_protein(hgvsp)
            # Gauchian catalogs pathogenic GBA1 alleles; mark P/LP so scoring/ranking keep them.
            rows.append(
                _blank_row(
                    sample_id=sample,
                    family_id=sample,
                    chrom="1",
                    pos=_GBA1_LOCUS_POS + 1 + offset,
                    end=_GBA1_LOCUS_POS + 1 + offset,
                    ref="N",
                    alt=hgvsc,
                    variant_id=f"GBA1-{hgvsc}-{sample}",
                    gene="GBA1",
                    consequence=consequence,
                    hgvsc=hgvsc,
                    hgvsp=hgvsp,
                    genotype="0/1",
                    zygosity="het",
                    depth=depth if depth is not None else pd.NA,
                    allele_balance=allele_balance if allele_balance is not None else pd.NA,
                    clinvar="Pathogenic",
                    clinvar_class="P/LP",
                    variant_class="SNV",
                    hit_source="gauchian",
                    copy_number=cn if cn is not None else pd.NA,
                    dosage_mechanism="other_variant",
                    dosage_event="other_variant",
                    scored_dosage=False,
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
