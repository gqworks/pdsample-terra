"""GATK-SV ingest for the panel, with dosage rules for SNCA and PRKN.

Records are streamed. A missing index is fine for the ~41 MB SV VCF; the
20 GB SNV VCF is never opened here. Samples present in the SNV callset but
absent from the SV header are reported as ``SV not assessed``.
"""

from __future__ import annotations

import re
from pathlib import Path

import pandas as pd

from .normalize import apply_gene_aliases, derive_zygosity
from .scoring import assign_tier
from .vcf_stream import IntervalIndex, chrom_bare, ensure_tabix, iter_vcf, parse_info

PREDICTED_FIELDS = {
    "PREDICTED_LOF": "LOF",
    "PREDICTED_COPY_GAIN": "COPY_GAIN",
    "PREDICTED_DUP_PARTIAL": "DUP_PARTIAL",
    "PREDICTED_INTRAGENIC_EXON_DUP": "INTRAGENIC_EXON_DUP",
    "PREDICTED_INTRONIC": "INTRONIC",
    "PREDICTED_INV_SPAN": "INV_SPAN",
    "PREDICTED_MSV_EXON_OVERLAP": "MSV_EXON_OVERLAP",
    "PREDICTED_PARTIAL_EXON_DUP": "PARTIAL_EXON_DUP",
    "PREDICTED_PROMOTER": "PROMOTER",
    "PREDICTED_TSS_DUP": "TSS_DUP",
    "PREDICTED_UTR": "UTR",
    "PREDICTED_BREAKEND_EXONIC": "BREAKEND_EXONIC",
    "PREDICTED_PARTIAL_DISPERSED_DUP": "PARTIAL_DISPERSED_DUP",
}

# Mechanism labels, not penetrance. Citations are for the dosage biology.
DEFAULT_DOSAGE_RULES = [
    {
        "gene": "SNCA",
        "effect": "copy_gain",
        "citation": "Singleton AB et al. Science. 2003;302:841. Chartier-Harlin MC et al. Lancet. 2004;364:1167.",
    },
    {
        "gene": "PRKN",
        "effect": "lof",
        "citation": "Kitada T et al. Nature. 1998;392:605. GeneReviews Parkinson Disease Overview (NBK1223).",
    },
]

DEFAULT_EXCLUDE_FILTERS = ["HIGH_ALGORITHM_FDR", "HIGH_NCR", "UNRESOLVED"]

LOF_PREDICTED = {"LOF", "INTRAGENIC_EXON_DUP", "BREAKEND_EXONIC", "PARTIAL_EXON_DUP"}
GAIN_PREDICTED = {"COPY_GAIN", "TSS_DUP"}

VARIANT_COLUMNS = [
    "sample_id", "family_id", "chrom", "pos", "end", "ref", "alt", "variant_id",
    "gene", "gene_id", "consequence", "genotype", "zygosity", "gq", "allele_balance",
    "depth", "cohort_af", "gnomad_af", "gnomad_popmax_af", "cadd", "revel", "spliceai",
    "clinvar", "clinvar_class", "variant_class", "hit_source", "svtype", "sv_id",
    "copy_number", "expected_cn", "dosage_mechanism", "dosage_event", "dosage_citation",
    "scored_dosage", "predicted_effects", "exons_overlapped", "exon_count", "filter",
    "cn_inferred_gt",
]


def _as_float(value: object) -> float | None:
    if value is None or value is True or value == "." or value == "":
        return None
    try:
        return float(str(value).split(",")[0])
    except ValueError:
        return None


def _as_int(value: object) -> int | None:
    number = _as_float(value)
    if number is None:
        return None
    return int(number)


def _gene_tokens(value: object) -> list[str]:
    if value is None or value is True:
        return []
    return [part.strip().upper() for part in re.split(r"[,|]", str(value)) if part.strip() and part.strip() != "."]


def _allele_index(gt: str) -> int:
    alleles = [a for a in re.split(r"[/|]", gt) if a not in {"", "."}]
    alts = [int(a) for a in alleles if a.isdigit() and int(a) > 0]
    return (min(alts) - 1) if alts else 0


def _number_for_allele(value: object, index: int) -> float | None:
    if value is None or value is True:
        return None
    parts = str(value).split(",")
    if index < len(parts):
        return _as_float(parts[index])
    return _as_float(parts[0]) if parts else None


def passes_filter(filt: str, exclude: set[str]) -> bool:
    if filt in {"", ".", "PASS", None}:
        return True
    tags = set(str(filt).split(";"))
    return not bool(tags & exclude)


def svtype_of(info: dict, alt: str) -> str:
    raw = info.get("SVTYPE")
    if isinstance(raw, str) and raw and raw != ".":
        return raw.split(":")[0]
    match = re.search(r"<([A-Z]+)", alt or "")
    return match.group(1) if match else "SV"


def gnomad_frequencies(info: dict) -> tuple[float | None, float | None]:
    af = _as_float(info.get("gnomad_v4.1_sv_AF"))
    pops = []
    for key, value in info.items():
        if key.startswith("gnomad_v4.1_sv_") and key.endswith("_AF") and key != "gnomad_v4.1_sv_AF":
            number = _as_float(value)
            if number is not None:
                pops.append(number)
    popmax = max(pops) if pops else af
    return af, popmax


def cohort_frequency(info: dict, svtype: str, allele_index: int) -> float | None:
    if svtype == "CNV" or "CN_NONREF_FREQ" in info:
        return _as_float(info.get("CN_NONREF_FREQ"))
    return _number_for_allele(info.get("AF"), allele_index)


def dosage_event(svtype: str, cn: int | None, ecn: int | None, predicted: set[str]) -> str:
    loss = cn is not None and ecn is not None and cn < ecn
    gain = cn is not None and ecn is not None and cn > ecn
    if svtype == "DEL" or loss or "LOF" in predicted:
        event = "deletion"
    elif svtype in {"DUP", "CNV"} or gain or predicted & GAIN_PREDICTED:
        event = "duplication"
    else:
        event = svtype.lower()
    if cn is not None and ecn is not None and cn >= ecn + 2 and event == "duplication":
        return "triplication_or_higher"
    return event


def _rule_matches(
    rule: dict,
    svtype: str,
    event: str,
    predicted: set[str],
    n_exons: int,
    n_gene_exons: int,
    cn: int | None,
    ecn: int | None,
) -> bool:
    if rule.get("effect") == "copy_gain":
        return event in {"duplication", "triplication_or_higher"} or bool(predicted & GAIN_PREDICTED)
    if rule.get("effect") != "lof":
        return False
    if predicted & {"LOF", "INTRAGENIC_EXON_DUP"}:
        return True
    if svtype == "DEL" and n_exons > 0:
        return True
    partial = n_gene_exons > 0 and 0 < n_exons < n_gene_exons
    if svtype in {"DUP", "CNV"} and partial:
        return True
    return bool(cn is not None and ecn is not None and cn < ecn and n_exons > 0)


def classify_dosage(
    gene: str,
    svtype: str,
    cn: int | None,
    ecn: int | None,
    predicted: set[str],
    n_exons: int,
    n_gene_exons: int,
    rules: list[dict] | None = None,
) -> tuple[str, str, bool, str]:
    """Return mechanism, event, whether to score it, and the citation.

    ``scored`` is true for configured dosage genes and for exon-disrupting
    deletions. Copy gains are scored only for genes whose rule says ``copy_gain``
    (SNCA duplication/triplication). The score itself reuses ``null_variant`` points.
    """
    rules = rules if rules is not None else DEFAULT_DOSAGE_RULES
    event = dosage_event(svtype, cn, ecn, predicted)
    for rule in rules:
        if str(rule.get("gene", "")).upper() != gene:
            continue
        if _rule_matches(rule, svtype, event, predicted, n_exons, n_gene_exons, cn, ecn):
            return str(rule["effect"]), event, True, str(rule.get("citation") or "")
    if predicted & LOF_PREDICTED or (svtype == "DEL" and n_exons > 0):
        return "lof", event if event == "deletion" else "deletion", True, "Exon disruption (GATK-SV PREDICTED_* or exon BED overlap)."
    if predicted & GAIN_PREDICTED:
        return "copy_gain", event, False, ""
    return "other", event, False, ""


def _consequence(mechanism: str, scored: bool, n_exons: int) -> str:
    if mechanism == "lof" and scored:
        return "transcript_ablation"
    if mechanism == "copy_gain":
        return "copy_number_gain"
    if n_exons > 0:
        return "structural_variant"
    return "intron_variant"


def _zygosity(gt: str, cn: int | None, ecn: int | None) -> tuple[str | None, str, bool]:
    """Return zygosity, genotype string, and whether the genotype was inferred from CN."""
    alleles = re.split(r"[/|]", gt) if gt else []
    called = [a for a in alleles if a not in {"", "."}]
    alts = [a for a in called if a != "0"]
    if called and alts:
        if len(alts) == len(called) and len(set(alts)) == 1 and len(called) > 1:
            return "hom", gt, False
        if len(called) == 1:
            return "hemi", gt, False
        return "het", gt, False
    if called and not alts:
        return "ref", gt, False
    if cn is None:
        return None, gt or "./.", False
    baseline = 2 if ecn is None else ecn
    if cn == baseline:
        return "ref", gt or "./.", False
    if cn == 0:
        return "hom", "1/1", True
    return "het", "0/1", True


def _is_carrier(gt: str, cn: int | None, ecn: int | None) -> bool:
    zyg, _, _ = _zygosity(gt, cn, ecn)
    return zyg in {"het", "hom", "hemi"}


def _record_end(info: dict, pos: int) -> int:
    end = _as_int(info.get("END"))
    return end if end is not None else pos


def _exon_hits(index: IntervalIndex, info: dict, chrom: str, pos: int, end: int) -> list[dict]:
    hits = list(index.overlap(chrom, pos - 1, end))
    raw = info.get("CPX_INTERVALS")
    if isinstance(raw, str):
        for match in re.finditer(r"(chr[\w]+):(\d+)-(\d+)", raw):
            hits.extend(index.overlap(match.group(1), int(match.group(2)) - 1, int(match.group(3))))
    return hits


def _effects_by_gene(info: dict) -> dict[str, set[str]]:
    found: dict[str, set[str]] = {}
    for field, label in PREDICTED_FIELDS.items():
        for gene in _gene_tokens(info.get(field)):
            found.setdefault(gene, set()).add(label)
    return found


def empty_variant_table() -> pd.DataFrame:
    return pd.DataFrame(columns=VARIANT_COLUMNS)


def coverage_table(
    snv_samples: list[str] | None,
    sv_samples: list[str],
    carrier_samples: list[str],
) -> pd.DataFrame:
    """``SV not assessed`` when the sample is absent from the SV VCF header."""
    header = {str(sample) for sample in sv_samples}
    carriers = {str(sample) for sample in carrier_samples}
    samples = [str(sample) for sample in (snv_samples if snv_samples else sv_samples)]
    rows = []
    for sample in samples:
        if sample not in header:
            status = "SV not assessed"
        elif sample in carriers:
            status = "SV carrier"
        else:
            status = "no SV"
        rows.append({"sample_id": sample, "sv_status": status})
    return pd.DataFrame(rows)


def ingest_sv(
    vcf: str | Path,
    gene_bed: str | Path,
    exon_bed: str | Path,
    panel: pd.DataFrame,
    *,
    aliases: dict[str, str] | None = None,
    sample_prefix: str | None = None,
    exclude_filters: list[str] | None = None,
    dosage_rules: list[dict] | None = None,
    snv_samples: list[str] | None = None,
    auto_index: bool = True,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Return ``(variant_table, sv_assessment)``.

    The variant table uses the same columns the SNV scorer expects, one row
    per carrier sample and panel gene.
    """
    vcf = Path(vcf)
    if auto_index:
        ensure_tabix(vcf)
    from .vcf_stream import read_samples

    header_samples = read_samples(vcf)
    exclude = set(exclude_filters if exclude_filters is not None else DEFAULT_EXCLUDE_FILTERS)
    rules = dosage_rules if dosage_rules is not None else DEFAULT_DOSAGE_RULES
    panel_genes = set(apply_gene_aliases(panel["gene"], aliases))
    gene_ids = {}
    if "gene_id" in panel.columns:
        gene_ids = dict(zip(apply_gene_aliases(panel["gene"], aliases), panel["gene_id"].astype(str)))
    exons = IntervalIndex.from_bed(exon_bed)
    genes = IntervalIndex.from_bed(gene_bed)
    exon_counts: dict[str, int] = {}
    for spans in exons.by_chrom.values():
        for _, _, payload in spans:
            if payload["gene"]:
                exon_counts[payload["gene"]] = exon_counts.get(payload["gene"], 0) + 1

    rows: list[dict] = []
    for record in iter_vcf(vcf, regions=gene_bed):
        if not passes_filter(record["filter"], exclude):
            continue
        info = record["info"]
        end = _record_end(info, record["pos"])
        svtype = svtype_of(info, record["alt"])
        raw_effects = _effects_by_gene(info)
        effects: dict[str, set[str]] = {}
        for symbol, labels in raw_effects.items():
            aliased = str(apply_gene_aliases(pd.Series([symbol]), aliases).iloc[0])
            effects.setdefault(aliased, set()).update(labels)
        exon_hit_rows = _exon_hits(exons, info, record["chrom"], record["pos"], end)
        exons_by_gene: dict[str, list[dict]] = {}
        for hit in exon_hit_rows:
            symbol = str(apply_gene_aliases(pd.Series([hit["gene"]]), aliases).iloc[0])
            exons_by_gene.setdefault(symbol, []).append(hit)
        gene_symbols = set(effects) | set(exons_by_gene)
        for hit in genes.overlap(record["chrom"], record["pos"] - 1, end):
            if hit["gene"]:
                gene_symbols.add(str(apply_gene_aliases(pd.Series([hit["gene"]]), aliases).iloc[0]))
        gene_symbols = {g for g in gene_symbols if g in panel_genes}
        if not gene_symbols:
            continue
        gnomad_af, popmax = gnomad_frequencies(info)
        alt = record["alt"].split(",")[0]
        sv_id = record["id"] if record["id"] not in {"", "."} else f"{chrom_bare(record['chrom'])}-{record['pos']}-{end}-{svtype}"
        for sample, fmt in record["samples"].items():
            if sample_prefix and not str(sample).startswith(sample_prefix):
                continue
            gt = fmt.get("GT", "./.")
            cn = _as_int(fmt.get("CN"))
            ecn = _as_int(fmt.get("ECN"))
            if not _is_carrier(gt, cn, ecn if ecn is not None else 2):
                continue
            zyg, genotype, inferred = _zygosity(gt, cn, ecn if ecn is not None else 2)
            gq = _as_float(fmt.get("GQ"))
            if gq is None:
                gq = _as_float(fmt.get("CNQ"))
            if gq is None:
                gq = _as_float(fmt.get("RD_GQ"))
            allele_index = _allele_index(gt)
            cohort_af = cohort_frequency(info, svtype, allele_index)
            for gene in sorted(gene_symbols):
                predicted = effects.get(gene, set())
                n_exons = len({(h["chrom"], h["start"], h["end"]) for h in exons_by_gene.get(gene, [])})
                mechanism, event, scored, citation = classify_dosage(
                    gene, svtype, cn, ecn if ecn is not None else 2, predicted,
                    n_exons, exon_counts.get(gene, 0), rules,
                )
                rows.append(
                    {
                        "sample_id": str(sample),
                        "family_id": str(sample),
                        "chrom": chrom_bare(record["chrom"]),
                        "pos": record["pos"],
                        "end": end,
                        "ref": record["ref"],
                        "alt": alt,
                        "variant_id": sv_id,
                        "gene": gene,
                        "gene_id": gene_ids.get(gene, ""),
                        "consequence": _consequence(mechanism, scored, n_exons),
                        "genotype": genotype,
                        "zygosity": zyg,
                        "gq": gq,
                        "allele_balance": pd.NA,
                        "depth": pd.NA,
                        "sv_panel_af": cohort_af,
                        "cohort_af": pd.NA,
                        "gnomad_af": gnomad_af,
                        "gnomad_popmax_af": popmax,
                        "cadd": pd.NA,
                        "revel": pd.NA,
                        "spliceai": pd.NA,
                        "clinvar": "",
                        "clinvar_class": "none",
                        "variant_class": "SV",
                        "hit_source": "sv",
                        "svtype": svtype,
                        "sv_id": sv_id,
                        "copy_number": cn if cn is not None else pd.NA,
                        "expected_cn": ecn if ecn is not None else 2,
                        "dosage_mechanism": mechanism,
                        "dosage_event": event,
                        "dosage_citation": citation,
                        "scored_dosage": scored,
                        "predicted_effects": ";".join(sorted(predicted)),
                        "exons_overlapped": n_exons,
                        "exon_count": exon_counts.get(gene, 0),
                        "filter": record["filter"],
                        "cn_inferred_gt": inferred,
                    }
                )
    table = pd.DataFrame(rows) if rows else empty_variant_table()
    if not table.empty:
        table = derive_zygosity(table)
        # derive_zygosity keeps a provided zygosity. CN-inferred genotypes are already set.
    carriers = table["sample_id"].astype(str).unique().tolist() if not table.empty else []
    assessment = coverage_table(snv_samples, header_samples, carriers)
    return table, assessment


def apply_dosage_evidence(df: pd.DataFrame, scfg: dict) -> pd.DataFrame:
    """Add the existing ``null_variant`` points once for scored dosage SVs.

    LoF SVs already receive those points via ``transcript_ablation``. Copy-gain
    events on configured genes (SNCA) receive the same configured point value,
    recorded as ``sv_dosage_gain``. No new point value is introduced.
    """
    if df.empty or "dosage_mechanism" not in df.columns or "evidence_score" not in df.columns:
        return df
    df = df.copy()
    points = scfg["points"]["null_variant"]
    reasons = df["evidence_reasons"].fillna("").astype(str)
    if "scored_dosage" in df.columns:
        scored = df["scored_dosage"].fillna(False).astype(bool)
    else:
        scored = df["dosage_mechanism"].isin(["lof", "copy_gain"])
    gain = scored & (df["dosage_mechanism"] == "copy_gain") & ~reasons.str.contains("null_variant|sv_dosage")
    df.loc[gain, "evidence_score"] = df.loc[gain, "evidence_score"] + points
    df.loc[gain, "evidence_reasons"] = (reasons.loc[gain] + ";sv_dosage_gain").str.strip(";")
    reasons = df["evidence_reasons"].fillna("").astype(str)
    label = scored & (df["dosage_mechanism"] == "lof") & ~reasons.str.contains("sv_dosage")
    df.loc[label, "evidence_reasons"] = (reasons.loc[label] + ";sv_dosage_lof").str.strip(";")
    df["tier"] = assign_tier(df["evidence_score"], scfg["tiers"])
    return df


def parse_info_field(text: str) -> dict[str, str | bool]:
    """Public alias so tests can target the INFO grammar without a whole VCF."""
    return parse_info(text)
