"""GATK-SV ingest for the panel, with dosage rules for SNCA and PRKN.

Records are streamed. A missing index is fine for the ~41 MB SV VCF; the
20 GB SNV VCF is never opened here. Samples present in the SNV callset but
absent from the SV header are reported as ``SV not assessed``.
"""

from __future__ import annotations

import re
from pathlib import Path

import pandas as pd

from .normalize import apply_gene_aliases, derive_zygosity, normalize_clinvar
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
        "min_span_bp": 10_000,
        "require_copy_gain_or_exon": True,
        "citation": "Singleton AB et al. Science. 2003;302:841. Chartier-Harlin MC et al. Lancet. 2004;364:1167.",
    },
    {
        "gene": "PRKN",
        "effect": "lof",
        "citation": "Kitada T et al. Nature. 1998;392:605. GeneReviews Parkinson Disease Overview (NBK1223).",
    },
]

# DEPTH_LT_5KB / SAMPLE_WGD_OUTLIER are GATK-SV depth-call QC tags; they mark
# choppy low-coverage dips that look like exon DELs in IGV but are not reliable.
DEFAULT_EXCLUDE_FILTERS = [
    "HIGH_ALGORITHM_FDR",
    "HIGH_NCR",
    "UNRESOLVED",
    "DEPTH_LT_5KB",
    "SAMPLE_WGD_OUTLIER",
]

LOF_PREDICTED = {"LOF", "INTRAGENIC_EXON_DUP", "BREAKEND_EXONIC", "PARTIAL_EXON_DUP"}
GAIN_PREDICTED = {"COPY_GAIN", "TSS_DUP"}

VARIANT_COLUMNS = [
    "sample_id", "family_id", "chrom", "pos", "end", "ref", "alt", "variant_id",
    "gene", "gene_id", "consequence", "genotype", "zygosity", "gq", "allele_balance",
    "depth", "cohort_af", "gnomad_af", "gnomad_popmax_af", "cadd", "revel",
    "sift", "sift_pred", "polyphen", "polyphen_pred", "eigen",
    "muttaster", "muttaster_pred", "fathmm", "fathmm_pred",
    "clinvar", "clinvar_class", "variant_class", "hit_source", "svtype", "svlen", "sv_id",
    "copy_number", "expected_cn", "dosage_mechanism", "dosage_event", "dosage_citation",
    "scored_dosage", "predicted_effects", "exons_overlapped", "exon_count", "filter",
    "cn_inferred_gt", "sv_algorithms", "sv_evidence", "rd_gq", "ogq",
    "high_sr_background", "freq_het", "sv_artifact",
]

# Default QC for depth-only IGV artifacts that remain FILTER=PASS.
DEFAULT_DEPTH_ONLY_ARTIFACT = {
    "enabled": True,
    "max_rd_gq": 5.0,
}

# Raw Manta (and similar PE/SR) priors lack FORMAT CN / RD_CN depth genotyping.
# Multi-Mb PASS DELs then exon-overlap into false transcript_ablation solves.
DEFAULT_PESR_CALLER_QC = {
    "enabled": True,
    # Absolute span cap for DEL/DUP when CN is absent (|SVLEN| or END-POS).
    "max_span_bp": 2_000_000,
    # DEL/DUP at or above this span require FORMAT CN or RD_CN (depth genotype). Raw Manta
    # junction-only calls in the 100 kb–2 Mb band (e.g. false 1.66 Mb TWNK DEL)
    # otherwise pass SR checks and score as transcript_ablation.
    "require_cn_min_span_bp": 100_000,
    # Drop INFO/IMPRECISE DEL/DUP when CN is absent.
    "drop_imprecise": True,
    # Above this span, require Manta-style SR alt support (FORMAT SR = ref,alt).
    "large_event_min_span_bp": 50_000,
    "large_event_require_sr": True,
    # When JUNCTION_QUAL is present, drop below this floor (covers JUNCTION_QUAL=1 junk).
    "min_junction_qual": 20.0,
    # Homozygous DEL/DUP without CN above this span are almost never real in adults.
    "max_hom_span_bp": 50_000,
}

_PESR_EVIDENCE = frozenset({"PE", "SR"})
_SPAN_FILTER_TYPES = frozenset({"DEL", "DUP"})


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


# Allowed FILTER values for SV ingest (GATK-SV keeps MULTIALLELIC as informative).
ALLOWED_SV_FILTERS = {"PASS", "MULTIALLELIC"}


def passes_filter(filt: str, exclude: set[str] | None = None) -> bool:
    """Keep ``FILTER=PASS`` or ``MULTIALLELIC`` SVs; drop all other tags.

    ``exclude`` still drops listed tags even if they co-occur with an allowed
    label (e.g. ``PASS;DEPTH_LT_5KB`` would fail via exclude if configured).
    """
    if filt is None:
        return False
    tags = {part for part in str(filt).split(";") if part and part not in {"", "."}}
    if not tags or not tags.issubset(ALLOWED_SV_FILTERS):
        return False
    if exclude and tags & set(exclude):
        return False
    return True


def _csv_tokens(value: object) -> set[str]:
    if value is None or value is True or value == "." or value == "":
        return set()
    return {part.strip().upper() for part in re.split(r"[,;]", str(value)) if part.strip() and part.strip() != "."}


def is_depth_only_algorithms(info: dict) -> bool:
    """True when INFO/ALGORITHMS is exclusively ``depth``."""
    algs = _csv_tokens(info.get("ALGORITHMS"))
    return bool(algs) and algs == {"DEPTH"}


def has_pesr_format_evidence(fmt: dict) -> bool:
    """True when FORMAT ``EV`` includes paired-end or split-read support."""
    return bool(_csv_tokens(fmt.get("EV")) & _PESR_EVIDENCE)


def depth_only_artifact_config(sv_cfg: dict | None) -> dict:
    """Resolve depth-only artifact QC knobs from ``config.sv``."""
    cfg = sv_cfg or {}
    return {
        "enabled": bool(cfg.get("drop_depth_only_artifacts", DEFAULT_DEPTH_ONLY_ARTIFACT["enabled"])),
        "max_rd_gq": float(cfg.get("depth_only_artifact_max_rd_gq", DEFAULT_DEPTH_ONLY_ARTIFACT["max_rd_gq"])),
    }


def pesr_caller_qc_config(sv_cfg: dict | None) -> dict:
    """Resolve raw PE/SR-caller (Manta) QC knobs from ``config.sv``."""
    cfg = sv_cfg or {}
    defaults = DEFAULT_PESR_CALLER_QC
    min_jq = cfg.get("min_junction_qual", defaults["min_junction_qual"])
    return {
        "enabled": bool(cfg.get("pesr_caller_qc", defaults["enabled"])),
        "max_span_bp": int(cfg.get("max_span_bp", defaults["max_span_bp"])),
        "require_cn_min_span_bp": int(
            cfg.get("require_cn_min_span_bp", defaults["require_cn_min_span_bp"])
        ),
        "drop_imprecise": bool(cfg.get("drop_imprecise", defaults["drop_imprecise"])),
        "large_event_min_span_bp": int(
            cfg.get("large_event_min_span_bp", defaults["large_event_min_span_bp"])
        ),
        "large_event_require_sr": bool(
            cfg.get("large_event_require_sr", defaults["large_event_require_sr"])
        ),
        "min_junction_qual": None if min_jq is None else float(min_jq),
        "max_hom_span_bp": int(cfg.get("max_hom_span_bp", defaults["max_hom_span_bp"])),
    }


def event_span_bp(info: dict, pos: int, end: int, svlen: int | None = None) -> int:
    """Positive genomic span for size filters (|SVLEN| preferred, else END-POS)."""
    length = svlen if svlen is not None else sv_length(info, pos, end)
    if length is not None and length != 0:
        return abs(int(length))
    return abs(int(end) - int(pos))


def manta_alt_count(fmt: dict, key: str = "SR") -> int | None:
    """Manta FORMAT ``SR``/``PR`` are ``ref,alt`` counts; return the alt count."""
    raw = fmt.get(key)
    if raw is None or raw is True or raw == "." or raw == "":
        return None
    parts = str(raw).split(",")
    if len(parts) < 2:
        return None
    try:
        return int(float(parts[1]))
    except ValueError:
        return None


def depth_copy_number(fmt: dict) -> int | None:
    """FORMAT ``CN``, else GATK-SV ``RD_CN``. Raw Manta priors have neither."""
    cn = _as_int(fmt.get("CN"))
    if cn is not None:
        return cn
    return _as_int(fmt.get("RD_CN"))


def lacks_depth_cn(fmt: dict) -> bool:
    """True when FORMAT lacks ``CN`` and ``RD_CN`` (raw Manta-like)."""
    return depth_copy_number(fmt) is None


def is_pesr_caller_artifact(
    info: dict,
    fmt: dict,
    *,
    svtype: str,
    span_bp: int,
    zygosity: str | None,
    max_span_bp: int = 2_000_000,
    require_cn_min_span_bp: int = 100_000,
    drop_imprecise: bool = True,
    large_event_min_span_bp: int = 50_000,
    large_event_require_sr: bool = True,
    min_junction_qual: float | None = 20.0,
    max_hom_span_bp: int = 50_000,
) -> bool:
    """True for raw PE/SR-caller DEL/DUP calls that should not enter triage.

    Applied only when FORMAT ``CN`` and ``RD_CN`` are both absent. GATK-SV
    records keep depth CN (``CN`` or ``RD_CN``) and skip this path. Rules
    (any match drops):

    - absolute span above ``max_span_bp``
    - span at or above ``require_cn_min_span_bp`` without depth ``CN``
    - ``INFO/IMPRECISE`` when ``drop_imprecise``
    - ``JUNCTION_QUAL`` present and below ``min_junction_qual``
    - span above ``large_event_min_span_bp`` without SR alt support
    - homozygous (or hemizygous) span above ``max_hom_span_bp``
    """
    if not lacks_depth_cn(fmt):
        return False
    if svtype not in _SPAN_FILTER_TYPES:
        return False
    if span_bp > max_span_bp:
        return True
    if require_cn_min_span_bp > 0 and span_bp >= require_cn_min_span_bp:
        return True
    if drop_imprecise and info.get("IMPRECISE") is True:
        return True
    if min_junction_qual is not None:
        jq = _as_float(info.get("JUNCTION_QUAL"))
        if jq is not None and jq < min_junction_qual:
            return True
    if large_event_require_sr and span_bp > large_event_min_span_bp:
        sr_alt = manta_alt_count(fmt, "SR")
        if sr_alt is None or sr_alt < 1:
            return True
    if zygosity in {"hom", "hemi"} and span_bp > max_hom_span_bp:
        return True
    return False


def is_depth_only_artifact(
    info: dict,
    fmt: dict,
    *,
    max_rd_gq: float = 5.0,
) -> bool:
    """True for depth-only DELs/CNVs that look like IGV coverage artifacts.

    Requires all of: ``ALGORITHMS=depth`` only, FORMAT ``EV`` without PE/SR,
    and ``RD_GQ`` or ``OGQ`` at or below ``max_rd_gq``. Cohort AF is not used.
    """
    if not is_depth_only_algorithms(info):
        return False
    if has_pesr_format_evidence(fmt):
        return False
    rd_gq = _as_float(fmt.get("RD_GQ"))
    ogq = _as_float(fmt.get("OGQ"))
    weak_rd = rd_gq is not None and rd_gq <= max_rd_gq
    weak_ogq = ogq is not None and ogq <= max_rd_gq
    if not (weak_rd or weak_ogq):
        return False
    return True


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
    span_bp: int = 0,
) -> bool:
    if rule.get("effect") == "copy_gain":
        require = bool(rule.get("require_copy_gain_or_exon", True))
        if require:
            copy_gain_or_exon = bool(predicted & GAIN_PREDICTED) or n_exons > 0
            if not copy_gain_or_exon:
                return False
        min_span = int(rule.get("min_span_bp") or 0)
        if min_span and int(span_bp or 0) < min_span:
            return False
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
    span_bp: int = 0,
) -> tuple[str, str, bool, str]:
    """Return mechanism, event, whether to score it, and the citation.

    ``scored`` is true for configured dosage genes and for exon-disrupting
    deletions. Copy gains are scored only for genes whose rule says ``copy_gain``
    (SNCA duplication/triplication) when the call looks like a gene
    multiplication (PREDICTED_COPY_GAIN / TSS / exon overlap and gene-scale span).
    LoF SVs use ``null_variant`` plus ``sv_dosage_lof``; copy-gain reuses
    ``null_variant`` as ``sv_dosage_gain``.
    """
    rules = rules if rules is not None else DEFAULT_DOSAGE_RULES
    event = dosage_event(svtype, cn, ecn, predicted)
    for rule in rules:
        if str(rule.get("gene", "")).upper() != gene:
            continue
        if _rule_matches(rule, svtype, event, predicted, n_exons, n_gene_exons, cn, ecn, span_bp):
            return str(rule["effect"]), event, True, str(rule.get("citation") or "")
    if predicted & LOF_PREDICTED or (svtype == "DEL" and n_exons > 0):
        return "lof", event if event == "deletion" else "deletion", True, "Exon disruption (GATK-SV PREDICTED_* or exon BED overlap)."
    if predicted & GAIN_PREDICTED:
        return "copy_gain", event, False, ""
    return "other", event, False, ""


def _consequence(mechanism: str, scored: bool, n_exons: int) -> str:
    if mechanism == "lof" and scored:
        return "transcript_ablation"
    if mechanism == "copy_gain" and scored:
        return "copy_number_gain"
    if n_exons > 0:
        return "structural_variant"
    return "intron_variant"


def _clinvar_raw(info: dict) -> str:
    """ClinVar CLNSIG from a panel-sliced INFO annotation, if present."""
    for key in ("ClinVar_CLNSIG", "CLNSIG", "clinvar_clnsig"):
        value = info.get(key)
        if value and value is not True:
            return str(value).split("|")[0]
    return ""


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


def sv_length(info: dict, pos: int, end: int) -> int | None:
    """INFO/SVLEN when the VCF has it, otherwise the POS-to-END span.

    The VCF sign is kept, so a deletion is negative when the caller followed
    the spec. A missing tag uses ``end - pos``.
    """
    tagged = _as_int(info.get("SVLEN"))
    if tagged is not None:
        return tagged
    if end != pos:
        return int(end) - int(pos)
    return None


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
    depth_artifact: dict | None = None,
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
    artifact_cfg = depth_only_artifact_config(depth_artifact)
    pesr_cfg = pesr_caller_qc_config(depth_artifact)
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
        svlen = sv_length(info, record["pos"], end)
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
            cn = depth_copy_number(fmt)
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
            span_bp = event_span_bp(info, record["pos"], end, svlen)
            artifact = bool(
                artifact_cfg["enabled"]
                and is_depth_only_artifact(
                    info,
                    fmt,
                    max_rd_gq=artifact_cfg["max_rd_gq"],
                )
            )
            if pesr_cfg["enabled"] and is_pesr_caller_artifact(
                info,
                fmt,
                svtype=svtype,
                span_bp=span_bp,
                zygosity=zyg,
                max_span_bp=pesr_cfg["max_span_bp"],
                require_cn_min_span_bp=pesr_cfg["require_cn_min_span_bp"],
                drop_imprecise=pesr_cfg["drop_imprecise"],
                large_event_min_span_bp=pesr_cfg["large_event_min_span_bp"],
                large_event_require_sr=pesr_cfg["large_event_require_sr"],
                min_junction_qual=pesr_cfg["min_junction_qual"],
                max_hom_span_bp=pesr_cfg["max_hom_span_bp"],
            ):
                continue
            algorithms = str(info.get("ALGORITHMS") or "") if info.get("ALGORITHMS") not in {True, None} else ""
            evidence = str(fmt.get("EV") or "")
            rd_gq = _as_float(fmt.get("RD_GQ"))
            ogq = _as_float(fmt.get("OGQ"))
            high_sr = info.get("HIGH_SR_BACKGROUND") is True or str(info.get("HIGH_SR_BACKGROUND") or "").upper() in {
                "TRUE",
                "YES",
                "1",
            }
            freq_het = _as_float(info.get("FREQ_HET"))
            clinvar_raw = _clinvar_raw(info)
            for gene in sorted(gene_symbols):
                predicted = effects.get(gene, set())
                n_exons = len({(h["chrom"], h["start"], h["end"]) for h in exons_by_gene.get(gene, [])})
                mechanism, event, scored, citation = classify_dosage(
                    gene, svtype, cn, ecn if ecn is not None else 2, predicted,
                    n_exons, exon_counts.get(gene, 0), rules, span_bp,
                )
                if artifact:
                    mechanism, scored, citation = "other", False, ""
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
                        # INFO AF / CN_NONREF_FREQ is the callset (cohort) frequency on
                        # joint or cohort-genotyped SV VCFs — annotate only. Do not store
                        # it as sv_panel_af; that field is reserved for a true GATK-SV
                        # reference-panel AF when available, and hard frequency filters
                        # must not treat in-cohort AF as population frequency.
                        "sv_panel_af": pd.NA,
                        "cohort_af": cohort_af,
                        "gnomad_af": gnomad_af,
                        "gnomad_popmax_af": popmax,
                        "cadd": pd.NA,
                        "revel": pd.NA,
                        "sift": pd.NA,
                        "sift_pred": "",
                        "polyphen": pd.NA,
                        "polyphen_pred": "",
                        "eigen": pd.NA,
                        "muttaster": pd.NA,
                        "muttaster_pred": "",
                        "fathmm": pd.NA,
                        "fathmm_pred": "",
                        "clinvar": clinvar_raw,
                        "clinvar_class": "none",
                        "variant_class": "SV",
                        "hit_source": "sv",
                        "svtype": svtype,
                        "svlen": svlen if svlen is not None else pd.NA,
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
                        "sv_algorithms": algorithms,
                        "sv_evidence": evidence,
                        "rd_gq": rd_gq if rd_gq is not None else pd.NA,
                        "ogq": ogq if ogq is not None else pd.NA,
                        "high_sr_background": high_sr,
                        "freq_het": freq_het if freq_het is not None else pd.NA,
                        "sv_artifact": artifact,
                    }
                )
    table = pd.DataFrame(rows) if rows else empty_variant_table()
    if not table.empty:
        table["clinvar_class"] = normalize_clinvar(table["clinvar"])
        table = derive_zygosity(table)
        # derive_zygosity keeps a provided zygosity. CN-inferred genotypes are already set.
    carriers = table["sample_id"].astype(str).unique().tolist() if not table.empty else []
    assessment = coverage_table(snv_samples, header_samples, carriers)
    return table, assessment


def apply_dosage_evidence(df: pd.DataFrame, scfg: dict) -> pd.DataFrame:
    """Adjust evidence for scored dosage SVs.

    LoF SVs already get ``null_variant`` via ``transcript_ablation``. They also
    receive ``sv_dosage_lof`` bonus points (SVs lack ClinVar/CADD/REVEL/AVI).
    Copy-gain events on configured genes (SNCA) receive ``null_variant`` points
    recorded as ``sv_dosage_gain``.
    """
    if df.empty or "dosage_mechanism" not in df.columns or "evidence_score" not in df.columns:
        return df
    df = df.copy()
    points = scfg["points"]
    null_pts = int(points["null_variant"])
    lof_bonus = int(points.get("sv_dosage_lof") or 0)
    reasons = df["evidence_reasons"].fillna("").astype(str)
    if "scored_dosage" in df.columns:
        scored = df["scored_dosage"].fillna(False).astype(bool)
    else:
        scored = df["dosage_mechanism"].isin(["lof", "copy_gain"])
    gain = scored & (df["dosage_mechanism"] == "copy_gain") & ~reasons.str.contains(
        "null_variant|sv_dosage", regex=True
    )
    df.loc[gain, "evidence_score"] = df.loc[gain, "evidence_score"] + null_pts
    df.loc[gain, "evidence_reasons"] = (reasons.loc[gain] + ";sv_dosage_gain").str.strip(";")
    reasons = df["evidence_reasons"].fillna("").astype(str)
    lof = scored & (df["dosage_mechanism"] == "lof") & ~reasons.str.contains("sv_dosage_lof", regex=False)
    if lof_bonus:
        df.loc[lof, "evidence_score"] = df.loc[lof, "evidence_score"] + lof_bonus
    df.loc[lof, "evidence_reasons"] = (reasons.loc[lof] + ";sv_dosage_lof").str.strip(";")
    df["tier"] = assign_tier(df["evidence_score"], scfg["tiers"])
    return df

