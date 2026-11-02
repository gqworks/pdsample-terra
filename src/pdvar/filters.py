"""Panel, frequency, quality and consequence filters."""

from __future__ import annotations

import math

import pandas as pd

from .incidental import ttn_truncation_mask
from .vcf_stream import chrom_bare


# Null / splice-disrupting terms eligible for near-miss depth rescue.
DEFAULT_DEPTH_RESCUE_HIGH_IMPACT = (
    "transcript_ablation",
    "splice_acceptor_variant",
    "splice_donor_variant",
    "stop_gained",
    "frameshift_variant",
    "stop_lost",
    "start_lost",
)


def _le_or_missing(series: pd.Series, threshold: float) -> pd.Series:
    """True where value <= threshold or missing (absent from gnomAD counts as rare)."""
    return series.isna() | (series <= threshold)


def _ge_or_missing(series: pd.Series, threshold: float) -> pd.Series:
    return series.isna() | (series >= threshold)


def panel_mask(df: pd.DataFrame, panel: pd.DataFrame) -> pd.Series:
    return df["gene"].astype(str).str.upper().isin(set(panel["gene"]))


def frequency_mask(df: pd.DataFrame, fcfg: dict) -> pd.Series:
    """Population frequency ceilings (gnomAD; optional true GATK-SV ref-panel AF).

    ``cohort_af`` (in-callset AF on a joint SV VCF) is never applied here — a
    pathogenic SV that happens to recur in a small PD cohort must not be dropped.
    ``sv_panel_af`` is filtered only when present and is reserved for reference-panel
    AF, not cohort INFO/AF.
    """
    mask = pd.Series(True, index=df.index)
    for col, key in (
        ("gnomad_af", "max_gnomad_af"),
        ("gnomad_popmax_af", "max_gnomad_popmax_af"),
        ("sv_panel_af", "max_gnomad_af"),
    ):
        if col in df.columns and key in fcfg:
            mask &= _le_or_missing(df[col], fcfg[key])
    return mask


def _gq_ab_mask(df: pd.DataFrame, fcfg: dict) -> pd.Series:
    """Quality gates other than depth (GQ and het allele balance)."""
    mask = pd.Series(True, index=df.index)
    if "gq" in df.columns and "min_gq" in fcfg:
        mask &= _ge_or_missing(df["gq"], fcfg["min_gq"])
    if "allele_balance" in df.columns and "min_allele_balance" in fcfg:
        het = df.get("zygosity", pd.Series("het", index=df.index)) == "het"
        mask &= ~het | _ge_or_missing(df["allele_balance"], fcfg["min_allele_balance"])
    return mask


def quality_mask(df: pd.DataFrame, fcfg: dict) -> pd.Series:
    mask = _gq_ab_mask(df, fcfg)
    if "depth" in df.columns and "min_depth" in fcfg:
        mask &= _ge_or_missing(df["depth"], fcfg["min_depth"])
    return mask


def consequence_mask(df: pd.DataFrame, keep: list[str]) -> pd.Series:
    """True if any '&'/','-separated consequence term is in ``keep``."""
    keep_set = set(keep)
    terms = df["consequence"].fillna("").astype(str).str.split(r"[&,]")
    return terms.map(lambda ts: any(t.strip() in keep_set for t in ts)).astype(bool)


def _is_lof_deletion(df: pd.DataFrame) -> pd.Series:
    """Exon-disrupting / LoF deletions that can hemizygose an overlapping SNV/indel."""
    svtype = df["svtype"].fillna("").astype(str) if "svtype" in df.columns else pd.Series("", index=df.index)
    event = df["dosage_event"].fillna("").astype(str) if "dosage_event" in df.columns else pd.Series("", index=df.index)
    mechanism = (
        df["dosage_mechanism"].fillna("").astype(str).str.lower()
        if "dosage_mechanism" in df.columns
        else pd.Series("", index=df.index)
    )
    consequence = df["consequence"].fillna("").astype(str) if "consequence" in df.columns else pd.Series("", index=df.index)
    scored = (
        df["scored_dosage"].fillna(False).astype(bool)
        if "scored_dosage" in df.columns
        else pd.Series(False, index=df.index)
    )
    exon_count = (
        pd.to_numeric(df["exon_count"], errors="coerce").fillna(0)
        if "exon_count" in df.columns
        else pd.Series(0, index=df.index)
    )
    hit_source = df["hit_source"].fillna("").astype(str) if "hit_source" in df.columns else pd.Series("", index=df.index)
    variant_class = (
        df["variant_class"].fillna("").astype(str) if "variant_class" in df.columns else pd.Series("", index=df.index)
    )
    is_sv = hit_source.eq("sv") | variant_class.eq("SV")
    is_del = svtype.eq("DEL") | event.eq("deletion")
    lof_signal = mechanism.eq("lof") | consequence.str.contains("transcript_ablation", regex=False) | (
        scored & (exon_count > 0)
    )
    return is_sv & is_del & lof_signal


def _interval_end(row: pd.Series) -> int | None:
    try:
        start = int(row["pos"])
    except (TypeError, ValueError, KeyError):
        return None
    try:
        if "end" in row and pd.notna(row["end"]):
            return max(start, int(row["end"]))
    except (TypeError, ValueError):
        pass
    return start


def _inside_lof_del_mask(df: pd.DataFrame, anchors: pd.DataFrame) -> pd.Series:
    """True for small-variant rows whose position falls inside a LoF DEL anchor."""
    out = pd.Series(False, index=df.index)
    if anchors.empty or "pos" not in df.columns or "chrom" not in df.columns:
        return out
    if "hit_source" in df.columns:
        small = df["hit_source"].fillna("").astype(str).eq("small_variant")
    elif "variant_class" in df.columns:
        small = ~df["variant_class"].fillna("").astype(str).eq("SV")
    else:
        small = pd.Series(True, index=df.index)

    for aidx, anchor in anchors.iterrows():
        a_chrom = chrom_bare(str(anchor.get("chrom") or ""))
        a_end = _interval_end(anchor)
        try:
            a_start = int(anchor["pos"])
        except (TypeError, ValueError, KeyError):
            continue
        if a_end is None or not a_chrom:
            continue
        a_gene = str(anchor.get("gene") or "").upper()
        a_sample = str(anchor.get("sample_id") or "")
        for idx in df.index[small & ~out]:
            row = df.loc[idx]
            if chrom_bare(str(row.get("chrom") or "")) != a_chrom:
                continue
            if a_sample and "sample_id" in row.index and str(row.get("sample_id") or "") not in {"", a_sample}:
                continue
            gene = str(row.get("gene") or "").upper()
            if a_gene and gene not in {"", ".", a_gene}:
                continue
            try:
                pos = int(row["pos"])
            except (TypeError, ValueError, KeyError):
                continue
            if a_start <= pos <= a_end:
                out.loc[idx] = True
    return out


def _depth_rescue_mask(df: pd.DataFrame, fcfg: dict) -> tuple[pd.Series, pd.Series]:
    """Rescue depth-only QC failures for HIGH-impact near-miss or LoF-DEL overlap.

    Returns ``(rescued, reason)`` where reason is ``near_miss``, ``lof_del``, both
    joined with ``+``, or empty.
    """
    false = pd.Series(False, index=df.index)
    reason = pd.Series("", index=df.index, dtype=object)
    if "depth" not in df.columns or "min_depth" not in fcfg:
        return false, reason

    min_depth = float(fcfg["min_depth"])
    depth = pd.to_numeric(df["depth"], errors="coerce")
    depth_fail = depth.notna() & (depth < min_depth)
    if not depth_fail.any():
        return false, reason

    non_depth_ok = _gq_ab_mask(df, fcfg)
    eligible = depth_fail & non_depth_ok

    near = false.copy()
    factor = fcfg.get("depth_rescue_near_miss_factor")
    if factor is not None:
        try:
            factor_f = float(factor)
        except (TypeError, ValueError):
            factor_f = 0.0
        if factor_f > 1.0:
            floor = min_depth / factor_f
            high_terms = list(fcfg.get("depth_rescue_high_impact_consequences") or DEFAULT_DEPTH_RESCUE_HIGH_IMPACT)
            high = consequence_mask(df, high_terms) if "consequence" in df.columns else false
            near = eligible & high & (depth >= floor)

    inside = false.copy()
    lof_floor = fcfg.get("depth_rescue_lof_del_min_depth")
    if lof_floor is not None:
        try:
            lof_floor_f = float(lof_floor)
        except (TypeError, ValueError):
            lof_floor_f = math.nan
        if not math.isnan(lof_floor_f):
            # Anchors must already clear panel/frequency/quality/consequence.
            ones = pd.Series(True, index=df.index)
            anchor_ok = (
                df["pass_panel"] if "pass_panel" in df.columns else ones
            ) & (
                df["pass_frequency"] if "pass_frequency" in df.columns else ones
            ) & (
                df["pass_quality"] if "pass_quality" in df.columns else ones
            ) & (
                df["pass_consequence"] if "pass_consequence" in df.columns else ones
            ) & _is_lof_deletion(df)
            anchors = df.loc[anchor_ok]
            inside = eligible & (depth >= lof_floor_f) & _inside_lof_del_mask(df, anchors)

    rescued = near | inside
    reason.loc[near & ~inside] = "near_miss"
    reason.loc[inside & ~near] = "lof_del"
    reason.loc[near & inside] = "near_miss+lof_del"
    return rescued, reason


def apply_filters(df: pd.DataFrame, panel: pd.DataFrame, fcfg: dict) -> pd.DataFrame:
    """Apply all filters and record which ones each retained variant passed.

    ClinVar P/LP variants are rescued from frequency/consequence filters when
    ``always_keep_clinvar_plp`` is set (they still need to be in the panel and pass QC).
    Panel genes flagged ``risk_allele_exception`` use the relaxed ``risk_exception_*`` AF ceilings.

    Depth-only failures can be rescued when a HIGH-impact panel allele is within
    ``depth_rescue_near_miss_factor`` of ``min_depth``, or when the allele sits
    inside a same-sample LoF deletion that itself passes filters (hemizygous
    coverage drop).
    """
    df = df.copy()
    df["pass_panel"] = panel_mask(df, panel) if fcfg.get("restrict_to_panel", True) else True
    df["pass_frequency"] = frequency_mask(df, fcfg)
    df["pass_quality"] = quality_mask(df, fcfg)
    df["pass_consequence"] = consequence_mask(df, fcfg["keep_consequences"])

    if "risk_allele_exception" in panel.columns and fcfg.get("risk_exception_max_af") is not None:
        exc_genes = set(panel.loc[panel["risk_allele_exception"].astype(bool), "gene"])
        relaxed = {
            **fcfg,
            "max_gnomad_af": fcfg["risk_exception_max_af"],
            "max_gnomad_popmax_af": fcfg["risk_exception_max_af"],
        }
        exc = df["gene"].astype(str).str.upper().isin(exc_genes) & frequency_mask(df, relaxed)
        df["risk_exception_rescued"] = exc & ~df["pass_frequency"]
        df["pass_frequency"] |= exc

    keep = df["pass_panel"] & df["pass_quality"] & df["pass_frequency"] & df["pass_consequence"]
    sources = fcfg.get("always_keep_sources") or []
    if sources and "hit_source" in df.columns:
        source_keep = df["hit_source"].isin(sources) & df["pass_quality"]
        df["source_rescued"] = source_keep & ~keep
        keep |= source_keep
    if fcfg.get("always_keep_clinvar_plp", True) and "clinvar_class" in df.columns:
        rescued = df["pass_panel"] & df["pass_quality"] & (df["clinvar_class"] == "P/LP")
        df["clinvar_rescued"] = rescued & ~keep
        keep |= rescued
    ttn = ttn_truncation_mask(df)
    df["ttn_truncation"] = ttn
    keep |= ttn

    depth_ok, depth_reason = _depth_rescue_mask(df, fcfg)
    depth_rescue = (
        df["pass_panel"]
        & df["pass_frequency"]
        & df["pass_consequence"]
        & ~df["pass_quality"]
        & depth_ok
    )
    df["depth_rescued"] = depth_rescue & ~keep
    df["depth_rescue_reason"] = ""
    df.loc[df["depth_rescued"], "depth_rescue_reason"] = depth_reason.loc[df["depth_rescued"]]
    keep |= depth_rescue

    return df.loc[keep].reset_index(drop=True)
