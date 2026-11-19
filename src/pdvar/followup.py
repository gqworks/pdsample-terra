"""Second-hit rescue, empty-shortlist module stubs, and early-onset keep rules."""

from __future__ import annotations

import pandas as pd

from .compound import (
    GAUCHIAN_NEGATIVE,
    NO_EXON_SV,
    NO_SECOND_SMALL,
    SINGLE_HET_REASON,
    _append_reason,
    _consequence_terms,
    _is_lof_sv,
    _text,
)
from .scoring import NULL_CONSEQUENCES

EARLY_ONSET_REASON = "early-onset recessive null below the tier floor; review for a second allele"
MODULE_STUB_STATUS = "module_not_assessed"

REQUIRED_MODULES = (
    ("gauchian", "gauchian not assessed"),
    ("expansionhunter", "expansionhunter not assessed"),
    ("mtdna", "mito not assessed"),
    ("introme", "introme not assessed"),
)

GBA1_PD_CLASSES = frozenset({"severe", "mild", "risk"})


def _numeric(value: object) -> float | None:
    if value is None or (isinstance(value, float) and pd.isna(value)):
        return None
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def _is_artifact(row: pd.Series) -> bool:
    raw = row.get("sv_artifact")
    if isinstance(raw, bool):
        return raw
    return _text(raw).lower() in {"true", "1", "yes"}


def _exon_sv(row: pd.Series) -> bool:
    if _is_lof_sv(row):
        return True
    try:
        return float(row.get("exons_overlapped") or 0) > 0
    except (TypeError, ValueError):
        return False


def _is_null(row: pd.Series) -> bool:
    return bool(_consequence_terms(row) & NULL_CONSEQUENCES)


def annotate_empty_searches(df: pd.DataFrame, *, gauchian_catalog_negative: bool = False) -> pd.DataFrame:
    """Append which second-hit checks were empty onto ``followup_reason``."""
    if df.empty or "followup_reason" not in df.columns:
        return df
    df = df.copy()
    follow = df["followup_reason"].fillna("").astype(str).str.strip().ne("")
    if not follow.any():
        return df
    for (_, _), grp in df[follow].groupby(["sample_id", "gene"], sort=False):
        has_second_small = False
        has_exon_sv = False
        statuses = set(grp["inheritance_status"].fillna("").astype(str))
        qualifying_small = 0
        for i in grp.index:
            row = grp.loc[i] if i in grp.index else df.loc[i]
            if _text(row.get("variant_class")) != "SV" and _text(row.get("hit_source")) != "sv":
                if _text(row.get("inheritance_status")) not in {"", "not_qualifying", "nan", "AR_intronic_candidate"}:
                    qualifying_small += 1
                if _is_null(df.loc[i]) or _text(df.loc[i].get("clinvar_class")) == "P/LP":
                    has_second_small = True
            if _exon_sv(df.loc[i]):
                has_exon_sv = True
        extras: list[str] = []
        if qualifying_small < 2 and not has_second_small:
            extras.append(NO_SECOND_SMALL)
        if not has_exon_sv:
            extras.append(NO_EXON_SV)
        gene = str(grp["gene"].iloc[0] or "").upper()
        if gene == "GBA1" and gauchian_catalog_negative:
            extras.append(GAUCHIAN_NEGATIVE)
        if not extras:
            continue
        for i in grp.index:
            reason = _text(df.at[i, "followup_reason"])
            for extra in extras:
                reason = _append_reason(reason, extra)
            df.at[i, "followup_reason"] = reason
        del statuses
    return df


def rescue_second_hits(
    scored: pd.DataFrame,
    unfiltered: pd.DataFrame | None,
    *,
    introme_variants: pd.DataFrame | None = None,
    cfg: dict | None = None,
    sample_id: str | None = None,
) -> pd.DataFrame:
    """Pull filter-dropped second-hit candidates in genes that already have followup_reason.

    Rescued rows keep ``followup_reason`` and are not auto-promoted to a solve.
    Introme injection is gene-scoped: high-score alleles in a follow-up gene that
    the caller missed.
    """
    if scored is None or scored.empty or "followup_reason" not in scored.columns:
        return scored
    follow_genes = {
        str(gene).upper()
        for gene in scored.loc[
            scored["followup_reason"].fillna("").astype(str).str.strip().ne(""),
            "gene",
        ].astype(str)
        if str(gene).strip()
    }
    if not follow_genes:
        return scored

    fcfg = (cfg or {}).get("filters") or {}
    min_depth = float(fcfg.get("second_hit_min_depth", 1))
    min_gq = float(fcfg.get("second_hit_min_gq", 20))
    scored_ids = set(scored["variant_id"].astype(str)) if "variant_id" in scored.columns else set()
    extras: list[pd.Series] = []

    def _keep_unfiltered(row: pd.Series) -> bool:
        gene = str(row.get("gene") or "").upper()
        if gene not in follow_genes:
            return False
        vid = str(row.get("variant_id") or "")
        if vid and vid in scored_ids:
            return False
        gq = _numeric(row.get("gq"))
        if gq is not None and gq < min_gq:
            return False
        if _is_null(row):
            depth = _numeric(row.get("depth"))
            if depth is not None and depth < min_depth:
                return False
            return True
        if _is_artifact(row) and _exon_sv(row):
            return True
        if _text(row.get("clinvar_class")) == "P/LP" and gene in follow_genes:
            return True
        return False

    if unfiltered is not None and not unfiltered.empty:
        for _, row in unfiltered.iterrows():
            if _keep_unfiltered(row):
                extras.append(row)

    if introme_variants is not None and not introme_variants.empty:
        spcfg = (cfg or {}).get("splicing") or {}
        high = float(spcfg.get("introme_high", 0.9))
        score_col = "introme_score" if "introme_score" in introme_variants.columns else None
        for _, row in introme_variants.iterrows():
            vid = str(row.get("variant_id") or "")
            if not vid or vid in scored_ids:
                continue
            genes = str(row.get("introme_genes") or row.get("gene") or "")
            tokens = {part.strip().upper() for part in genes.replace(",", ";").split(";") if part.strip()}
            if not tokens & follow_genes:
                continue
            score = _numeric(row.get(score_col) if score_col else None)
            if score is None or score < high:
                continue
            extras.append(_introme_row(row, tokens & follow_genes, sample_id))

    if not extras:
        return scored

    extra_df = pd.DataFrame(extras).copy()

    def _as_text(column: str, default: str = "") -> pd.Series:
        if column not in extra_df.columns:
            return pd.Series(default, index=extra_df.index, dtype="object")
        return extra_df[column].astype("string").fillna(default).astype(str).replace({"nan": default, "<NA>": default, "None": default})

    extra_df["followup_reason"] = _as_text("followup_reason").map(
        lambda value: _append_reason(value or SINGLE_HET_REASON, "rescued second-hit search")
    )
    status = _as_text("inheritance_status")
    status = status.where(~status.isin(["", "nan", "<NA>"]), "not_qualifying")
    extra_df["inheritance_status"] = status
    extra_df["evidence_score"] = pd.to_numeric(_as_text("evidence_score", "0"), errors="coerce").fillna(0)
    extra_df["tier"] = pd.to_numeric(_as_text("tier", "4"), errors="coerce").fillna(4).astype(int)
    extra_df["evidence_reasons"] = _as_text("evidence_reasons")
    if sample_id:
        sid = _as_text("sample_id")
        extra_df["sample_id"] = sid.where(~sid.isin(["", "nan", "<NA>"]), sample_id)
    if "variant_id" in extra_df.columns:
        extra_df = extra_df[~extra_df["variant_id"].astype(str).isin(scored_ids)]
    if extra_df.empty:
        return scored
    return pd.concat([scored, extra_df], ignore_index=True, sort=False)


def _introme_row(row: pd.Series, genes: set[str], sample_id: str | None) -> pd.Series:
    vid = str(row.get("variant_id") or "")
    chrom, pos, ref, alt = "", pd.NA, "", ""
    parts = vid.split("-")
    if len(parts) >= 4:
        chrom, pos_text, ref, alt = parts[0], parts[1], parts[2], "-".join(parts[3:])
        try:
            pos = int(pos_text)
        except ValueError:
            pos = pd.NA
    payload = {
        "sample_id": sample_id or "",
        "family_id": sample_id or "",
        "chrom": chrom,
        "pos": pos,
        "ref": ref,
        "alt": alt,
        "variant_id": vid,
        "gene": sorted(genes)[0] if genes else "",
        "consequence": "intron_variant",
        "clinvar_class": "none",
        "variant_class": "SNV" if ref and alt and len(str(ref)) == 1 and len(str(alt)) == 1 else "INDEL",
        "hit_source": "small_variant",
        "inheritance_status": "not_qualifying",
        "introme_score": row.get("introme_score"),
        "introme_class": "high",
        "zygosity": "het",
        "genotype": "0/1",
    }
    return pd.Series(payload)


def early_onset_keep(
    scored: pd.DataFrame,
    *,
    age_onset: object,
    classification: str,
    cfg: dict | None = None,
) -> pd.DataFrame:
    """Force-keep recessive nulls that missed the tier floor on early-onset unexplained samples."""
    if scored is None or scored.empty:
        return scored
    cfg = cfg or {}
    threshold = (cfg.get("phenotype") or {}).get("early_onset_age", 50)
    try:
        onset = float(age_onset)
    except (TypeError, ValueError):
        return scored
    if onset >= float(threshold):
        return scored
    statuses = scored["inheritance_status"].fillna("").astype(str) if "inheritance_status" in scored.columns else pd.Series("", index=scored.index)
    mendelian = statuses.isin(
        {
            "AR_hom",
            "AR_comphet_candidate",
            "AR_hemizygous",
            "AR_second_hit_splice",
            "AD_candidate",
            "XL_hemi",
        }
    )
    if classification == "solved" or mendelian.any():
        return scored
    if classification not in {"unexplained", "partially_explained"}:
        return scored
    genes = {
        str(item).upper()
        for item in ((cfg.get("biomarkers") or {}).get("gene_classes") or {}).get("recessive_eopd") or []
    }
    if not genes:
        return scored
    gene = scored["gene"].astype(str).str.upper() if "gene" in scored.columns else pd.Series("", index=scored.index)
    lof_sv = scored.apply(_is_lof_sv, axis=1) if not scored.empty else pd.Series(False, index=scored.index)
    nulls = scored.apply(_is_null, axis=1) if "consequence" in scored.columns else pd.Series(False, index=scored.index)
    target = gene.isin(genes) & (nulls | lof_sv) & statuses.eq("not_qualifying")
    if not target.any():
        return scored
    scored = scored.copy()
    if "followup_reason" not in scored.columns:
        scored["followup_reason"] = ""
    for i in scored.index[target]:
        scored.at[i, "followup_reason"] = _append_reason(scored.at[i, "followup_reason"], EARLY_ONSET_REASON)
    return scored


def module_stub_rows(
    sample_id: str,
    notes: dict,
    *,
    sv_status: str = "",
    avi_requested: bool = False,
    pdvar_git_commit: str = "",
) -> pd.DataFrame:
    """One ``candidates.tsv`` row per required module that was not assessed."""
    rows: list[dict] = []
    commit = str(pdvar_git_commit or "").strip()

    def _row(module: str, reason: str) -> dict:
        text = reason
        if commit:
            text = f"{text}; pdvar_git_commit {commit}"
        return {
            "rank": len(rows) + 1,
            "sample_id": sample_id,
            "gene": module,
            "moi": "",
            "inheritance_status": MODULE_STUB_STATUS,
            "variant_id": f"module-{module}",
            "followup_reason": text,
            "hit_source": "module",
            "evidence_score": 0,
            "tier": 4,
            "consequence": "",
        }

    for key, reason in REQUIRED_MODULES:
        if str(notes.get(key) or "") == "not_assessed":
            rows.append(_row(key, reason))
    gatk = str(notes.get("gatk_sv") or sv_status or "")
    if gatk in {"", "not_assessed"}:
        rows.append(_row("sv", "sv not assessed"))
    if str(sv_status or notes.get("gatk_sv") or "") == "prior":
        rows.append(_row("sv", "sv_status prior from an older callset"))
    if str(notes.get("introme") or "") == "prior":
        rows.append(_row("introme", "introme prior from an older callset"))
    if avi_requested and str(notes.get("avi") or "not_assessed") != "called":
        rows.append(_row("avi", "avi not assessed"))
    if not rows:
        return pd.DataFrame()
    for index, row in enumerate(rows, start=1):
        row["rank"] = index
    return pd.DataFrame(rows)
