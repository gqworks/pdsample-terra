"""Per-sample candidate lists and a short HTML/Markdown report."""

from __future__ import annotations

import html
import json
from pathlib import Path

import pandas as pd

STATUS_PRIORITY = {
    "AR_hom": 0,
    "AR_comphet_candidate": 0,
    "AR_hemizygous": 0,
    "AR_second_hit_splice": 0,
    # AD / XL / risk compete with AR solves so ClinVar P/LP alleles (e.g. LRRK2
    # R1441C) are not buried under intron AR Introme noise or depth-artifact SVs.
    # Intron-only / Introme AD_candidate rows are bumped to priority 1 in
    # rank_candidates so they do not share the Mendelian shortlist band.
    "AD_candidate": 0,
    "XL_hemi": 0,
    "risk_biallelic": 0,
    "risk_carrier": 0,
    "mtdna_pathogenic_candidate": 0,
    "AR_single_het": 3,
    "review_conflict": 3,
    # Haplogroup / unannotated mtDNA stays off the shortlist (see patient_candidates).
    "mtdna_variant": 8,
    "str_expanded": 1,
    "MOI_unknown_candidate": 4,
    # Introme-only intron support demoted below coding / SV AR solves.
    "AR_intronic_candidate": 5,
    # GoF / DN AD genes: intron splice LoF is not the disease mechanism.
    "AD_gof_intronic": 6,
    "not_qualifying": 9,
    "module_not_assessed": 8,
}

FOLLOWUP_STATUSES = frozenset({"AR_single_het", "review_conflict"})

# Genomic id, then the c./p. HGVS and SV length, so the report table reads
# coordinates and annotations together.
LEAD_COLUMNS = [
    "rank", "sample_id", "gene", "moi", "inheritance_status", "followup_reason",
    "variant_id", "chrom", "pos", "ref", "alt", "end",
    "hgvsc", "hgvsp", "svtype", "svlen",
]


def _intronic_rank_penalty(df: pd.DataFrame) -> pd.Series:
    """Push intron-only Introme AR labels below coding/SV rows of the same status."""
    cons = df["consequence"].fillna("").astype(str) if "consequence" in df.columns else pd.Series("", index=df.index)
    clin = df["clinvar_class"].fillna("").astype(str) if "clinvar_class" in df.columns else pd.Series("", index=df.index)
    intron = cons.str.contains("intron_variant", regex=False)
    codingish = cons.str.contains(
        r"missense|frameshift|stop_|start_lost|splice_acceptor|splice_donor|transcript_ablation|"
        r"inframe_|protein_altering|structural_variant|copy_number",
        regex=True,
        na=False,
    )
    return (intron & ~codingish & (clin != "P/LP")).astype(int)


def _clinvar_plp_rank_boost(df: pd.DataFrame) -> pd.Series:
    """Prefer ClinVar P/LP coding alleles over same-priority intron/SV noise (0 beats 1)."""
    clin = df["clinvar_class"].fillna("").astype(str) if "clinvar_class" in df.columns else pd.Series("", index=df.index)
    cons = df["consequence"].fillna("").astype(str) if "consequence" in df.columns else pd.Series("", index=df.index)
    codingish = cons.str.contains(
        r"missense|frameshift|stop_|start_lost|splice_acceptor|splice_donor|"
        r"inframe_|protein_altering",
        regex=True,
        na=False,
    )
    return (~(clin.eq("P/LP") & codingish)).astype(int)


def _clinvar_vus_coding_rank_boost(df: pd.DataFrame) -> pd.Series:
    """Prefer ClinVar VUS coding alleles after P/LP (0 beats 1)."""
    from .inheritance import clinvar_vus_coding_mask

    return (~clinvar_vus_coding_mask(df)).astype(int)


def _clinvar_vus_noncoding_rank_boost(df: pd.DataFrame) -> pd.Series:
    """Prefer other ClinVar VUS after coding VUS (0 beats 1)."""
    from .inheritance import clinvar_vus_coding_mask, clinvar_vus_mask

    vus = clinvar_vus_mask(df)
    return (~(vus & ~clinvar_vus_coding_mask(df))).astype(int)


def rank_candidates(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    df["status_priority"] = df["inheritance_status"].map(STATUS_PRIORITY).fillna(9)
    df["_intronic_penalty"] = _intronic_rank_penalty(df)
    # Introme / intron-only AD_candidate rows share priority 0 with coding Mendelian
    # hits and flood the shortlist; demote them to priority 1 (with STR expansions).
    ad = df["inheritance_status"].fillna("").astype(str).eq("AD_candidate")
    df.loc[ad & df["_intronic_penalty"].eq(1), "status_priority"] = 1
    df["_plp_boost"] = _clinvar_plp_rank_boost(df)
    df["_vus_coding_boost"] = _clinvar_vus_coding_rank_boost(df)
    df["_vus_noncoding_boost"] = _clinvar_vus_noncoding_rank_boost(df)
    df["_rank_score"] = df["evidence_score"]
    # ClinVar order: P/LP coding, then VUS coding, then other VUS, then the rest.
    # plp_boost before intronic_penalty so LRRK2 R1441C-style hits beat intron AR and
    # transcript_ablation SVs that share status_priority 0.
    sort_cols = [
        "sample_id",
        "status_priority",
        "_plp_boost",
        "_vus_coding_boost",
        "_vus_noncoding_boost",
        "_intronic_penalty",
        "tier",
        "_rank_score",
    ]
    asc = [True, True, True, True, True, True, True, False]
    drop_cols = ["_rank_score", "_intronic_penalty", "_plp_boost", "_vus_coding_boost", "_vus_noncoding_boost"]
    if "avi_phred" in df.columns:
        df["_avi"] = pd.to_numeric(df["avi_phred"], errors="coerce").fillna(-1)
        sort_cols.append("_avi")
        asc.append(False)
        drop_cols.append("_avi")
    if "phenotype_score" in df.columns:
        sort_cols.append("phenotype_score")
        asc.append(False)
    df = df.sort_values(sort_cols, ascending=asc, kind="mergesort")
    df["rank"] = df.groupby("sample_id").cumcount() + 1
    return df.drop(columns=drop_cols)


def _must_keep_mask(df: pd.DataFrame) -> pd.Series:
    """Rows that must not be crowded out of the per-sample shortlist.

    ClinVar P/LP coding alleles, Gauchian GBA1 calls, risk statuses, scored
    LoF SVs, and AR follow-up genes (single-het / review_conflict) are the
    findings a rereview must see even when Introme fills the top-N with intron
    AR labels. Follow-up alleles sort after Mendelian solves via
    ``STATUS_PRIORITY``.
    """
    status = df["inheritance_status"].fillna("").astype(str) if "inheritance_status" in df.columns else pd.Series("", index=df.index)
    clin = df["clinvar_class"].fillna("").astype(str) if "clinvar_class" in df.columns else pd.Series("", index=df.index)
    hit = df["hit_source"].fillna("").astype(str) if "hit_source" in df.columns else pd.Series("", index=df.index)
    mech = (
        df["dosage_mechanism"].fillna("").astype(str).str.lower()
        if "dosage_mechanism" in df.columns
        else pd.Series("", index=df.index)
    )
    scored = (
        df["scored_dosage"].fillna(False).astype(bool)
        if "scored_dosage" in df.columns
        else pd.Series(False, index=df.index)
    )
    follow = (
        df["followup_reason"].fillna("").astype(str).str.strip().ne("")
        if "followup_reason" in df.columns
        else pd.Series(False, index=df.index)
    )
    force_status = status.isin(
        {
            "AR_comphet_candidate",
            "AR_hemizygous",
            "AR_second_hit_splice",
            "str_expanded",
            "risk_carrier",
            "risk_biallelic",
            "AD_candidate",
            "XL_hemi",
            # High-evidence mtDNA only — common mtdna_variant haplogroup rows
            # must stay under the top-N cap (Bugbot #58).
            "mtdna_pathogenic_candidate",
            # Single-het AR / review_conflict: keep for same-gene second-hit review.
            *FOLLOWUP_STATUSES,
        }
    )
    # Introme / intron-only AD_candidate noise is ranked at priority 1 and must
    # not force-keep past the top-N cut (coding / P/LP AD still force-kept).
    intron_ad = status.eq("AD_candidate") & _intronic_rank_penalty(df).eq(1)
    mechanism = (
        df["mechanism"].fillna("").astype(str).str.lower()
        if "mechanism" in df.columns
        else pd.Series("", index=df.index)
    )
    # MAPT exon-10 splice_gof stays visible even though it is intron-like.
    splice_gof_ad = status.eq("AD_candidate") & mechanism.eq("splice_gof")
    force_status = (force_status & ~intron_ad) | splice_gof_ad
    plp = clin.eq("P/LP")
    gauchian = hit.eq("gauchian")
    lof_sv = hit.eq("sv") & (mech.eq("lof") | scored)
    return force_status | plp | gauchian | lof_sv | follow


def order_candidate_columns(df: pd.DataFrame) -> pd.DataFrame:
    """Put the genomic id, c./p. HGVS, and SV length at the front of the table."""
    if df is None or df.empty:
        return df
    lead = [column for column in LEAD_COLUMNS if column in df.columns]
    rest = [column for column in df.columns if column not in lead]
    return df.loc[:, lead + rest]


def patient_candidates(df: pd.DataFrame, max_per_patient: int) -> pd.DataFrame:
    """Top-N per sample among qualifying rows, plus forced clinical keepers.

    Generic ``mtdna_variant`` rows (score 0 / ClinVar B/LB haplogroup alleles) are
    not shortlisted. Only ``mtdna_pathogenic_candidate`` mtDNA appears, via the
    must-keep mask.
    """
    ranked = rank_candidates(df)
    status = ranked["inheritance_status"].fillna("").astype(str)
    follow = (
        ranked["followup_reason"].fillna("").astype(str).str.strip().ne("")
        if "followup_reason" in ranked.columns
        else pd.Series(False, index=ranked.index)
    )
    reportable = ranked[~status.eq("mtdna_variant") & (~status.eq("not_qualifying") | follow)].copy()
    if reportable.empty:
        return reportable
    # Re-rank after dropping haplogroup mtDNA so rank<=N is among reportable rows.
    # not_qualifying rows with followup_reason stay available for force-keep.
    qualifying = reportable[~reportable["inheritance_status"].fillna("").astype(str).eq("not_qualifying")].copy()
    if qualifying.empty:
        qualifying = reportable.copy()
    qualifying["rank"] = qualifying.groupby("sample_id", sort=False).cumcount() + 1
    top = qualifying[qualifying["rank"] <= max_per_patient]
    must = reportable[_must_keep_mask(reportable)]
    out = pd.concat([top, must], ignore_index=False).drop_duplicates()
    out = out.sort_values(["sample_id", "rank"], kind="mergesort")
    out["rank"] = out.groupby("sample_id", sort=False).cumcount() + 1
    return out.reset_index(drop=True)


def _jsonable(value):
    if isinstance(value, dict):
        return {str(key): _jsonable(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_jsonable(item) for item in value]
    if value is None or isinstance(value, (str, int, bool)):
        return value
    if isinstance(value, float):
        if value != value:
            return None
        return value
    if hasattr(value, "item"):
        try:
            return _jsonable(value.item())
        except (ValueError, AttributeError):
            return str(value)
    try:
        if pd.isna(value):
            return None
    except (TypeError, ValueError):
        pass
    return str(value)


def render_markdown(sample_id: str, candidates: pd.DataFrame, liability: pd.DataFrame, qc: dict) -> str:
    """Narrative report. This is a triage list, not an ACMG classification."""
    classification = "unexplained"
    explained = None
    note = ""
    if liability is not None and not liability.empty:
        row = liability.iloc[0]
        classification = str(row.get("classification") or "unexplained")
        explained = row.get("explained_liability")
        note = str(row.get("explained_liability_note") or "")
    lines = [
        f"# {sample_id}",
        "",
        "Per-sample Parkinson's disease variant triage. Scores are ACMG-style points, not an ACMG/AMP classification.",
        "",
        f"Classification: **{classification}**.",
        f"Explained liability: {explained if explained not in (None, '') else 'null'}. {note}",
        "",
    ]
    lines.extend(["## Modules", ""])
    modules = qc.get("modules") or {}
    for name, status in modules.items():
        lines.append(f"- {name}: {status}")
    if qc.get("sv_status"):
        lines.append(f"- gatk_sv status field: {qc.get('sv_status')}")
    if qc.get("pdvar_version"):
        lines.append(f"- pdvar: {qc.get('pdvar_version')}")
    if qc.get("pdvar_git_commit"):
        lines.append(f"- pdvar_git_commit: {qc.get('pdvar_git_commit')}")
    lines.extend(["", "## Callsets", ""])
    callsets = qc.get("callsets") or []
    sv_rows = qc.get("sv") or []
    if not callsets and not sv_rows:
        lines.append("No STR, GBA1, mtDNA, or SV assessment rows.")
    else:
        for row in callsets:
            lines.append(f"- {row.get('module')}: {row.get('callset_status')}")
        for row in sv_rows:
            lines.append(f"- gatk_sv: {row.get('sv_status')}")
    lines.extend(_phenotype_lines(candidates, qc))
    lines.extend(["", "## Ranked candidates", ""])
    if candidates is None or candidates.empty:
        lines.append("No qualifying candidate variants.")
    else:
        # heteroplasmy / MITOMAP clinical disease text are blank on nuclear rows;
        # they matter for mtdna_pathogenic_candidate shortlist hits.
        show = [
            column
            for column in (
                "rank",
                "gene",
                "inheritance_status",
                "variant_id",
                "hgvsc",
                "hgvsp",
                "svtype",
                "svlen",
                "consequence",
                "zygosity",
                "heteroplasmy",
                "mitomap_clinical_info",
                "tier",
                "evidence_score",
                "splice_rule",
                "splice_subgroup",
                "splice_sav_pct",
                "splice_heuristic_n",
                "clinical_match",
                "phenotype_score",
                "avi_phred",
                "evidence_reasons",
            )
            if column in candidates.columns
        ]
        lines.append(candidates[show].to_csv(sep="\t", index=False).rstrip())
    lines.append("")
    return "\n".join(lines)


def _risk_rows(df: pd.DataFrame) -> pd.DataFrame:
    status = df["inheritance_status"].astype(str) if "inheritance_status" in df.columns else pd.Series("", index=df.index)
    moi = df["moi"].astype(str).str.upper() if "moi" in df.columns else pd.Series("", index=df.index)
    if "risk_allele_exception" in df.columns:
        flag = df["risk_allele_exception"]
        if flag.dtype == bool:
            exception = flag.fillna(False)
        else:
            exception = flag.astype(str).str.strip().str.lower().isin(["true", "1", "yes"])
    else:
        exception = pd.Series(False, index=df.index)
    mask = status.isin(["risk_carrier", "risk_biallelic"]) | moi.eq("RISK") | exception
    return df.loc[mask]


def _show(value: object) -> str:
    if value is None:
        return "blank"
    text = str(value).strip()
    return text or "blank"


def _phenotype_detail_lines(info: dict) -> list[str]:
    lines = [
        "Age at onset: {age}. Symptoms: {symptoms}. Ancestry: {ancestry}.".format(
            age=_show(info.get("age_onset")),
            symptoms=_show(info.get("symptoms")),
            ancestry=_show(info.get("ethnicity")),
        )
    ]
    if info.get("phenotype_source") == "report":
        lines.append(
            "Phenotype report {gp2} / {study}. Reported sex: {sex}. PD status: {pd}. "
            "Age at sample collection: {sampling}. Age at diagnosis: {diagnosis_age}. Diagnosis type: {diagnosis}. "
            "Levodopa: {levodopa}. Age at levodopa start: {levodopa_age}.".format(
                gp2=_show(info.get("gp2_id")),
                study=_show(info.get("study_id")),
                sex=_show(info.get("reported_sex")),
                pd=_show(info.get("pd_status")),
                sampling=_show(info.get("age_sampling")),
                diagnosis_age=_show(info.get("age_diagnosis")),
                diagnosis=_show(info.get("diagnosis_type")),
                levodopa=_show(info.get("levodopa")),
                levodopa_age=_show(info.get("age_levodopa_start")),
            )
        )
        lines.append(
            "Levodopa use is stored as a biomarker note on early-onset recessive genes. "
            "Gene-profile keywords still come from symptom text. "
            "The comorbidity list is recorded in this section."
        )
    grouped: dict[str, list[str]] = {}
    for item in info.get("comorbidities") or []:
        if not isinstance(item, dict):
            continue
        relevance = str(item.get("relevance") or "unspecified")
        label = str(item.get("condition") or "").strip()
        if item.get("age") is not None:
            label = f"{label} (age {item['age']})"
        grouped.setdefault(relevance, []).append(label)
    for relevance, labels in grouped.items():
        lines.append(f"Comorbidities ({relevance}): {'; '.join(labels)}.")
    return lines


def _phenotype_lines(candidates: pd.DataFrame, qc: dict) -> list[str]:
    info = qc.get("phenotype") or {}
    lines = ["", "## Phenotype concordance", ""]
    if not info.get("provided"):
        lines.append(
            "No phenotype report was supplied, so clinical concordance was not assessed. "
            "Provide a GP2-style phenotype_report text file on the sample table when the patient has one."
        )
        return lines
    lines.extend(_phenotype_detail_lines(info))
    lines.append("")
    lines.append(
        "Candidates whose gene has a profile in config/gene_phenotypes.yaml are compared with those fields. "
        "clinical_match is matched when onset is inside the profile range, a symptom keyword hits a key or atypical feature, "
        "or the ancestry is one the profile lists. unmatched means the onset sits outside that range and no feature or ancestry hit. "
        "uninformative means the gene has a profile and these fields neither confirm nor refute it. "
        "no_profile means that gene has no phenotype profile yet. The comparison is a triage note."
    )
    lines.append("")
    lines.append("Risk candidates:")
    lines.append("")
    if candidates is None or candidates.empty:
        lines.append("No ranked candidates.")
        return lines
    risk = _risk_rows(candidates)
    if risk.empty:
        lines.append("No risk-labelled candidates were ranked.")
        return lines
    show = [
        column
        for column in (
            "rank",
            "gene",
            "moi",
            "inheritance_status",
            "variant_id",
            "clinical_match",
            "onset_match",
            "matched_features",
            "matched_atypical",
            "ancestry_enriched",
            "phenotype_score",
        )
        if column in risk.columns
    ]
    lines.append(risk[show].to_csv(sep="\t", index=False).rstrip())
    return lines


def render_html(markdown: str, candidates: pd.DataFrame) -> str:
    table = "<p>No ranked candidates.</p>"
    if candidates is not None and not candidates.empty:
        table = candidates.to_html(index=False, border=0, escape=True)
    body = html.escape(markdown)
    return (
        "<!DOCTYPE html><html><head><meta charset=\"utf-8\">"
        "<title>Variant assessment</title>"
        "<style>body{font-family:sans-serif;margin:2rem;} table{border-collapse:collapse;} "
        "td,th{border:1px solid #ccc;padding:0.3rem 0.5rem;}</style></head><body>"
        f"<pre>{body}</pre>{table}</body></html>\n"
    )


def write_sample_outputs(
    candidates: pd.DataFrame,
    liability: pd.DataFrame,
    qc: dict,
    outdir: Path,
    sample_id: str,
    summary: dict | None = None,
) -> dict[str, Path]:
    """Write Assess outputs: ``candidates.tsv``, ``report.html`` / ``report.md``, ``qc_summary.json``.

    Liability classification still feeds the report narrative. Follow-up AR
    alleles are rows in ``candidates.tsv`` with ``followup_reason`` set.
    ``sample_summary.tsv`` is the fourth ``candidates`` array member (TSV only;
    not copied into the HTML report or ``qc_summary.json``).
    """
    from .sample_summary import SUMMARY_TSV_COLUMNS, summary_to_row

    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    candidates = order_candidate_columns(candidates if candidates is not None else pd.DataFrame())
    paths = {
        "candidates_tsv": outdir / "candidates.tsv",
        "report_md": outdir / "report.md",
        "report_html": outdir / "report.html",
        "qc_json": outdir / "qc_summary.json",
        "sample_summary_tsv": outdir / "sample_summary.tsv",
    }
    candidates.to_csv(paths["candidates_tsv"], sep="\t", index=False)
    markdown = render_markdown(sample_id, candidates, liability, qc)
    paths["report_md"].write_text(markdown, encoding="utf-8")
    paths["report_html"].write_text(render_html(markdown, candidates), encoding="utf-8")
    paths["qc_json"].write_text(json.dumps(_jsonable(qc), indent=2) + "\n", encoding="utf-8")
    row = summary_to_row(summary or {"sample_id": sample_id})
    pd.DataFrame([row], columns=SUMMARY_TSV_COLUMNS).to_csv(paths["sample_summary_tsv"], sep="\t", index=False)
    return paths
