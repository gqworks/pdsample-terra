"""Per-sample candidate lists and a short HTML/Markdown report."""

from __future__ import annotations

import html
import json
from pathlib import Path

import pandas as pd

STATUS_PRIORITY = {
    "AR_hom": 0,
    "AR_comphet_candidate": 0,
    "AR_comphet_phased": 0,
    "AR_hemizygous": 0,
    "AR_second_hit_splice": 0,
    "AD_candidate": 1,
    "XL_hemi": 1,
    "risk_biallelic": 1,
    "risk_carrier": 2,
    "AR_single_het": 3,
    "review_conflict": 3,
    "mtdna_variant": 2,
    "str_expanded": 1,
    "MOI_unknown_candidate": 4,
    "incidental": 5,
    "not_qualifying": 9,
}

CURATION_COLUMNS = [
    "sample_id", "family_id", "gene", "moi", "panel_confidence", "risk_allele_exception", "variant_id", "hgvsc", "hgvsp", "consequence",
    "zygosity", "tier", "evidence_score", "evidence_reasons", "inheritance_status", "hit_source",
    "introme_class", "introme_score", "introme_cohort_af", "introme_n_carriers", "introme_gene_region",
    "clinvar_class", "gnomad_af", "gnomad_popmax_af",
    "cadd", "revel", "spliceai", "age_onset", "onset_category", "sex", "ethnicity",
    "onset_match", "matched_features", "phenotype_score", "symptoms",
    "variant_class", "dosage_mechanism", "phase_relation", "phase_source",
    "lof_weight", "weighted_evidence_score", "svtype", "biomarker_notes",
    "str_locus", "str_expanded", "heteroplasmy", "incidental", "disposition",
]


def rank_candidates(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    df["status_priority"] = df["inheritance_status"].map(STATUS_PRIORITY).fillna(9)
    if "weighted_evidence_score" in df.columns:
        df["_rank_score"] = df["weighted_evidence_score"].where(df["weighted_evidence_score"].notna(), df["evidence_score"])
    else:
        df["_rank_score"] = df["evidence_score"]
    sort_cols = ["sample_id", "status_priority", "tier", "_rank_score"]
    asc = [True, True, True, False]
    if "phenotype_score" in df.columns:
        sort_cols.append("phenotype_score")
        asc.append(False)
    df = df.sort_values(sort_cols, ascending=asc, kind="mergesort")
    df["rank"] = df.groupby("sample_id").cumcount() + 1
    return df.drop(columns="_rank_score")


def patient_candidates(df: pd.DataFrame, max_per_patient: int) -> pd.DataFrame:
    ranked = rank_candidates(df)
    ranked = ranked[ranked["inheritance_status"] != "not_qualifying"]
    return ranked[ranked["rank"] <= max_per_patient]


def gene_summary(df: pd.DataFrame) -> pd.DataFrame:
    """Per gene: number of carriers per inheritance status."""
    q = df[df["inheritance_status"] != "not_qualifying"]
    if q.empty:
        return pd.DataFrame(columns=["gene"])
    return (
        q.groupby(["gene", "inheritance_status"])["sample_id"].nunique()
        .unstack(fill_value=0)
        .assign(total_carriers=lambda t: t.sum(axis=1))
        .sort_values("total_carriers", ascending=False)
        .reset_index()
    )


def write_reports(
    df: pd.DataFrame,
    results_dir: Path,
    rcfg: dict,
    extras: dict[str, pd.DataFrame] | None = None,
) -> dict[str, Path]:
    tables = results_dir / "tables"
    reports = results_dir / "reports"
    tables.mkdir(parents=True, exist_ok=True)
    reports.mkdir(parents=True, exist_ok=True)

    cands = patient_candidates(df, rcfg["max_candidates_per_patient"])
    genes = gene_summary(df)
    cols = ["rank"] + [c for c in CURATION_COLUMNS if c in cands.columns]
    curation = cands[cols].copy()
    curation["curator_decision"] = ""
    curation["curator_notes"] = ""

    out = {
        "candidates": tables / "patient_candidates.tsv",
        "genes": tables / "gene_summary.tsv",
        "excel": reports / rcfg["excel_filename"],
    }
    cands.to_csv(out["candidates"], sep="\t", index=False)
    genes.to_csv(out["genes"], sep="\t", index=False)
    extras = extras or {}
    for name, table in extras.items():
        if table is None or table.empty:
            continue
        path = tables / f"{name}.tsv"
        table.to_csv(path, sep="\t", index=False)
        out[name] = path
    with pd.ExcelWriter(out["excel"]) as xl:
        curation.to_excel(xl, sheet_name="candidates", index=False)
        genes.to_excel(xl, sheet_name="gene_summary", index=False)
        for name, table in extras.items():
            if table is None or table.empty:
                continue
            table.to_excel(xl, sheet_name=name[:31], index=False)
    return out


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
        "## Modules",
        "",
    ]
    modules = qc.get("modules") or {}
    for name, status in modules.items():
        lines.append(f"- {name}: {status}")
    if qc.get("sv_status"):
        lines.append(f"- gatk_sv status field: {qc.get('sv_status')}")
    lines.extend(["", "## Callsets", ""])
    callsets = qc.get("callsets") or []
    if not callsets:
        lines.append("No STR, GBA1, or mtDNA assessment rows.")
    else:
        for row in callsets:
            lines.append(f"- {row.get('module')}: {row.get('callset_status')}")
    incidental = pd.DataFrame()
    if candidates is not None and not candidates.empty and "incidental" in candidates.columns:
        incidental = candidates[candidates["incidental"].fillna("").astype(str).str.len() > 0]
    if not incidental.empty:
        lines.extend(
            [
                "",
                "## Incidental",
                "",
                "TTN truncating variants are tagged incidental. They are not classified as benign and they are not treated as the Parkinson's disease diagnosis.",
                "",
            ]
        )
    lines.extend(["", "## Ranked candidates", ""])
    if candidates is None or candidates.empty:
        lines.append("No qualifying candidate variants.")
    else:
        show = [column for column in ("rank", "gene", "inheritance_status", "variant_id", "consequence", "zygosity", "tier", "evidence_score", "incidental") if column in candidates.columns]
        lines.append(candidates[show].to_csv(sep="\t", index=False).rstrip())
    lines.append("")
    return "\n".join(lines)


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
    followup: pd.DataFrame,
    outdir: Path,
    sample_id: str,
) -> dict[str, Path]:
    """Write the files the Assess task returns to the Terra row."""
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    candidates = candidates if candidates is not None else pd.DataFrame()
    paths = {
        "candidates_tsv": outdir / "candidates.tsv",
        "report_md": outdir / "report.md",
        "report_html": outdir / "report.html",
        "qc_json": outdir / "qc_summary.json",
        "liability_tsv": outdir / "liability.tsv",
        "followup_tsv": outdir / "followup.tsv",
    }
    candidates.to_csv(paths["candidates_tsv"], sep="\t", index=False)
    if liability is None or liability.empty:
        pd.DataFrame(columns=["sample_id", "classification"]).to_csv(paths["liability_tsv"], sep="\t", index=False)
    else:
        liability.to_csv(paths["liability_tsv"], sep="\t", index=False)
    if followup is None or followup.empty:
        pd.DataFrame(columns=["sample_id", "gene", "reason"]).to_csv(paths["followup_tsv"], sep="\t", index=False)
    else:
        followup.to_csv(paths["followup_tsv"], sep="\t", index=False)
    markdown = render_markdown(sample_id, candidates, liability, qc)
    paths["report_md"].write_text(markdown, encoding="utf-8")
    paths["report_html"].write_text(render_html(markdown, candidates), encoding="utf-8")
    paths["qc_json"].write_text(json.dumps(_jsonable(qc), indent=2) + "\n", encoding="utf-8")
    return paths
