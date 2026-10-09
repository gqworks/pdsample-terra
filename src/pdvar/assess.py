"""Per-sample assessment from the callsets PdSample writes.

The sample JSON is an internal handoff between the WDL Assess task and this
function. It is not a cohort manifest.
"""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

from . import (
    biomarkers,
    callsets,
    compound,
    filters,
    inheritance,
    io,
    isoform,
    liability,
    normalize,
    phenotype,
    report,
    scoring,
    splicing,
    sv,
)
from .incidental import apply_incidental_tags
from .vep import parse_vep

def _path(value: object) -> Path | None:
    if value is None:
        return None
    text = str(value).strip()
    if text in {"", "null", "None"}:
        return None
    return Path(text)


def _exists(value: object) -> Path | None:
    path = _path(value)
    if path is None or not path.exists():
        return None
    return path


def _relabel(table: pd.DataFrame, sample_id: str) -> pd.DataFrame:
    if table.empty or "sample_id" not in table.columns:
        return table
    ids = set(table["sample_id"].astype(str))
    if sample_id in ids:
        return table[table["sample_id"].astype(str) == sample_id].copy()
    if len(ids) == 1:
        table = table.copy()
        table["sample_id"] = sample_id
        if "family_id" in table.columns:
            table["family_id"] = sample_id
        return table
    return table.iloc[0:0].copy()


def _male_x(df: pd.DataFrame, sex: str) -> pd.DataFrame:
    """A homozygous X genotype in a male is hemizygous for XL logic."""
    if sex != "male" or df.empty or "chrom" not in df.columns:
        return df
    df = df.copy()
    mask = df["chrom"].astype(str).eq("X") & df["zygosity"].astype(str).eq("hom")
    if "hit_source" in df.columns:
        mask &= df["hit_source"].astype(str).eq("small_variant")
    df.loc[mask, "zygosity"] = "hemi"
    return df


def _integrate_special_calls(df: pd.DataFrame) -> pd.DataFrame:
    """Keep mtDNA carriers and cited STR expansions visible after MOI logic."""
    if df.empty or "hit_source" not in df.columns:
        return df
    df = df.copy()
    status = df["inheritance_status"].astype(str) if "inheritance_status" in df.columns else pd.Series("", index=df.index)
    mt = df["hit_source"].eq("mtdna") & status.isin(["not_qualifying", "nan", ""])
    df.loc[mt, "inheritance_status"] = "mtdna_variant"
    if "str_expanded" in df.columns:
        expanded = df["hit_source"].eq("expansionhunter") & df["str_expanded"].fillna(False).astype(bool)
        df.loc[expanded, "inheritance_status"] = "str_expanded"
    return df


def _load_panel(path: Path, cfg: dict) -> pd.DataFrame:
    raw = io.read_table(path)
    return normalize.normalize_panel(raw, cfg["columns"]["panel"], cfg.get("gene_aliases"))


def _phenotype_frame(sample_id: str, sex: str, payload: dict | None) -> pd.DataFrame:
    row = {"sample_id": sample_id, "sex": sex}
    if payload:
        row.update({key: value for key, value in payload.items() if key != "sample_id"})
    return pd.DataFrame([row])


def _empty_candidates() -> pd.DataFrame:
    return pd.DataFrame(columns=["rank", "sample_id", "gene", "inheritance_status", "incidental"])


def assess_sample(
    sample: dict,
    cfg: dict,
    *,
    panel_csv: str | Path,
    outdir: str | Path,
    gene_bed: str | Path | None = None,
    exon_bed: str | Path | None = None,
) -> dict[str, Path]:
    """Rank one sample and write the Terra output files."""
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    sample_id = str(sample["sample_id"])
    sex = str(sample.get("sex") or "unknown")
    panel = _load_panel(Path(panel_csv), cfg)
    frames: list[pd.DataFrame] = []
    notes: dict[str, str] = {}

    vep_path = _exists(sample.get("vep_vcf"))
    if vep_path is None:
        raise FileNotFoundError("sample JSON vep_vcf is missing")
    small = _relabel(parse_vep(vep_path, sample_id), sample_id)
    if not small.empty:
        small["gene"] = normalize.apply_gene_aliases(small["gene"], cfg.get("gene_aliases"))
        small = _male_x(small, sex)
    frames.append(small)
    notes["small_variants"] = "called"

    sv_path = _exists(sample.get("sv_vcf"))
    if sv_path is None:
        raise FileNotFoundError("sample JSON sv_vcf is missing")
    bed = _exists(gene_bed) or _exists(cfg.get("paths", {}).get("gene_bed"))
    exons = _exists(exon_bed) or _exists(cfg.get("paths", {}).get("exon_bed"))
    if bed is None or exons is None:
        raise FileNotFoundError("gene_bed and exon_bed are required to assess the SV VCF")
    table, sv_assessment = sv.ingest_sv(
        sv_path,
        bed,
        exons,
        panel,
        aliases=cfg.get("gene_aliases"),
        exclude_filters=(cfg.get("sv") or {}).get("exclude_filters"),
        dosage_rules=(cfg.get("sv") or {}).get("dosage_rules"),
        snv_samples=[sample_id],
    )
    frames.append(_relabel(table, sample_id))
    notes["gatk_sv"] = "called"
    sv_status = "called"

    callset_table, callset_assessment = callsets.ingest_sample_callsets(
        sample_id,
        expansionhunter=_exists(sample.get("expansionhunter_vcf")),
        gauchian=_exists(sample.get("gauchian_json")),
        mtdna=_exists(sample.get("mtdna_vcf")),
        cfg=cfg,
    )
    frames.append(callset_table)
    for module, label in (
        ("expansionhunter", "expansionhunter"),
        ("gauchian", "gauchian"),
        ("mtdna", "mtdna"),
    ):
        rows = callset_assessment[callset_assessment["module"] == module] if not callset_assessment.empty else callset_assessment
        notes[label] = "not_assessed" if rows.empty or rows.iloc[0]["state"] == "not_assessed" else "called"

    variants = pd.concat([frame for frame in frames if frame is not None and not frame.empty], ignore_index=True) if any(not frame.empty for frame in frames) else pd.DataFrame()
    if not variants.empty and "gene" in variants.columns:
        variants["gene"] = normalize.apply_gene_aliases(variants["gene"], cfg.get("gene_aliases"))
    variants = splicing.ensure_hit_source(variants, default="small_variant") if not variants.empty else variants

    introme_path = _exists(sample.get("introme_tsv"))
    introme_variants = pd.DataFrame()
    introme_carriers = pd.DataFrame()
    if introme_path is None:
        notes["introme"] = "not_assessed"
    else:
        raw_introme = io.read_table(introme_path)
        if not raw_introme.empty:
            introme_variants, introme_carriers = normalize.normalize_introme(
                raw_introme, cfg["columns"]["introme"], cfg.get("gene_aliases")
            )
            introme_carriers = _relabel(introme_carriers, sample_id)
        notes["introme"] = "called"
    notes["vep"] = "called"

    follow = pd.DataFrame()
    if variants.empty:
        scored = variants
    else:
        scored = filters.apply_filters(variants, panel, cfg["filters"])
        if scored.empty:
            follow = pd.DataFrame()
        else:
            scored = scoring.score_variants(scored, cfg["scoring"])
            scored = sv.apply_dosage_evidence(scored, cfg["scoring"])
            if not introme_variants.empty:
                scored = splicing.merge_introme(scored, introme_variants)
            scored = splicing.apply_introme_evidence(scored, cfg["splicing"], cfg["scoring"])
            vault = _exists(sample.get("splicevault"))
            if vault is not None:
                scored = splicing.merge_splicevault(scored, io.read_table(vault))
            hits = splicing.find_second_hits(scored, introme_carriers, introme_variants, panel, cfg)
            if not hits.empty:
                hits["sample_id"] = sample_id
                scored = pd.concat([scored, hits], ignore_index=True)
            scored = splicing.summarize_splice_scores(scored)
            scored = inheritance.annotate_moi(scored, panel)
            scored = inheritance.evaluate_inheritance(scored, cfg["inheritance"])
            phased = _exists(sample.get("phased_vcf"))
            if phased is not None:
                try:
                    scored = compound.attach_phase_from_vcf(scored, phased)
                except (RuntimeError, FileNotFoundError):
                    pass
            scored, follow = compound.refine(scored, cfg.get("compound"))
            scored = _integrate_special_calls(scored)
            scored = apply_incidental_tags(scored)
            pheno = _phenotype_frame(sample_id, sex, sample.get("phenotype") if isinstance(sample.get("phenotype"), dict) else None)
            scored = phenotype.phenotype_concordance(scored, pheno, cfg.get("gene_phenotypes") or {}, cfg["phenotype"])
            scored = biomarkers.apply_biomarker_priors(
                scored,
                cfg.get("biomarkers"),
                early_onset_age=(cfg.get("phenotype") or {}).get("early_onset_age"),
            )
            expression = None
            expression_path = _exists(sample.get("expression"))
            if expression_path is not None:
                expression = io.read_table(expression_path)
            scored = isoform.weight_lof(scored, expression, cfg.get("expression") or {})

    if scored.empty or "inheritance_status" not in scored.columns:
        candidates = _empty_candidates()
        classes = liability.classify_cohort(
            pd.DataFrame(columns=["sample_id", "inheritance_status"]),
            [sample_id],
            _phenotype_frame(sample_id, sex, sample.get("phenotype") if isinstance(sample.get("phenotype"), dict) else None),
            None,
            cfg.get("liability") or {},
        )
    else:
        candidates = report.patient_candidates(scored, cfg["report"]["max_candidates_per_patient"])
        prs = None
        weights_path = _exists(sample.get("prs_weights"))
        if weights_path is not None:
            weights = liability.load_prs_weights(weights_path)
            if weights is not None:
                prs = liability.score_prs(vep_path, weights, None)
                if prs is not None and not prs.empty and sample_id not in set(prs["sample_id"].astype(str)) and len(prs) == 1:
                    prs = prs.copy()
                    prs["sample_id"] = sample_id
        classes = liability.classify_cohort(
            scored,
            [sample_id],
            _phenotype_frame(sample_id, sex, sample.get("phenotype") if isinstance(sample.get("phenotype"), dict) else None),
            prs,
            cfg.get("liability") or {},
        )

    qc = dict(sample.get("qc") or {})
    qc["sample_id"] = sample_id
    qc["modules"] = notes
    qc["sv_status"] = sv_status
    if not callset_assessment.empty:
        qc["callsets"] = callset_assessment.to_dict(orient="records")
    if not sv_assessment.empty:
        qc["sv"] = sv_assessment.to_dict(orient="records")
    return report.write_sample_outputs(candidates, classes, qc, follow, outdir, sample_id)


def load_sample_json(path: str | Path) -> dict:
    return json.loads(Path(path).read_text(encoding="utf-8"))
