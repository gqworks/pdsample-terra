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
    followup,
    gba1_pd,
    inheritance,
    io,
    liability,
    mtdna,
    normalize,
    phenotype,
    report,
    scoring,
    splicing,
    sv,
)
from .avi import merge_avi
from .sample_summary import build_sample_summary
from .splice_heuristics import IndexedFasta, annotate_splice_heuristics, prepare_exons
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


def _integrate_special_calls(df: pd.DataFrame, cfg: dict | None = None) -> pd.DataFrame:
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
        if "followup_reason" not in df.columns:
            df["followup_reason"] = ""
        locus = (
            df["str_locus"].astype("string").fillna("").astype(str)
            if "str_locus" in df.columns
            else pd.Series("", index=df.index, dtype=str)
        )
        gene = (
            df["gene"].astype("string").fillna("").astype(str)
            if "gene" in df.columns
            else pd.Series("", index=df.index, dtype=str)
        )
        locus = locus.mask(locus.str.lower().isin(["", "nan", "<na>", "none"]), "")
        gene = gene.mask(gene.str.lower().isin(["", "nan", "<na>", "none"]), "")
        keys = locus.where(locus.str.strip().ne(""), gene)
        ar = expanded & keys.map(lambda loc: callsets.str_moi_for(str(loc or ""), cfg or {}) == "AR")
        n_hit = pd.to_numeric(df.get("str_n_expanded_alleles"), errors="coerce")
        motif_ok = df["str_motif_pathogenic"].fillna(False).astype(bool) if "str_motif_pathogenic" in df.columns else pd.Series(False, index=df.index)
        incomplete = ar & ~((n_hit.fillna(0) >= 2) & motif_ok)
        for idx in df.index[incomplete]:
            n_val = n_hit.at[idx] if idx in n_hit.index else None
            if pd.notna(n_val) and int(n_val) >= 2:
                extra = "both alleles expanded, pathogenic motif not confirmed"
            else:
                extra = "single expanded allele in recessive STR, motif not confirmed"
            df.at[idx, "followup_reason"] = compound._append_reason(df.at[idx, "followup_reason"], extra)
    max_af = None
    if cfg is not None:
        max_af = (cfg.get("filters") or {}).get("max_gnomad_af")
    return mtdna.promote_mtdna_status(df, max_af)


def _load_panel(path: Path, cfg: dict) -> pd.DataFrame:
    raw = io.read_table(path)
    return normalize.normalize_panel(raw, cfg["columns"]["panel"], cfg.get("gene_aliases"))


def _phenotype_payload(sample: dict) -> dict | None:
    payload = sample.get("phenotype")
    return payload if isinstance(payload, dict) else None


def _phenotype_frame(sample_id: str, sex: str, payload: dict | None) -> pd.DataFrame:
    row = {"sample_id": sample_id, "sex": sex}
    row.update(phenotype.phenotype_for_frame(payload))
    return pd.DataFrame([row])


def _annotate_splice_heuristics(
    scored: pd.DataFrame,
    notes: dict,
    fasta_path: Path | None,
    fai_path: Path | None,
    gtf_path: Path | None,
    panel: pd.DataFrame,
    cfg: dict,
) -> pd.DataFrame:
    """Label Figure 3/4 splice windows. A missing FASTA or GTF leaves the columns blank."""
    if notes.get("splice_heuristics") != "called":
        return annotate_splice_heuristics(scored, None, None)
    fasta = None
    try:
        fasta = IndexedFasta(fasta_path, fai_path)
        genes = panel["gene"] if "gene" in panel.columns else None
        exons = prepare_exons(gtf_path, genes, cfg.get("gene_aliases"))
        return annotate_splice_heuristics(scored, exons, fasta)
    except (OSError, ValueError, UnicodeError):
        notes["splice_heuristics"] = "not_assessed"
        return annotate_splice_heuristics(scored, None, None)
    finally:
        if fasta is not None:
            fasta.close()


def _empty_candidates() -> pd.DataFrame:
    return pd.DataFrame(columns=["rank", "sample_id", "gene", "inheritance_status"])


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
    sv_cfg = cfg.get("sv") or {}
    table, sv_assessment = sv.ingest_sv(
        sv_path,
        bed,
        exons,
        panel,
        aliases=cfg.get("gene_aliases"),
        exclude_filters=sv_cfg.get("exclude_filters"),
        dosage_rules=sv_cfg.get("dosage_rules"),
        snv_samples=[sample_id],
        depth_artifact=sv_cfg,
    )
    frames.append(_relabel(table, sample_id))
    # Workflow may set sv_status to "prior" (reuse this.sv_vcf), "joint"
    # (subset workspace.sv_joint_vcf), or "called" (GATKSVPipelineSingleSample).
    sv_status = str(sample.get("sv_status") or "called")
    notes["gatk_sv"] = sv_status

    gba1_browser = _exists(sample.get("gba1_pd_browser")) or gba1_pd.default_browser_path(cfg)
    callset_table, callset_assessment = callsets.ingest_sample_callsets(
        sample_id,
        expansionhunter=_exists(sample.get("expansionhunter_vcf")),
        gauchian=_exists(sample.get("gauchian_json")),
        mtdna=_exists(sample.get("mtdna_vcf")),
        cfg=cfg,
        gba1_pd_browser=gba1_browser,
    )
    frames.append(callset_table)
    for module, label in (
        ("expansionhunter", "expansionhunter"),
        ("gauchian", "gauchian"),
        ("mtdna", "mtdna"),
    ):
        rows = callset_assessment[callset_assessment["module"] == module] if not callset_assessment.empty else callset_assessment
        notes[label] = "not_assessed" if rows.empty or rows.iloc[0]["state"] == "not_assessed" else "called"
    notes["mtdna_annotation"] = str(sample.get("mtdna_annotation") or "not_assessed")

    variants = pd.concat([frame for frame in frames if frame is not None and not frame.empty], ignore_index=True) if any(not frame.empty for frame in frames) else pd.DataFrame()
    if not variants.empty and "gene" in variants.columns:
        variants["gene"] = normalize.apply_gene_aliases(variants["gene"], cfg.get("gene_aliases"))
    variants = splicing.ensure_hit_source(variants, default="small_variant") if not variants.empty else variants
    if not variants.empty:
        variants = gba1_pd.annotate_gauchian(variants, gba1_pd.load_browser(gba1_browser))
        variants = callsets.drop_redundant_gauchian_snvs(variants)
    unfiltered = variants.copy() if not variants.empty else variants

    introme_path = _exists(sample.get("introme_tsv"))
    introme_variants = pd.DataFrame()
    if introme_path is None:
        notes["introme"] = "not_assessed"
    else:
        raw_introme = io.read_table(introme_path)
        if not raw_introme.empty:
            # Scores attach to callset alleles. Gene-scoped injection of uncalled
            # Introme hits runs later, only in genes tagged with followup_reason.
            introme_variants, _carriers = normalize.normalize_introme(
                raw_introme, cfg["columns"]["introme"], cfg.get("gene_aliases")
            )
        notes["introme"] = str(sample.get("introme_status") or "called")
    notes["vep"] = "called"

    avi_path = _exists(sample.get("avi_tsv"))
    avi_table = pd.DataFrame()
    if avi_path is not None:
        avi_table = io.read_table(avi_path)
    notes["avi"] = "called" if not avi_table.empty else "not_assessed"
    fasta_path = _exists(sample.get("ref_fasta"))
    fai_path = _exists(sample.get("ref_fasta_index"))
    gtf_path = _exists(sample.get("gtf"))
    notes["splice_heuristics"] = "called" if fasta_path and fai_path and gtf_path else "not_assessed"

    if variants.empty:
        scored = variants
    else:
        scored = filters.apply_filters(variants, panel, cfg["filters"])
        if not scored.empty:
            scored = merge_avi(scored, avi_table)
            scored = scoring.score_variants(scored, cfg["scoring"])
            scored = sv.apply_dosage_evidence(scored, cfg["scoring"])
            scored = mtdna.apply_mtdna_evidence(scored, cfg["scoring"])
            scored = gba1_pd.annotate_and_score(scored, cfg, browser_path=gba1_browser)
            if not introme_variants.empty:
                scored = splicing.merge_introme(scored, introme_variants)
            scored = splicing.apply_introme_evidence(scored, cfg["splicing"], cfg["scoring"])
            scored = splicing.summarize_splice_scores(scored)
            scored = _annotate_splice_heuristics(scored, notes, fasta_path, fai_path, gtf_path, panel, cfg)
            scored = inheritance.annotate_moi(scored, panel)
            scored = inheritance.evaluate_inheritance(scored, cfg["inheritance"])
            scored = compound.refine(scored)
            gauchian_neg = notes.get("gauchian") == "called" and (
                unfiltered.empty
                or "hit_source" not in unfiltered.columns
                or not unfiltered["hit_source"].fillna("").astype(str).eq("gauchian").any()
            )
            scored = followup.rescue_second_hits(
                scored,
                unfiltered,
                introme_variants=introme_variants,
                cfg=cfg,
                sample_id=sample_id,
            )
            if not introme_variants.empty and "introme_score" in scored.columns:
                scored = splicing.apply_introme_evidence(scored, cfg["splicing"], cfg["scoring"])
            scored = followup.annotate_empty_searches(scored, gauchian_catalog_negative=gauchian_neg)
            scored = _integrate_special_calls(scored, cfg)
            pheno = _phenotype_frame(sample_id, sex, _phenotype_payload(sample))
            scored = phenotype.phenotype_concordance(scored, pheno, cfg.get("gene_phenotypes") or {}, cfg["phenotype"])
            scored = biomarkers.apply_biomarker_priors(
                scored,
                cfg.get("biomarkers"),
                early_onset_age=(cfg.get("phenotype") or {}).get("early_onset_age"),
            )

    liability_cfg = dict(cfg.get("liability") or {})
    for key in ("str_ar_loci", "str_xl_loci", "str_loci", "str_pathogenic_motifs"):
        if key in cfg:
            liability_cfg[key] = cfg[key]
    pheno_frame = _phenotype_frame(sample_id, sex, _phenotype_payload(sample))
    commit = str((sample.get("qc") or {}).get("pdvar_git_commit") or sample.get("pdvar_git_commit") or "")
    avi_requested = bool(sample.get("avi_requested"))

    def _maybe_stubs(candidates: pd.DataFrame, classification: str) -> pd.DataFrame:
        if classification != "unexplained" or (candidates is not None and not candidates.empty):
            return candidates
        stubs = followup.module_stub_rows(
            sample_id,
            notes,
            sv_status=sv_status,
            avi_requested=avi_requested,
            pdvar_git_commit=commit,
        )
        return stubs if not stubs.empty else candidates

    if scored.empty or "inheritance_status" not in scored.columns:
        candidates = _empty_candidates()
        classes = liability.classify_cohort(
            pd.DataFrame(columns=["sample_id", "inheritance_status"]),
            [sample_id],
            pheno_frame,
            liability_cfg,
        )
        classification = str(classes.iloc[0]["classification"]) if not classes.empty else "unexplained"
        candidates = _maybe_stubs(candidates, classification)
    else:
        classes = liability.classify_cohort(
            scored,
            [sample_id],
            pheno_frame,
            liability_cfg,
        )
        classification = str(classes.iloc[0]["classification"]) if not classes.empty else "unexplained"
        onset = None
        if not pheno_frame.empty and "age_onset" in pheno_frame.columns:
            onset = pheno_frame.iloc[0].get("age_onset")
        scored = followup.early_onset_keep(
            scored,
            age_onset=onset,
            classification=classification,
            cfg=cfg,
        )
        candidates = report.patient_candidates(scored, cfg["report"]["max_candidates_per_patient"])
        candidates = _maybe_stubs(candidates, classification)

    qc = dict(sample.get("qc") or {})
    qc["sample_id"] = sample_id
    qc["modules"] = notes
    qc["sv_status"] = sv_status
    patient_phenotype = phenotype.clean_phenotype(phenotype.prepare_phenotype(_phenotype_payload(sample)))
    qc["phenotype"] = {"provided": bool(patient_phenotype), **patient_phenotype}
    if not callset_assessment.empty:
        qc["callsets"] = callset_assessment.to_dict(orient="records")
    if not sv_assessment.empty:
        qc["sv"] = sv_assessment.to_dict(orient="records")
    summary = build_sample_summary(sample, cfg)
    qc["sample_summary"] = summary
    notes["haplogrep"] = str(summary.get("haplogroup_status") or "not_assessed")
    notes["panel_prs"] = str(summary.get("prs_status") or "not_assessed")
    notes["apoe"] = str(summary.get("apoe_status") or "not_assessed")
    notes["mt_2158"] = str(summary.get("mt_2158_status") or "not_assessed")
    notes["mean_mlc"] = str(summary.get("mlc_status") or "not_assessed")
    qc["modules"] = notes
    return report.write_sample_outputs(candidates, classes, qc, outdir, sample_id, summary=summary)


def load_sample_json(path: str | Path) -> dict:
    return json.loads(Path(path).read_text(encoding="utf-8"))
