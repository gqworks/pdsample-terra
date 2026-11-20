"""Per-sample explained-liability report.

Classification is qualitative:

- ``solved``: a mode-of-inheritance-complete finding in a gene that is not on
  the reduced-penetrance list (default GBA1 and LRRK2, from the panel notes)
- ``partially_explained``: a risk allele, a reduced-penetrance allele, or a
  single heterozygous hit in a recessive gene
- ``unexplained``: none of those

``explained_liability`` is numeric only for penetrance values that are present
in config with a citation. It is ``1 - product(1 - p)`` over those values.
Sex and age are reported beside it. They are folded in only when the user
supplies a cited coefficient. No penetrance or allele effect is shipped in
the default config.
"""

from __future__ import annotations

import math

import pandas as pd

from .callsets import DEFAULT_STR_AR, str_moi_for
from .inheritance import weak_intronic_mask

SOLVED_STATUSES = {
    "AR_hom",
    "AR_comphet_candidate",
    "AR_hemizygous",
    "AR_second_hit_splice",
    "AD_candidate",
    "XL_hemi",
}
PARTIAL_STATUSES = {
    "AR_single_het",
    "risk_carrier",
    "risk_biallelic",
    "MOI_unknown_candidate",
    "review_conflict",
    "mtdna_pathogenic_candidate",
}


def _truthy(value: object) -> bool:
    if value is None or (isinstance(value, float) and pd.isna(value)):
        return False
    if isinstance(value, bool):
        return value
    text = str(value).strip().lower()
    return text in {"1", "true", "yes", "t"}


def _str_n_expanded(rec: dict) -> int | None:
    raw = rec.get("str_n_expanded_alleles")
    if raw is None or (isinstance(raw, float) and pd.isna(raw)):
        return None
    try:
        return int(raw)
    except (TypeError, ValueError):
        return None


def _classify_str_expanded(rec: dict, cfg: dict) -> str:
    """``solved`` / ``partial`` for a cited STR expansion.

    AD and XL expansions solve. AR expansions (RFC1, CSTB, FXN) solve only when
    both alleles meet the floor and the pathogenic motif is confirmed. Until the
    motif is parsed, keep AR expansions partial.
    """
    gene = str(rec.get("gene") or rec.get("str_locus") or "").upper()
    locus = str(rec.get("str_locus") or gene)
    moi = str(rec.get("str_moi") or rec.get("moi") or "").strip().upper()
    if moi not in {"AD", "AR", "XL"}:
        moi = str_moi_for(locus, cfg)
    n_hit = _str_n_expanded(rec)
    motif_ok = rec.get("str_motif_pathogenic")
    ar_loci = {str(item).upper() for item in (cfg.get("str_ar_loci") or DEFAULT_STR_AR)}
    is_ar = moi == "AR" or gene in ar_loci or locus in ar_loci
    if is_ar:
        biallelic = n_hit is not None and n_hit >= 2
        if biallelic and _truthy(motif_ok):
            return "solved"
        return "partial"
    return "solved"
BIALLELIC = {"AR_hom", "AR_comphet_candidate", "AR_hemizygous", "risk_biallelic"}


def _penetrance(gene: str, status: str, table: dict) -> tuple[float | None, str | None]:
    rec = table.get(gene) or table.get(gene.upper()) or {}
    if not isinstance(rec, dict):
        return None, None
    key = "biallelic" if status in BIALLELIC else "heterozygous"
    value = rec.get(key)
    citation = rec.get("citation")
    if value is None or not citation:
        return None, citation if citation else None
    return float(value), str(citation)


def _combine(penetrances: list[float]) -> float | None:
    if not penetrances:
        return None
    remaining = 1.0
    for value in penetrances:
        remaining *= 1.0 - value
    return 1.0 - remaining


def classify_cohort(
    variants: pd.DataFrame,
    samples: list[str],
    pheno: pd.DataFrame | None = None,
    cfg: dict | None = None,
) -> pd.DataFrame:
    cfg = cfg or {}
    reduced = {gene.upper() for gene in cfg.get("reduced_penetrance_genes") or []}
    penetrance_table = cfg.get("penetrance") or {}
    pheno_map = {}
    if pheno is not None and not pheno.empty and "sample_id" in pheno.columns:
        pheno_map = pheno.set_index(pheno["sample_id"].astype(str)).to_dict(orient="index")

    rows = []
    grouped = variants.groupby(variants["sample_id"].astype(str)) if not variants.empty and "sample_id" in variants.columns else []
    by_sample = {str(sample): grp for sample, grp in grouped} if not isinstance(grouped, list) else {}
    for sample in samples:
        grp = by_sample.get(str(sample), variants.iloc[0:0])
        solved, partial, components = [], [], []
        penetrances = []
        intron_ad = pd.Series(False, index=grp.index)
        if not grp.empty and "inheritance_status" in grp.columns:
            # MAPT splice_gof skips GoF demotion (stays AD_candidate / force-kept)
            # but intron-only AD still is not a sample-level solve.
            intron_ad = grp["inheritance_status"].fillna("").astype(str).eq("AD_candidate") & weak_intronic_mask(grp)
        intron_ids = set(grp.index[intron_ad]) if intron_ad.any() else set()
        if not grp.empty and "inheritance_status" in grp.columns:
            for rec, idx in zip(grp.to_dict(orient="records"), grp.index):
                status = str(rec.get("inheritance_status") or "")
                gene = str(rec.get("gene") or "").upper()
                if status in {"", "not_qualifying", "nan"}:
                    continue
                value, citation = _penetrance(gene, status, penetrance_table)
                component = {
                    "gene": gene,
                    "status": status,
                    "variant_id": rec.get("variant_id"),
                    "variant_class": rec.get("variant_class"),
                    "penetrance": value,
                    "citation": citation,
                    "age_at_onset_median": (penetrance_table.get(gene) or {}).get("age_at_onset_median")
                    if isinstance(penetrance_table.get(gene), dict)
                    else None,
                }
                components.append(component)
                use_solved = status in SOLVED_STATUSES and gene not in reduced
                if status == "AR_second_hit_splice":
                    use_solved = gene not in reduced
                if status == "AD_candidate" and idx in intron_ids:
                    # Intron-only AD (no ClinVar P/LP, no coding/canonical splice) is not a solve.
                    use_solved = False
                if status == "str_expanded":
                    bucket = _classify_str_expanded(rec, cfg)
                    if bucket == "solved" and gene not in reduced:
                        solved.append(gene)
                    else:
                        partial.append(gene)
                elif use_solved:
                    solved.append(gene)
                elif status in PARTIAL_STATUSES or gene in reduced or (
                    status == "AD_candidate" and idx in intron_ids
                ):
                    partial.append(gene)
                if value is not None:
                    penetrances.append(value)
        if solved:
            classification = "solved"
        elif partial:
            classification = "partially_explained"
        else:
            classification = "unexplained"
        contributing = set(solved + partial)
        ph = pheno_map.get(str(sample), {})
        rows.append(
            {
                "sample_id": str(sample),
                "classification": classification,
                "oligogenic": classification == "partially_explained" and len(contributing) >= 2,
                "solved_genes": ";".join(sorted(set(solved))),
                "partial_genes": ";".join(sorted(set(partial))),
                "explained_liability": _combine(penetrances),
                "explained_liability_note": (
                    "1 - product(1 - p) over cited penetrances"
                    if penetrances
                    else "null until a cited penetrance is set in config liability.penetrance"
                ),
                "sex": ph.get("sex", pd.NA),
                "age_onset": ph.get("age_onset", pd.NA),
                "sex_coefficient": cfg.get("sex_coefficient"),
                "age_coefficient": cfg.get("age_coefficient"),
                "components_json": pd.Series([components]).astype(str).iloc[0] if False else _json(components),
            }
        )
    return pd.DataFrame(rows)


def _json(components: list[dict]) -> str:
    import json

    def convert(value):
        if value is None or (isinstance(value, float) and math.isnan(value)):
            return None
        if isinstance(value, float):
            return value
        return value if not hasattr(value, "item") else value.item()

    clean = [{key: convert(val) for key, val in comp.items()} for comp in components]
    return json.dumps(clean, default=str)
