"""Per-sample explained-liability report.

Classification is qualitative:

- ``solved``: a mode-of-inheritance-complete finding in a gene that is not on
  the reduced-penetrance list (default GBA1 and LRRK2, from the panel notes)
- ``partially_explained``: a risk allele, a reduced-penetrance allele, or a
  single heterozygous hit in a recessive gene
- ``unexplained``: none of those

``explained_liability`` is numeric only for penetrance values that are present
in config with a citation. It is ``1 - product(1 - p)`` over those values.
PRS, sex, and age are reported beside it. They are folded in only when the
user supplies a cited coefficient. No penetrance, PRS weight, or allele effect
is shipped in the default config.
"""

from __future__ import annotations

import math
from pathlib import Path

import pandas as pd

from .vcf_stream import chrom_bare, iter_vcf

SOLVED_STATUSES = {
    "AR_hom",
    "AR_comphet_candidate",
    "AR_comphet_phased",
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
}
BIALLELIC = {"AR_hom", "AR_comphet_candidate", "AR_comphet_phased", "AR_hemizygous", "risk_biallelic"}

PRS_COLUMNS = ["chrom", "pos", "effect_allele", "other_allele", "weight", "rsid", "citation"]


def load_prs_weights(path: str | Path | None) -> pd.DataFrame | None:
    if not path:
        return None
    path = Path(path)
    if not path.exists():
        return None
    table = pd.read_csv(path, sep="\t", comment="#")
    if table.empty:
        return None
    missing = [col for col in ("chrom", "pos", "effect_allele", "weight", "citation") if col not in table.columns]
    if missing:
        raise ValueError(
            f"PRS weights file is missing {missing}. Expected columns: {', '.join(PRS_COLUMNS)}. "
            "See resources/prs/weights_format.tsv."
        )
    if table["citation"].isna().any() or table["citation"].astype(str).str.strip().isin(["", "nan"]).any():
        raise ValueError("Every PRS weight row needs a citation. This pipeline does not ship unsourced weights.")
    if table["weight"].isna().any():
        raise ValueError("Every PRS weight row needs a numeric weight.")
    return table


def _dosage(gt: str, effect: str, other: str | None, ref: str, alt: str) -> int | None:
    if not gt or "." in gt:
        return None
    alleles = gt.replace("|", "/").split("/")
    try:
        coded = []
        for allele in alleles:
            base = ref if allele == "0" else alt.split(",")[int(allele) - 1]
            coded.append(base)
    except (ValueError, IndexError):
        return None
    return sum(base == effect for base in coded)


def score_prs(vcf: str | Path, weights: pd.DataFrame, samples: list[str] | None = None) -> pd.DataFrame:
    """Region-query ``weights`` positions and return a raw dosage x weight sum per sample."""
    import tempfile

    totals = {sample: 0.0 for sample in samples} if samples else {}
    used = {sample: 0 for sample in totals}
    seen: set[tuple] = set()
    with tempfile.TemporaryDirectory() as tmp:
        bed = Path(tmp) / "prs.bed"
        with bed.open("w", encoding="utf-8") as fh:
            for row in weights.itertuples(index=False):
                pos = int(row.pos)
                fh.write(f"chr{chrom_bare(str(row.chrom))}\t{pos - 1}\t{pos}\n")
        for record in iter_vcf(vcf, regions=bed):
            for row in weights.itertuples(index=False):
                if chrom_bare(str(row.chrom)) != chrom_bare(record["chrom"]) or int(row.pos) != int(record["pos"]):
                    continue
                if str(row.effect_allele) not in {record["ref"], *str(record["alt"]).split(",")}:
                    continue
                seen.add((chrom_bare(str(row.chrom)), int(row.pos), str(row.effect_allele)))
                for sample, fmt in record["samples"].items():
                    if samples is not None and sample not in totals:
                        continue
                    dosage = _dosage(fmt.get("GT", ""), str(row.effect_allele), getattr(row, "other_allele", None), record["ref"], record["alt"])
                    if dosage is None:
                        continue
                    totals[sample] = totals.get(sample, 0.0) + dosage * float(row.weight)
                    used[sample] = used.get(sample, 0) + 1
    n_weights = len(weights)
    rows = []
    for sample, score in totals.items():
        rows.append(
            {
                "sample_id": sample,
                "prs": score,
                "prs_sites_used": used.get(sample, 0),
                "prs_sites_missing": n_weights - used.get(sample, 0),
                "prs_citation": ";".join(sorted(set(weights["citation"].astype(str)))),
            }
        )
    # Samples with no genotypes still appear when ``samples`` was provided.
    # Sites that were not in the VCF count as missing for everyone.
    if not rows and samples:
        rows = [
            {"sample_id": sample, "prs": 0.0, "prs_sites_used": 0, "prs_sites_missing": n_weights, "prs_citation": ""}
            for sample in samples
        ]
    return pd.DataFrame(rows)


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
    prs: pd.DataFrame | None = None,
    cfg: dict | None = None,
) -> pd.DataFrame:
    cfg = cfg or {}
    reduced = {gene.upper() for gene in cfg.get("reduced_penetrance_genes") or []}
    penetrance_table = cfg.get("penetrance") or {}
    prs_map = {}
    if prs is not None and not prs.empty:
        prs_map = prs.set_index(prs["sample_id"].astype(str)).to_dict(orient="index")
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
        if not grp.empty and "inheritance_status" in grp.columns:
            for rec in grp.to_dict(orient="records"):
                status = str(rec.get("inheritance_status") or "")
                gene = str(rec.get("gene") or "").upper()
                if status in {"", "not_qualifying", "nan", "incidental"}:
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
                if use_solved:
                    solved.append(gene)
                elif status in PARTIAL_STATUSES or gene in reduced:
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
        pr = prs_map.get(str(sample), {})
        prs_beta = cfg.get("prs_beta")
        prs_citation = cfg.get("prs_beta_citation")
        prs_value = pr.get("prs")
        prs_contrib = None
        if prs_value is not None and prs_beta is not None and prs_citation:
            prs_contrib = float(prs_beta) * float(prs_value)
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
                "prs": prs_value if prs_value is not None else pd.NA,
                "prs_sites_used": pr.get("prs_sites_used", pd.NA),
                "prs_sites_missing": pr.get("prs_sites_missing", pd.NA),
                "prs_liability_contribution": pd.NA if prs_contrib is None else prs_contrib,
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
