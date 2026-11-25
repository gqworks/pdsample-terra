"""Panel-gene Parkinson's PRS from a cited PGS Catalog weights table.

Weights are PGS000902 / Nalls 2019 (see ``resources/prs/README.md``). This is
not the unused empty stub removed in PR 96: every row must carry a citation
and a numeric weight. ``panel_prs`` sums only ``score_in_panel_prs`` rows.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from .config import PROJECT_ROOT, resolve
from .vcf_stream import chrom_bare, iter_vcf

PRS_COLUMNS = [
    "chrom",
    "pos",
    "effect_allele",
    "other_allele",
    "weight",
    "rsid",
    "citation",
]

DEFAULT_WEIGHTS = PROJECT_ROOT / "resources" / "prs" / "pgs000902_nalls2019.tsv"
CITATION = (
    "PGS000902; Nalls MA et al. Lancet Neurol 2019;18:1091-1102. "
    "doi:10.1016/S1474-4422(19)30320-5"
)


def default_weights_path(cfg: dict | None = None) -> Path | None:
    cfg = cfg or {}
    configured = (cfg.get("inputs") or {}).get("prs_weights")
    candidates: list[Path] = []
    if configured:
        raw = Path(str(configured))
        if raw.is_absolute():
            candidates.append(raw)
        else:
            root = Path(str(cfg.get("_root") or PROJECT_ROOT))
            candidates.extend([root / raw, resolve(raw), Path.cwd() / raw])
    candidates.append(DEFAULT_WEIGHTS)
    for path in candidates:
        if path.is_file():
            return path
    return None


def load_prs_weights(path: str | Path | None, *, panel_only: bool = True) -> pd.DataFrame | None:
    """Load a cited weights TSV. Empty / missing path returns None."""
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
            f"PRS weights file is missing {missing}. Expected columns: {', '.join(PRS_COLUMNS)}."
        )
    if table["citation"].isna().any() or table["citation"].astype(str).str.strip().isin(["", "nan"]).any():
        raise ValueError("Every PRS weight row needs a citation. This pipeline does not ship unsourced weights.")
    if table["weight"].isna().any():
        raise ValueError("Every PRS weight row needs a numeric weight.")
    table = table.copy()
    table["chrom"] = table["chrom"].map(lambda value: chrom_bare(str(value)))
    table["pos"] = pd.to_numeric(table["pos"], errors="coerce").astype("Int64")
    table["effect_allele"] = table["effect_allele"].astype(str).str.upper()
    table["weight"] = pd.to_numeric(table["weight"], errors="coerce")
    if panel_only and "score_in_panel_prs" in table.columns:
        flag = table["score_in_panel_prs"].astype(str).str.strip().str.lower()
        table = table.loc[flag.isin(["1", "true", "yes", "t"])].copy()
    table = table.loc[table["pos"].notna() & table["weight"].notna()]
    return table.reset_index(drop=True)


def _dosage(gt: str, effect: str, ref: str, alt: str) -> int | None:
    """Count copies of ``effect`` in a GT string. Missing alleles return None."""
    if not gt or gt in {".", "./.", ".|.", ".|", "|."}:
        return None
    alleles = gt.replace("|", "/").split("/")
    alts = [item.strip().upper() for item in str(alt).split(",") if item.strip()]
    try:
        bases = []
        for allele in alleles:
            token = str(allele).strip()
            if token in {"", "."}:
                return None
            if token == "0":
                bases.append(str(ref).upper())
            else:
                bases.append(alts[int(token) - 1])
    except (ValueError, IndexError):
        return None
    effect = str(effect).upper()
    return sum(base == effect for base in bases)


def index_genotypes(vcf: str | Path | None) -> dict[tuple[str, int], dict]:
    """Map ``(chrom, pos)`` to the first VCF record at that site."""
    if not vcf:
        return {}
    path = Path(vcf)
    if not path.exists():
        return {}
    index: dict[tuple[str, int], dict] = {}
    for record in iter_vcf(path, allow_full_scan=True):
        key = (chrom_bare(record["chrom"]), int(record["pos"]))
        index.setdefault(key, record)
    return index


def merge_genotype_indexes(indexes: list[dict[tuple[str, int], dict]]) -> dict[tuple[str, int], dict]:
    """Earlier indexes win (summary all-sites VCF before variant-only callsets)."""
    merged: dict[tuple[str, int], dict] = {}
    for index in reversed(indexes):
        merged.update(index)
    return merged


def score_prs(
    weights: pd.DataFrame | None,
    genotypes: dict[tuple[str, int], dict],
    sample_id: str,
) -> dict:
    """Raw dosage × weight sum. Missing sites are not imputed as hom-ref."""
    empty = {
        "panel_prs": None,
        "prs_sites_used": 0,
        "prs_sites_missing": 0,
        "prs_sites_total": 0,
        "prs_citation": CITATION,
        "prs_status": "not_assessed",
        "prs_used_rsids": "",
        "prs_missing_rsids": "",
    }
    if weights is None or weights.empty:
        return empty
    total = 0.0
    used = 0
    used_rsids: list[str] = []
    missing_rsids: list[str] = []
    sample_key = str(sample_id)
    for row in weights.itertuples(index=False):
        key = (chrom_bare(str(row.chrom)), int(row.pos))
        rsid = str(getattr(row, "rsid", "") or "")
        record = genotypes.get(key)
        dosage = None
        if record is not None:
            samples = record.get("samples") or {}
            fmt = samples.get(sample_key)
            if fmt is None and len(samples) == 1:
                fmt = next(iter(samples.values()))
            gt = str((fmt or {}).get("GT") or "")
            dosage = _dosage(gt, str(row.effect_allele), str(record["ref"]), str(record["alt"]))
        if dosage is None:
            missing_rsids.append(rsid)
            continue
        total += dosage * float(row.weight)
        used += 1
        used_rsids.append(rsid)
    n_weights = len(weights)
    citations = sorted({str(value) for value in weights["citation"].astype(str) if str(value).strip()})
    return {
        "panel_prs": round(total, 6) if used else None,
        "prs_sites_used": used,
        "prs_sites_missing": n_weights - used,
        "prs_sites_total": n_weights,
        "prs_citation": "; ".join(citations) if citations else CITATION,
        "prs_status": "called" if used else "not_assessed",
        "prs_used_rsids": ";".join(used_rsids),
        "prs_missing_rsids": ";".join(missing_rsids),
    }
