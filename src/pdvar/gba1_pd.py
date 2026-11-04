"""GBA1-PD Browser lookup for Gauchian catalog SNVs.

Source: Parlar et al. Mov Disord. 2023 (GBA1-PD Browser). Used as triage
annotation and a score bump on ``hit_source=gauchian`` alleles only.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from .callsets import gba_c_dot
from .config import PROJECT_ROOT, resolve
from .scoring import assign_tier

ANNOTATION_COLUMNS = (
    "gba1_pd_class",
    "gba1_gd_acmg",
    "gba1_pd_variant_name",
    "gba1_pd_rsid",
    "gba1_pd_matched_cdna",
)

CLASS_REASON = {
    "severe": "gba1_pd_severe",
    "mild": "gba1_pd_mild",
    "risk": "gba1_pd_risk",
}


def default_browser_path(cfg: dict | None = None) -> Path | None:
    """Resolve the shipped browser TSV from config or the repo default."""
    cfg = cfg or {}
    configured = (
        (cfg.get("inputs") or {}).get("gba1_pd_browser")
        or (cfg.get("paths") or {}).get("gba1_pd_browser")
    )
    candidates: list[Path] = []
    if configured:
        raw = Path(str(configured))
        if raw.is_absolute():
            candidates.append(raw)
        else:
            root = Path(str(cfg.get("_root") or PROJECT_ROOT))
            candidates.extend([root / raw, resolve(raw), Path.cwd() / raw])
    candidates.append(PROJECT_ROOT / "resources" / "gba1" / "gba1_pd_browser.tsv")
    for path in candidates:
        if path.is_file():
            return path
    return None


def load_browser(path: str | Path | None) -> pd.DataFrame:
    """Load and index the cleaned browser table by normalized c. HGVS."""
    if path is None:
        return pd.DataFrame()
    path = Path(path)
    if not path.is_file():
        return pd.DataFrame()
    table = pd.read_csv(path, sep="\t", dtype=str).fillna("")
    if table.empty or "cdna" not in table.columns:
        return pd.DataFrame()
    table = table.copy()
    table["cdna_key"] = table["cdna"].map(gba_c_dot)
    table = table[table["cdna_key"].astype(bool)]
    table["pd_class"] = (
        table.get("pd_class", pd.Series("", index=table.index))
        .astype(str)
        .str.strip()
        .str.lower()
    )
    return table.drop_duplicates(subset=["cdna_key"], keep="first")


def annotate_gauchian(df: pd.DataFrame, browser: pd.DataFrame | None) -> pd.DataFrame:
    """Left-join browser fields onto Gauchian GBA1 SNVs matched by c. HGVS."""
    df = df.copy()
    for column in ANNOTATION_COLUMNS:
        if column not in df.columns:
            df[column] = pd.NA
    if df.empty or browser is None or browser.empty:
        return df
    if "hit_source" not in df.columns or "gene" not in df.columns:
        return df
    gene = df["gene"].astype(str).str.upper()
    hit = df["hit_source"].fillna("").astype(str)
    gauchian = gene.eq("GBA1") & hit.eq("gauchian")
    if not gauchian.any():
        return df
    hgvsc = df["hgvsc"] if "hgvsc" in df.columns else pd.Series("", index=df.index)
    keys = hgvsc.map(gba_c_dot)
    lookup = browser.set_index("cdna_key", drop=False)
    for idx in df.index[gauchian]:
        key = keys.at[idx]
        if not key or key not in lookup.index:
            continue
        row = lookup.loc[key]
        if isinstance(row, pd.DataFrame):
            row = row.iloc[0]
        df.at[idx, "gba1_pd_class"] = str(row.get("pd_class") or "") or pd.NA
        df.at[idx, "gba1_gd_acmg"] = str(row.get("gd_acmg") or "") or pd.NA
        df.at[idx, "gba1_pd_variant_name"] = str(row.get("variant_name") or "") or pd.NA
        df.at[idx, "gba1_pd_rsid"] = str(row.get("rsid") or "") or pd.NA
        df.at[idx, "gba1_pd_matched_cdna"] = str(row.get("cdna") or key) or pd.NA
    return df


def apply_evidence(df: pd.DataFrame, scfg: dict) -> pd.DataFrame:
    """Add configured score bumps for Gauchian rows with a browser PD class."""
    if df.empty or "evidence_score" not in df.columns:
        return df
    if "gba1_pd_class" not in df.columns or "hit_source" not in df.columns:
        return df
    df = df.copy()
    points = scfg.get("points") or {}
    hit = df["hit_source"].fillna("").astype(str)
    classes = df["gba1_pd_class"].fillna("").astype(str).str.strip().str.lower()
    reasons = df["evidence_reasons"].fillna("").astype(str)
    for pd_class, reason in CLASS_REASON.items():
        pts = int(points.get(reason) or 0)
        if not pts:
            continue
        mask = hit.eq("gauchian") & classes.eq(pd_class) & ~reasons.str.contains(reason, regex=False)
        if not mask.any():
            continue
        df.loc[mask, "evidence_score"] = df.loc[mask, "evidence_score"] + pts
        df.loc[mask, "evidence_reasons"] = (reasons.loc[mask] + f";{reason}").str.strip(";")
        reasons = df["evidence_reasons"].fillna("").astype(str)
    if "tiers" in scfg:
        df["tier"] = assign_tier(df["evidence_score"], scfg["tiers"])
    return df


def annotate_and_score(
    df: pd.DataFrame,
    cfg: dict | None = None,
    *,
    browser_path: str | Path | None = None,
) -> pd.DataFrame:
    """Resolve the browser table, annotate Gauchian rows, and apply score bumps."""
    cfg = cfg or {}
    path = Path(browser_path) if browser_path else default_browser_path(cfg)
    browser = load_browser(path)
    scored = annotate_gauchian(df, browser)
    return apply_evidence(scored, cfg.get("scoring") or {})
