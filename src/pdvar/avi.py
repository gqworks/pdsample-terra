"""Join AlphaGenome Atlas AVI scores onto small-variant rows.

The Atlas lookup itself is ``workflow/scripts/annotate_avi.py``. This module
only attaches the resulting table. Missing scores add no evidence points.
"""

from __future__ import annotations

import pandas as pd

from .normalize import add_variant_id

AVI_COLUMNS = ("avi_raw", "avi_quantile", "avi_phred", "avi_url")


def merge_avi(variants: pd.DataFrame, avi: pd.DataFrame | None) -> pd.DataFrame:
    """Left-join AVI columns on ``variant_id``. Unmatched rows stay empty."""
    out = variants.copy()
    for column in AVI_COLUMNS:
        if column in out.columns:
            out = out.drop(columns=column)
    if "variant_id" not in out.columns and {"chrom", "pos", "ref", "alt"}.issubset(out.columns):
        out = add_variant_id(out)
    if avi is None or avi.empty or "pos" not in avi.columns or "variant_id" not in out.columns:
        for column in AVI_COLUMNS:
            out[column] = pd.NA
        return out
    table = avi.copy()
    needed = {"chrom", "pos", "ref", "alt"}
    missing = needed.difference(table.columns)
    if missing:
        raise ValueError(f"AVI table is missing {sorted(missing)}")
    table = add_variant_id(table)
    keep = ["variant_id", *[column for column in AVI_COLUMNS if column in table.columns]]
    table = table[keep].drop_duplicates("variant_id", keep="first")
    merged = out.merge(table, on="variant_id", how="left")
    for column in ("avi_raw", "avi_quantile", "avi_phred"):
        if column in merged.columns:
            merged[column] = pd.to_numeric(merged[column], errors="coerce")
    for column in AVI_COLUMNS:
        if column not in merged.columns:
            merged[column] = pd.NA
    return merged
