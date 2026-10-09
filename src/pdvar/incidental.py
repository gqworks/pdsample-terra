"""Tags that are neither a Parkinson's diagnosis nor a benign dismissal.

TTN truncating variants are common incidental findings. They stay visible and
are not relabelled benign.
"""

from __future__ import annotations

import pandas as pd

from .scoring import is_null_variant

TTN_GENES = {"TTN"}
INCIDENTAL_TTN = "TTN_truncation"


def ttn_truncation_mask(df: pd.DataFrame) -> pd.Series:
    """True for null (truncating or canonical-splice) alleles in TTN."""
    if df.empty or "gene" not in df.columns or "consequence" not in df.columns:
        return pd.Series(False, index=df.index)
    gene = df["gene"].astype(str).str.upper().isin(TTN_GENES)
    return gene & is_null_variant(df["consequence"])


def apply_incidental_tags(df: pd.DataFrame) -> pd.DataFrame:
    """Set ``incidental`` and an inheritance status that is not benign.

    ClinVar significance is left as reported. A TTN truncation is not forced
    to B/LB and is not left as a panel disease call.
    """
    df = df.copy()
    if "incidental" not in df.columns:
        df["incidental"] = pd.NA
    mask = ttn_truncation_mask(df)
    df.loc[mask, "incidental"] = INCIDENTAL_TTN
    if "inheritance_status" in df.columns:
        df.loc[mask, "inheritance_status"] = "incidental"
    if "disposition" not in df.columns:
        df["disposition"] = pd.NA
    df.loc[mask, "disposition"] = "incidental"
    return df
