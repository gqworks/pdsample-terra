"""Readers for pipeline inputs and tables."""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from .config import resolve


def read_table(path: str | Path) -> pd.DataFrame:
    """Read a TSV/CSV/Excel/Parquet file based on its extension."""
    path = Path(path)
    suffixes = "".join(path.suffixes).lower()
    if suffixes.endswith(".parquet"):
        return pd.read_parquet(path)
    if suffixes.endswith((".xlsx", ".xls")):
        return pd.read_excel(path)
    if suffixes.endswith((".csv", ".csv.gz")):
        return pd.read_csv(path, low_memory=False)
    return pd.read_csv(path, sep="\t", low_memory=False)


def find_inputs(cfg: dict, kind: str) -> list[Path]:
    """Return files matching the configured glob (relative to project root, or absolute)."""
    pattern = Path(cfg["inputs"][kind])
    return sorted(resolve(pattern.parent).glob(pattern.name))


def load_raw(cfg: dict, kind: str) -> pd.DataFrame:
    """Load and concatenate all raw files for an input kind (``introme``, ``panel``, ...)."""
    files = find_inputs(cfg, kind)
    if not files:
        raise FileNotFoundError(
            f"No {kind} files matching '{cfg['inputs'][kind]}'"
        )
    frames = []
    for f in files:
        df = read_table(f)
        df["source_file"] = f.name
        frames.append(df)
    return pd.concat(frames, ignore_index=True)
