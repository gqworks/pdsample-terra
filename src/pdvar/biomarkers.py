"""Phenotype and biomarker priors over gene classes.

Numeric multipliers are applied only when the rule sets both ``multiplier`` and
``citation``. Otherwise the match is recorded and the prior stays neutral
(multiplier 1). Missing phenotype columns are neutral.
"""

from __future__ import annotations

import pandas as pd

_TRUE = {"1", "true", "yes", "y", "positive", "pos", "abnormal"}
_FALSE = {"0", "false", "no", "n", "negative", "neg", "normal"}


def _field(df: pd.DataFrame, names: list[str]) -> pd.Series | None:
    for name in names:
        if name in df.columns:
            return df[name]
    return None


def _matches(value: object, rule: dict) -> bool:
    if value is None or (isinstance(value, float) and pd.isna(value)):
        return False
    op = rule.get("op", "equals")
    if op == "below":
        try:
            return float(value) < float(rule["value"])
        except (TypeError, ValueError):
            return False
    allowed = {str(item).lower() for item in rule.get("equals", [])}
    text = str(value).strip().lower()
    if text in _TRUE and "true" in allowed:
        return True
    if text in _FALSE and "false" in allowed:
        return True
    return text in allowed


def apply_biomarker_priors(df: pd.DataFrame, cfg: dict | None, early_onset_age: float | None = None) -> pd.DataFrame:
    """Add prior notes. Scores change only when a cited multiplier is configured."""
    df = df.copy()
    cfg = cfg or {}
    classes: dict[str, set[str]] = {
        name: {str(gene).upper() for gene in genes}
        for name, genes in (cfg.get("gene_classes") or {}).items()
    }
    fields = cfg.get("fields") or {}
    rules = cfg.get("rules") or []
    notes = [[] for _ in range(len(df))]
    multipliers = [1.0] * len(df)
    adjusted = [False] * len(df)
    genes = df["gene"].astype(str).str.upper() if "gene" in df.columns else pd.Series([""] * len(df), index=df.index)

    resolved_fields: dict[str, pd.Series | None] = {}
    for rule in rules:
        field_name = rule.get("field")
        if field_name not in resolved_fields:
            aliases = list(fields.get(field_name, []))
            if field_name and field_name not in aliases:
                aliases = [field_name, *aliases]
            resolved_fields[field_name] = _field(df, aliases) if aliases else None

    for rule in rules:
        series = resolved_fields.get(rule.get("field"))
        if series is None:
            continue
        prepared = dict(rule)
        if prepared.get("use_early_onset_age") and early_onset_age is not None:
            prepared["value"] = early_onset_age
            prepared["op"] = "below"
        members = classes.get(prepared.get("gene_class"), set())
        for pos, (idx, gene, value) in enumerate(zip(df.index, genes, series)):
            if members and gene not in members:
                continue
            if not _matches(value, prepared):
                continue
            notes[pos].append(f"{prepared.get('id')}:{prepared.get('direction', 'note')}")
            multiplier = prepared.get("multiplier")
            citation = prepared.get("citation")
            if multiplier is not None and citation:
                multipliers[pos] *= float(multiplier)
                adjusted[pos] = True
            elif multiplier is not None:
                notes[pos].append(f"{prepared.get('id')}:multiplier_ignored_missing_citation")

    df["biomarker_notes"] = [";".join(items) for items in notes]
    df["biomarker_prior_multiplier"] = multipliers
    df["biomarker_prior_adjusted"] = adjusted
    for canonical, aliases in fields.items():
        series = _field(df, list(aliases))
        if series is not None:
            df[f"biomarker_{canonical}"] = series
    return df
