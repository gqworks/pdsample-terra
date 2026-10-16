"""Parse the GP2-style phenotype text used as a per-sample report."""

from __future__ import annotations

import re

_NA = {"", "na", "n/a", "nan", "none", "null", "."}

_FIELD_MAP = {
    "gp2id": "gp2_id",
    "studycode_id": "study_id",
    "sex": "reported_sex",
    "age at sample collection": "age_sampling",
    "pd status": "pd_status",
    "levodopa": "levodopa",
    "age at levodopa start": "age_levodopa_start",
    "age at diagnosis": "age_diagnosis",
    "diagnosis type": "diagnosis_type",
    "age at onset": "age_onset",
}

_BULLET = re.compile(r"^-\s+(.+?)\s*\(age:\s*([^)]+)\)\s*$", re.IGNORECASE)


def looks_like_phenotype_report(text: str) -> bool:
    """True for the labelled sheet, so it is not searched as free-text symptoms."""
    lowered = text.lower()
    if "exposures and comorbidities" in lowered:
        return True
    return "gp2id:" in lowered and "age at onset:" in lowered


def parse_phenotype_report(text: str | None) -> dict:
    """Read age at onset, levodopa, sampling age, and the comorbidity list.

    ``NA`` values are omitted. Levodopa is kept as recorded. It is not rewritten
    to the profile phrase ``levodopa responsive``.
    """
    if text is None or not str(text).strip():
        return {}
    parsed: dict = {}
    comorbidities: list[dict] = []
    relevance: str | None = None
    for raw_line in str(text).splitlines():
        line = raw_line.strip()
        if not line:
            continue
        header = line.lower().rstrip(":")
        if header == "exposures and comorbidities":
            continue
        if header in {"high clinical relevance", "high relevance", "higher relevance"}:
            relevance = "high"
            continue
        if header in {"lower relevance", "low relevance", "lower clinical relevance"}:
            relevance = "lower"
            continue
        bullet = _BULLET.match(line)
        if bullet:
            item = {"condition": bullet.group(1).strip(), "relevance": relevance or "unspecified"}
            age = _number(bullet.group(2))
            if age is not None:
                item["age"] = age
            comorbidities.append(item)
            continue
        if ":" not in line:
            continue
        for label, value in _pairs(line):
            field = _FIELD_MAP.get(label)
            if field is None:
                continue
            token = _token(value)
            if token is None:
                continue
            if field == "reported_sex":
                token = _sex(token)
            parsed[field] = token
    if comorbidities:
        parsed["comorbidities"] = comorbidities
    if parsed:
        parsed["phenotype_source"] = "report"
    return parsed


def _pairs(line: str) -> list[tuple[str, str]]:
    pairs: list[tuple[str, str]] = []
    for part in line.split("|"):
        if ":" not in part:
            continue
        label, value = part.split(":", 1)
        pairs.append((label.strip().lower(), value.strip()))
    return pairs


def _token(value: str) -> str | None:
    text = value.strip()
    if text.lower() in _NA:
        return None
    return text


def _number(value: str) -> int | float | None:
    token = _token(value)
    if token is None:
        return None
    try:
        number = float(token)
    except ValueError:
        return None
    if number.is_integer():
        return int(number)
    return number


def _sex(value: str) -> str:
    lowered = value.strip().lower()
    if lowered in {"m", "male"}:
        return "male"
    if lowered in {"f", "female"}:
        return "female"
    return value.strip()
