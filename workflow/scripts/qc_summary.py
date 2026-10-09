"""Single-sample QC summary for the PdSample workflow.

Reads VerifyBamID2 selfSM, a somalier samples table, a mosdepth summary, and
the reference-check JSON. Relatedness pairs are ignored. The standard library
only, so the WDL task can run this on python:3.11-slim.
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path


def _flag(value: str) -> bool:
    return str(value).strip().lower() in {"1", "true", "yes"}


def normalize_sex(value: object) -> str:
    """Map free text to male, female, or unknown."""
    text = "" if value is None else str(value).strip().lower()
    if text in {"male", "m", "xy", "1"}:
        return "male"
    if text in {"female", "f", "xx", "2"}:
        return "female"
    return "unknown"


def _rows(path: Path) -> list[dict[str, str]]:
    """Tab-separated rows. A leading ``#`` on the header is stripped.

    VerifyBamID selfSM headers look like ``#SEQ_ID\\tFREEMIX``. Comment lines
    without a tab are ignored.
    """
    lines: list[str] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line.strip():
            continue
        if line.startswith("#"):
            stripped = line[1:].lstrip()
            if not lines and "\t" in stripped:
                lines.append(stripped)
            continue
        lines.append(line)
    if not lines:
        return []
    return list(csv.DictReader(lines, delimiter="\t"))


def _column(row: dict[str, str], *names: str) -> str | None:
    lookup = {key.lower(): key for key in row}
    for name in names:
        key = lookup.get(name.lower())
        if key is not None:
            return row[key]
    return None


def parse_freemix(path: Path) -> float | None:
    """Largest FREEMIX value in a VerifyBamID2 selfSM table."""
    values = []
    for row in _rows(path):
        raw = _column(row, "FREEMIX")
        if raw in (None, "", "."):
            continue
        try:
            values.append(float(raw))
        except ValueError:
            continue
    return max(values) if values else None


def parse_somalier_sex(path: Path) -> str:
    """Predicted sex from somalier, else an X-heterozygosity fallback."""
    rows = _rows(path)
    if not rows:
        return "unknown"
    row = rows[0]
    for name in ("predicted_sex", "sex"):
        raw = _column(row, name)
        if raw is None:
            continue
        sex = normalize_sex(raw)
        if sex != "unknown":
            return sex
    raw_n = _column(row, "x_n", "x_count")
    raw_het = _column(row, "x_het")
    if raw_n is None or raw_het is None:
        return "unknown"
    try:
        total = float(raw_n)
        het = float(raw_het)
    except ValueError:
        return "unknown"
    if total < 10:
        return "unknown"
    rate = het / total
    if rate > 0.25:
        return "female"
    if rate < 0.05:
        return "male"
    return "unknown"


def parse_mean_coverage(path: Path) -> float | None:
    """Mean depth from a mosdepth summary. Prefer the total line."""
    rows = _rows(path)
    if not rows:
        return None
    chosen = rows[0]
    for row in rows:
        chrom = (_column(row, "chrom", "chromosome") or "").lower()
        if chrom in {"total", "total_region"}:
            chosen = row
            break
    raw = _column(chosen, "mean")
    if raw in (None, "", "."):
        return None
    try:
        return float(raw)
    except ValueError:
        return None


def build_summary(
    *,
    sample_id: str,
    selfsm: Path | None,
    somalier: Path | None,
    mosdepth: Path | None,
    reference: Path | None,
    provided_sex: str,
    max_freemix: float,
    min_mean_coverage: float,
    hard_fail_contamination: bool,
    hard_fail_sex: bool,
    hard_fail_coverage: bool,
) -> tuple[dict, str]:
    """Return ``(qc_object, sex_used)``. ``qc_object['pass']`` is false when a hard fail trips."""
    reference_doc: dict = {}
    if reference and reference.exists():
        reference_doc = json.loads(reference.read_text(encoding="utf-8"))
    freemix = parse_freemix(selfsm) if selfsm and selfsm.exists() else None
    somalier_sex = parse_somalier_sex(somalier) if somalier and somalier.exists() else "unknown"
    provided = normalize_sex(provided_sex)
    if provided != "unknown":
        sex, source = provided, "provided"
    elif somalier_sex != "unknown":
        sex, source = somalier_sex, "somalier"
    else:
        sex, source = "unknown", "unknown"
    mean = parse_mean_coverage(mosdepth) if mosdepth and mosdepth.exists() else None
    contamination_pass = freemix is None or freemix <= max_freemix
    coverage_pass = mean is None or min_mean_coverage <= 0 or mean >= min_mean_coverage
    sex_discordant = provided != "unknown" and somalier_sex != "unknown" and provided != somalier_sex
    sex_pass = not sex_discordant and not (hard_fail_sex and sex == "unknown")
    failures = []
    if hard_fail_contamination and not contamination_pass:
        failures.append(f"FREEMIX {freemix} exceeds {max_freemix}")
    if hard_fail_coverage and not coverage_pass:
        failures.append(f"mean coverage {mean} is below {min_mean_coverage}")
    if hard_fail_sex and not sex_pass:
        failures.append(f"sex check failed (provided={provided}, somalier={somalier_sex})")
    summary = {
        "sample_id": sample_id,
        "pass": not failures,
        "failures": failures,
        "sex": {"value": sex, "source": source, "somalier": somalier_sex, "provided": provided},
        "contamination": {
            "freemix": freemix,
            "max_freemix": max_freemix,
            "pass": contamination_pass,
        },
        "coverage": {
            "mean": mean,
            "min_mean_coverage": min_mean_coverage,
            "pass": coverage_pass,
        },
        "reference": {
            "ok": bool(reference_doc.get("ok", True)) if reference_doc else None,
            "cram_contigs": reference_doc.get("cram_contigs"),
            "reference_contigs": reference_doc.get("reference_contigs"),
        },
        "relatedness": "not_run",
    }
    return summary, sex


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--selfsm", required=True)
    parser.add_argument("--somalier", required=True)
    parser.add_argument("--mosdepth", required=True)
    parser.add_argument("--reference", required=True)
    parser.add_argument("--provided-sex", default="")
    parser.add_argument("--max-freemix", type=float, default=0.02)
    parser.add_argument("--min-mean-coverage", type=float, default=0)
    parser.add_argument("--hard-fail-contamination", default="true")
    parser.add_argument("--hard-fail-sex", default="false")
    parser.add_argument("--hard-fail-coverage", default="false")
    parser.add_argument("--out", required=True)
    parser.add_argument("--sex-out", required=True)
    args = parser.parse_args(argv)
    summary, sex = build_summary(
        sample_id=args.sample_id,
        selfsm=Path(args.selfsm),
        somalier=Path(args.somalier),
        mosdepth=Path(args.mosdepth),
        reference=Path(args.reference),
        provided_sex=args.provided_sex,
        max_freemix=args.max_freemix,
        min_mean_coverage=args.min_mean_coverage,
        hard_fail_contamination=_flag(args.hard_fail_contamination),
        hard_fail_sex=_flag(args.hard_fail_sex),
        hard_fail_coverage=_flag(args.hard_fail_coverage),
    )
    Path(args.out).write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    Path(args.sex_out).write_text(sex + "\n", encoding="utf-8")
    if not summary["pass"]:
        print("; ".join(summary["failures"]))
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
