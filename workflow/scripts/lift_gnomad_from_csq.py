#!/usr/bin/env python3
"""Lift gnomAD AF from dbNSFP CSQ fields into INFO/gnomAD_AF and INFO/AF_popmax.

Academic dbNSFP 5.4a annotates ``gnomAD4.1_joint_AF`` / ``gnomAD4.1_joint_POPMAX_AF``
in VEP CSQ. Introme's rare filter looks for INFO keys, so this copies those CSQ
values (preferring already-present INFO) onto each record.
"""

from __future__ import annotations

import argparse
import gzip
import re
import sys
from pathlib import Path

_CSQ_FORMAT = re.compile(r"Format:\s*([^\"\s]+)")
_AF_CSQ = ("gnomAD4.1_joint_AF", "gnomAD_AF", "gnomADe_AF", "gnomADg_AF", "MAX_AF")
_POP_CSQ = (
    "gnomAD4.1_joint_POPMAX_AF",
    "gnomAD_AF_popmax",
    "gnomADe_AF_popmax",
    "gnomADg_AF_popmax",
    "MAX_AF_POPS",
)
_AF_HEADER = (
    '##INFO=<ID=gnomAD_AF,Number=1,Type=Float,'
    'Description="gnomAD AF from dbNSFP CSQ (gnomAD4.1_joint_AF) or prior INFO">'
)
_POP_HEADER = (
    '##INFO=<ID=AF_popmax,Number=1,Type=Float,'
    'Description="gnomAD popmax AF from dbNSFP CSQ (gnomAD4.1_joint_POPMAX_AF) or prior INFO">'
)


def _open(path: Path, mode: str = "rt"):
    if str(path).endswith(".gz"):
        return gzip.open(path, mode, encoding="utf-8" if "t" in mode else None)
    return path.open(mode, encoding="utf-8" if "t" in mode else None)


def _csq_fields_from_header(header_lines: list[str]) -> list[str]:
    for line in header_lines:
        if line.startswith("##INFO=<ID=CSQ,"):
            match = _CSQ_FORMAT.search(line)
            if match:
                return match.group(1).split("|")
    return []


def _info_map(info: str) -> dict[str, str]:
    out: dict[str, str] = {}
    if info in {".", ""}:
        return out
    for part in info.split(";"):
        if "=" in part:
            key, value = part.split("=", 1)
            out[key] = value
        elif part:
            out[part] = "."
    return out


def _first_float_token(value: str) -> str | None:
    if not value or value in {".", "NA", "na", "None"}:
        return None
    for token in re.split(r"[&;,/]", value):
        text = token.strip()
        if not text or text in {".", "NA", "na"}:
            continue
        try:
            float(text)
            return text
        except ValueError:
            continue
    return None


def _csq_value(csq_blob: str, fields: list[str], names: tuple[str, ...]) -> str | None:
    if not csq_blob or not fields:
        return None
    item = csq_blob.split(",")[0]
    parts = item.split("|")
    mapped = {name: (parts[i] if i < len(parts) else "") for i, name in enumerate(fields)}
    for name in names:
        token = _first_float_token(mapped.get(name, ""))
        if token is not None:
            return token
    return None


def _ensure_headers(header_lines: list[str]) -> list[str]:
    has_af = any(line.startswith("##INFO=<ID=gnomAD_AF,") for line in header_lines)
    has_pop = any(line.startswith("##INFO=<ID=AF_popmax,") for line in header_lines)
    chrom_idx = next((i for i, line in enumerate(header_lines) if line.startswith("#CHROM")), len(header_lines))
    insert_at = chrom_idx
    extras: list[str] = []
    if not has_af:
        extras.append(_AF_HEADER)
    if not has_pop:
        extras.append(_POP_HEADER)
    return header_lines[:insert_at] + extras + header_lines[insert_at:]


def lift(in_vcf: Path, out_vcf: Path) -> None:
    with _open(in_vcf, "rt") as src:
        header_lines: list[str] = []
        body_lines: list[str] = []
        for line in src:
            if line.startswith("#"):
                header_lines.append(line.rstrip("\n"))
            else:
                body_lines.append(line)
                body_lines.extend(src)
                break

    fields = _csq_fields_from_header(header_lines)
    if not fields:
        print("WARNING: no CSQ Format in header; copying VCF unchanged aside from INFO headers", file=sys.stderr)
    out_header = _ensure_headers(header_lines)

    with _open(out_vcf, "wt") as dest:
        for line in out_header:
            dest.write(line + "\n")
        for line in body_lines:
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 8:
                dest.write(line if line.endswith("\n") else line + "\n")
                continue
            info = _info_map(cols[7])
            af = info.get("gnomAD_AF")
            pop = info.get("AF_popmax")
            if af in {None, "", "."}:
                af = _csq_value(info.get("CSQ", ""), fields, _AF_CSQ)
            if pop in {None, "", "."}:
                pop = _csq_value(info.get("CSQ", ""), fields, _POP_CSQ)
            parts = [] if cols[7] in {".", ""} else cols[7].split(";")
            rebuilt: list[str] = []
            seen_af = seen_pop = False
            for part in parts:
                key = part.split("=", 1)[0]
                if key == "gnomAD_AF" and af not in {None, "", "."}:
                    rebuilt.append(f"gnomAD_AF={af}")
                    seen_af = True
                elif key == "AF_popmax" and pop not in {None, "", "."}:
                    rebuilt.append(f"AF_popmax={pop}")
                    seen_pop = True
                else:
                    rebuilt.append(part)
            if not seen_af and af not in {None, "", "."}:
                rebuilt.append(f"gnomAD_AF={af}")
            if not seen_pop and pop not in {None, "", "."}:
                rebuilt.append(f"AF_popmax={pop}")
            cols[7] = ";".join(rebuilt) if rebuilt else "."
            dest.write("\t".join(cols) + "\n")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--in-vcf", required=True)
    parser.add_argument("--out-vcf", required=True)
    args = parser.parse_args()
    lift(Path(args.in_vcf), Path(args.out_vcf))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
