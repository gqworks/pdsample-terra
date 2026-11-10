#!/usr/bin/env python3
"""Annotate a sample VCF with gnomAD AF from a panel-sliced sites VCF.

Looks up each sample allele by CHROM/POS/REF/ALT in the sites file and writes
``INFO/gnomAD_AF`` / ``INFO/AF_popmax`` from joint (``AF_joint`` /
``AF_grpmax_joint``) or standard (``AF`` / ``AF_grpmax``) INFO fields.

When a sites allele matches, those INFO keys are set (overwriting prior values)
so true gnomAD AF wins over dbNSFP CSQ lifts. Alleles absent from the sites
file are left unchanged for ``LiftGnomadInfo`` / dbNSFP fallback.
"""

from __future__ import annotations

import argparse
import gzip
import sys
from pathlib import Path

_AF_KEYS = ("AF_joint", "AF", "gnomAD_AF")
_POP_KEYS = ("AF_grpmax_joint", "AF_grpmax", "AF_popmax", "gnomAD_AF_popmax")
_AF_HEADER = (
    '##INFO=<ID=gnomAD_AF,Number=A,Type=Float,'
    'Description="gnomAD AF from panel-sliced sites (AF_joint or AF)">'
)
_POP_HEADER = (
    '##INFO=<ID=AF_popmax,Number=A,Type=Float,'
    'Description="gnomAD grpmax/popmax AF from panel-sliced sites">'
)


def _open_text(path: Path, mode: str = "rt"):
    if str(path).endswith(".gz"):
        return gzip.open(path, mode, encoding="utf-8" if "t" in mode else None)
    return path.open(mode, encoding="utf-8" if "t" in mode else None)


def bare_chrom(chrom: str) -> str:
    text = str(chrom).strip()
    rest = text[3:] if text.lower().startswith("chr") else text
    return "M" if rest.upper() in {"M", "MT"} else rest


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


def _allele_token(value: str, allele_index: int) -> str | None:
    if not value or value in {".", "NA", "na"}:
        return None
    parts = value.split(",")
    if allele_index < len(parts):
        token = parts[allele_index].strip()
    else:
        token = parts[0].strip()
    if not token or token in {".", "NA", "na"}:
        return None
    try:
        float(token)
    except ValueError:
        return None
    return token


def _pick_af(info: dict[str, str], keys: tuple[str, ...], allele_index: int) -> str | None:
    for key in keys:
        token = _allele_token(info.get(key, ""), allele_index)
        if token is not None:
            return token
    return None


def load_sites(path: Path) -> dict[tuple[str, int, str, str], tuple[str | None, str | None]]:
    """Index sites alleles as (bare_chrom, pos, ref, alt) -> (af, popmax)."""
    table: dict[tuple[str, int, str, str], tuple[str | None, str | None]] = {}
    with _open_text(path, "rt") as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 8:
                continue
            chrom, pos_text, _vid, ref, alts = cols[0], cols[1], cols[2], cols[3], cols[4]
            try:
                pos = int(pos_text)
            except ValueError:
                continue
            info = _info_map(cols[7])
            ref_u = ref.upper()
            chrom_key = bare_chrom(chrom)
            for allele_i, alt in enumerate(alts.split(",")):
                alt_u = alt.upper()
                if not alt_u or alt_u == ".":
                    continue
                af = _pick_af(info, _AF_KEYS, allele_i)
                pop = _pick_af(info, _POP_KEYS, allele_i)
                if af is None and pop is None:
                    continue
                table[(chrom_key, pos, ref_u, alt_u)] = (af, pop)
    return table


def _ensure_headers(header_lines: list[str]) -> list[str]:
    has_af = any(line.startswith("##INFO=<ID=gnomAD_AF,") for line in header_lines)
    has_pop = any(line.startswith("##INFO=<ID=AF_popmax,") for line in header_lines)
    chrom_idx = next(
        (i for i, line in enumerate(header_lines) if line.startswith("#CHROM")),
        len(header_lines),
    )
    extras: list[str] = []
    if not has_af:
        extras.append(_AF_HEADER)
    if not has_pop:
        extras.append(_POP_HEADER)
    return header_lines[:chrom_idx] + extras + header_lines[chrom_idx:]


def _set_info_fields(info: str, updates: dict[str, str]) -> str:
    parts = [] if info in {".", ""} else info.split(";")
    rebuilt: list[str] = []
    seen = set()
    for part in parts:
        key = part.split("=", 1)[0]
        if key in updates:
            rebuilt.append(f"{key}={updates[key]}")
            seen.add(key)
        else:
            rebuilt.append(part)
    for key, value in updates.items():
        if key not in seen:
            rebuilt.append(f"{key}={value}")
    return ";".join(rebuilt) if rebuilt else "."


def annotate(in_vcf: Path, sites_vcf: Path, out_vcf: Path) -> int:
    sites = load_sites(sites_vcf)
    matched = 0
    with _open_text(in_vcf, "rt") as src:
        header_lines: list[str] = []
        body_lines: list[str] = []
        for line in src:
            if line.startswith("#"):
                header_lines.append(line.rstrip("\n"))
            else:
                body_lines.append(line)
                body_lines.extend(src)
                break

    out_header = _ensure_headers(header_lines)
    with _open_text(out_vcf, "wt") as dest:
        for line in out_header:
            dest.write(line + "\n")
        for line in body_lines:
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 8:
                dest.write(line if line.endswith("\n") else line + "\n")
                continue
            chrom, pos_text, ref, alts = cols[0], cols[1], cols[3], cols[4]
            try:
                pos = int(pos_text)
            except ValueError:
                dest.write("\t".join(cols) + "\n")
                continue
            chrom_key = bare_chrom(chrom)
            ref_u = ref.upper()
            alt_list = [a.upper() for a in alts.split(",")]
            af_tokens: list[str] = []
            pop_tokens: list[str] = []
            any_hit = False
            for alt_u in alt_list:
                hit = sites.get((chrom_key, pos, ref_u, alt_u))
                if hit is None:
                    af_tokens.append(".")
                    pop_tokens.append(".")
                    continue
                any_hit = True
                af, pop = hit
                af_tokens.append(af if af is not None else ".")
                pop_tokens.append(pop if pop is not None else ".")
            if not any_hit:
                dest.write("\t".join(cols) + "\n")
                continue
            matched += 1
            updates: dict[str, str] = {}
            if any(tok != "." for tok in af_tokens):
                updates["gnomAD_AF"] = ",".join(af_tokens)
            if any(tok != "." for tok in pop_tokens):
                updates["AF_popmax"] = ",".join(pop_tokens)
            if updates:
                cols[7] = _set_info_fields(cols[7], updates)
            dest.write("\t".join(cols) + "\n")
    return matched


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--in-vcf", required=True)
    parser.add_argument("--sites-vcf", required=True)
    parser.add_argument("--out-vcf", required=True)
    args = parser.parse_args()
    matched = annotate(Path(args.in_vcf), Path(args.sites_vcf), Path(args.out_vcf))
    print(f"annotate_gnomad_af: matched {matched} site(s)", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
