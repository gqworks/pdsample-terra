#!/usr/bin/env python3
"""Annotate a sample VCF with panel-sliced ClinVar P/LP sites.

Looks up each sample allele by CHROM/POS/REF/ALT (and, for SVs, ID) in the
sites file and writes ``INFO/ClinVar_CLNSIG`` from ClinVar ``CLNSIG``. Matching
sites overwrite a prior ClinVar field so the panel slice wins over dbNSFP CSQ
for indels and SVs that table does not cover.
"""

from __future__ import annotations

import argparse
import gzip
import sys
from pathlib import Path

_CLINVAR_HEADER = (
    '##INFO=<ID=ClinVar_CLNSIG,Number=1,Type=String,'
    'Description="ClinVar CLNSIG from panel-sliced P/LP sites">'
)


def _is_gzip_path(path: Path) -> bool:
    name = str(path).lower()
    if name.endswith((".gz", ".bgz")):
        return True
    if not path.is_file():
        return False
    try:
        with path.open("rb") as handle:
            return handle.read(2) == b"\x1f\x8b"
    except OSError:
        return False


def _open_text(path: Path, mode: str = "rt"):
    encoding = "utf-8" if "t" in mode else None
    use_gzip = (
        str(path).lower().endswith((".gz", ".bgz"))
        if "w" in mode
        else _is_gzip_path(path)
    )
    if use_gzip:
        return gzip.open(path, mode, encoding=encoding)
    return path.open(mode, encoding=encoding)


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


def is_plp(clnsig: str) -> bool:
    text = str(clnsig or "").replace("_", " ").lower()
    if not text or "conflicting" in text:
        return False
    if "benign" in text and "pathogenic" not in text:
        return False
    return "pathogenic" in text


def load_sites(path: Path) -> tuple[dict[tuple[str, int, str, str], str], dict[str, str]]:
    """Index P/LP sites by allele and by VCF ID."""
    alleles: dict[tuple[str, int, str, str], str] = {}
    ids: dict[str, str] = {}
    with _open_text(path, "rt") as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 8:
                continue
            chrom, pos_text, vid, ref, alts = cols[0], cols[1], cols[2], cols[3], cols[4]
            try:
                pos = int(pos_text)
            except ValueError:
                continue
            info = _info_map(cols[7])
            clnsig = info.get("CLNSIG") or info.get("ClinVar_CLNSIG") or ""
            if not is_plp(clnsig):
                continue
            chrom_key = bare_chrom(chrom)
            ref_u = ref.upper()
            for alt in alts.split(","):
                alt_u = alt.upper()
                if not alt_u or alt_u == ".":
                    continue
                alleles[(chrom_key, pos, ref_u, alt_u)] = clnsig
            if vid and vid not in {".", ""}:
                ids[vid] = clnsig
    return alleles, ids


def _ensure_headers(header_lines: list[str]) -> list[str]:
    has = any(line.startswith("##INFO=<ID=ClinVar_CLNSIG,") for line in header_lines)
    chrom_idx = next(
        (i for i, line in enumerate(header_lines) if line.startswith("#CHROM")),
        len(header_lines),
    )
    extras = [] if has else [_CLINVAR_HEADER]
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
    alleles, ids = load_sites(sites_vcf)
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
            chrom, pos_text, vid, ref, alts = cols[0], cols[1], cols[2], cols[3], cols[4]
            try:
                pos = int(pos_text)
            except ValueError:
                dest.write("\t".join(cols) + "\n")
                continue
            chrom_key = bare_chrom(chrom)
            clnsig = ids.get(vid) if vid and vid not in {".", ""} else None
            if clnsig is None:
                for alt in alts.split(","):
                    hit = alleles.get((chrom_key, pos, ref.upper(), alt.upper()))
                    if hit:
                        clnsig = hit
                        break
            if not clnsig:
                dest.write("\t".join(cols) + "\n")
                continue
            matched += 1
            cols[7] = _set_info_fields(cols[7], {"ClinVar_CLNSIG": clnsig})
            dest.write("\t".join(cols) + "\n")
    return matched


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--in-vcf", required=True)
    parser.add_argument("--sites-vcf", required=True)
    parser.add_argument("--out-vcf", required=True)
    args = parser.parse_args()
    matched = annotate(Path(args.in_vcf), Path(args.sites_vcf), Path(args.out_vcf))
    print(f"annotate_clinvar_plp: matched {matched} site(s)", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
