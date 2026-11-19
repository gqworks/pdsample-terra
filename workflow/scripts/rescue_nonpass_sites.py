#!/usr/bin/env python3
"""Restore non-PASS genotypes at ClinVar P/LP or GBA1-PD Browser sites.

Rewrites FILTER=PASS and sets INFO/PDVAR_RESCUE on carrier records that hard
filters dropped, so NormPass / VEP / Introme still see them. Alleles the
caller never emitted stay absent.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import sys
from pathlib import Path

_RESCUE_HEADER = (
    '##INFO=<ID=PDVAR_RESCUE,Number=0,Type=Flag,'
    'Description="Non-PASS site restored for ClinVar P/LP or GBA1-PD Browser genotype">'
)
_GBA1_CHROM = "1"


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


def load_clinvar(path: Path | None) -> set[tuple[str, int, str, str]]:
    sites: set[tuple[str, int, str, str]] = set()
    if path is None or not path.is_file():
        return sites
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
            clnsig = info.get("CLNSIG") or info.get("ClinVar_CLNSIG") or ""
            if not is_plp(clnsig):
                continue
            chrom_key = bare_chrom(chrom)
            for alt in alts.split(","):
                alt_u = alt.upper()
                if not alt_u or alt_u == ".":
                    continue
                sites.add((chrom_key, pos, ref.upper(), alt_u))
    return sites


def load_gba1_positions(path: Path | None) -> set[tuple[str, int]]:
    positions: set[tuple[str, int]] = set()
    if path is None or not path.is_file():
        return positions
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            raw = str((row or {}).get("pos_hg38") or "").strip()
            if not raw:
                continue
            try:
                pos = int(float(raw))
            except ValueError:
                continue
            positions.add((_GBA1_CHROM, pos))
    return positions


def _is_carrier(gt: str) -> bool:
    if not gt or gt in {".", "./.", ".|."}:
        return False
    alleles = [part for part in gt.replace("|", "/").split("/") if part not in {"", "."}]
    return any(allele not in {"0"} for allele in alleles)


def _sample_is_carrier(cols: list[str], format_field: str) -> bool:
    if len(cols) < 10:
        return True
    keys = format_field.split(":")
    if "GT" not in keys:
        return True
    gt_i = keys.index("GT")
    for sample in cols[9:]:
        parts = sample.split(":")
        gt = parts[gt_i] if gt_i < len(parts) else ""
        if _is_carrier(gt):
            return True
    return False


def _ensure_headers(header_lines: list[str]) -> list[str]:
    has = any(line.startswith("##INFO=<ID=PDVAR_RESCUE,") for line in header_lines)
    chrom_idx = next(
        (i for i, line in enumerate(header_lines) if line.startswith("#CHROM")),
        len(header_lines),
    )
    extras = [] if has else [_RESCUE_HEADER]
    return header_lines[:chrom_idx] + extras + header_lines[chrom_idx:]


def _set_flag(info: str, flag: str) -> str:
    parts = [] if info in {".", ""} else info.split(";")
    keys = {part.split("=", 1)[0] for part in parts if part}
    if flag in keys:
        return info if info not in {".", ""} else flag
    if info in {".", ""}:
        return flag
    return info + ";" + flag


def rescue(
    in_vcf: Path,
    out_vcf: Path,
    *,
    clinvar_vcf: Path | None = None,
    gba1_browser: Path | None = None,
) -> int:
    alleles = load_clinvar(clinvar_vcf)
    gba1 = load_gba1_positions(gba1_browser)
    rescued = 0
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
    with _open_text(out_vcf, "wt") as dest:
        for line in _ensure_headers(header_lines):
            dest.write(line + "\n")
        for line in body_lines:
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 8:
                dest.write(line if line.endswith("\n") else line + "\n")
                continue
            filt = cols[6]
            if filt in {"PASS", "."}:
                dest.write("\t".join(cols) + "\n")
                continue
            chrom, pos_text, ref, alts = cols[0], cols[1], cols[3], cols[4]
            try:
                pos = int(pos_text)
            except ValueError:
                dest.write("\t".join(cols) + "\n")
                continue
            chrom_key = bare_chrom(chrom)
            fmt = cols[8] if len(cols) > 8 else "GT"
            if not _sample_is_carrier(cols, fmt):
                dest.write("\t".join(cols) + "\n")
                continue
            hit = (chrom_key, pos) in gba1
            if not hit:
                for alt in alts.split(","):
                    if (chrom_key, pos, ref.upper(), alt.upper()) in alleles:
                        hit = True
                        break
            if not hit:
                dest.write("\t".join(cols) + "\n")
                continue
            cols[6] = "PASS"
            cols[7] = _set_flag(cols[7], "PDVAR_RESCUE")
            dest.write("\t".join(cols) + "\n")
            rescued += 1
    return rescued


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--in-vcf", required=True)
    parser.add_argument("--out-vcf", required=True)
    parser.add_argument("--clinvar-vcf", default="")
    parser.add_argument("--gba1-browser", default="")
    args = parser.parse_args()
    clinvar = Path(args.clinvar_vcf) if args.clinvar_vcf.strip() else None
    browser = Path(args.gba1_browser) if args.gba1_browser.strip() else None
    count = rescue(Path(args.in_vcf), Path(args.out_vcf), clinvar_vcf=clinvar, gba1_browser=browser)
    print(f"rescue_nonpass_sites: rescued {count} site(s)", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
