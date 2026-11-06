#!/usr/bin/env python3
"""Join MitImpact 4 CSV annotations onto a MitochondriaPipeline VCF.

Matches ``chrM`` / ``MT`` / ``M`` alleles by POS + REF + ALT (first ALT).
When ``--table`` is unset, copies the VCF through and writes status
``not_assessed``. With a table, writes ``called`` even if no alleles match.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import sys
from pathlib import Path

# MitImpact CSV column -> VCF INFO tag (namespaced MI_*).
_FIELD_MAP = (
    ("Gene_symbol", "MI_gene"),
    ("class", "MI_moltype"),
    ("ref_aa", "MI_aa_ref"),
    ("alt_aa", "MI_aa_alt"),
    ("aa_pos", "MI_aa_pos"),
    ("HGVS_genomic", "MI_HGVS"),
    ("APOGEE2", "MI_APOGEE2"),
    ("APOGEE2_score", "MI_APOGEE2_score"),
    ("Clinvar_CLNSIG", "MI_ClinVar_CLNSIG"),
    ("Clinvar_CLNDN", "MI_ClinVar_CLNDN"),
    ("MITOMAP_Disease_Status", "MI_MITOMAP_status"),
    ("MITOMAP_Disease_GenBank_Freq", "MI_MITOMAP_GB_freq"),
    ("MLC_score", "MI_MLC_score"),
    ("AF_hom_gnomad", "MI_gnomAD_AF_hom"),
    ("AF_het_gnomad", "MI_gnomAD_AF_het"),
    ("MitoTIP_Score", "MI_MitoTIP"),
    ("Mitotip_Quartile", "MI_MitoTIP_quartile"),
)

_INFO_HEADERS = [
    '##INFO=<ID=MI_gene,Number=1,Type=String,Description="MitImpact gene symbol">',
    '##INFO=<ID=MI_moltype,Number=1,Type=String,Description="MitImpact molecule class">',
    '##INFO=<ID=MI_aa_ref,Number=1,Type=String,Description="MitImpact reference amino acid">',
    '##INFO=<ID=MI_aa_alt,Number=1,Type=String,Description="MitImpact alternate amino acid">',
    '##INFO=<ID=MI_aa_pos,Number=1,Type=String,Description="MitImpact amino-acid position">',
    '##INFO=<ID=MI_HGVS,Number=1,Type=String,Description="MitImpact HGVS genomic">',
    '##INFO=<ID=MI_APOGEE2,Number=1,Type=String,Description="APOGEE2 class from MitImpact">',
    '##INFO=<ID=MI_APOGEE2_score,Number=1,Type=Float,Description="APOGEE2 score from MitImpact">',
    '##INFO=<ID=MI_ClinVar_CLNSIG,Number=1,Type=String,Description="ClinVar CLNSIG from MitImpact">',
    '##INFO=<ID=MI_ClinVar_CLNDN,Number=1,Type=String,Description="ClinVar disease name from MitImpact">',
    '##INFO=<ID=MI_MITOMAP_status,Number=1,Type=String,Description="MITOMAP disease status from MitImpact">',
    '##INFO=<ID=MI_MITOMAP_GB_freq,Number=1,Type=String,Description="MITOMAP disease GenBank frequency from MitImpact">',
    '##INFO=<ID=MI_MLC_score,Number=1,Type=Float,Description="MitImpact mitochondrial local constraint score">',
    '##INFO=<ID=MI_gnomAD_AF_hom,Number=1,Type=Float,Description="gnomAD v3.1 mtDNA AF_hom from MitImpact">',
    '##INFO=<ID=MI_gnomAD_AF_het,Number=1,Type=Float,Description="gnomAD v3.1 mtDNA AF_het from MitImpact">',
    '##INFO=<ID=MI_MitoTIP,Number=1,Type=Float,Description="MitoTIP score from MitImpact">',
    '##INFO=<ID=MI_MitoTIP_quartile,Number=1,Type=String,Description="MitoTIP quartile from MitImpact">',
]


def _open_text(path: Path, mode: str = "rt"):
    if str(path).endswith(".gz"):
        return gzip.open(path, mode, encoding="utf-8" if "t" in mode else None)
    return path.open(mode, encoding="utf-8" if "t" in mode else None)


def _is_chr_m(chrom: str) -> bool:
    text = str(chrom).strip()
    rest = text[3:] if text.lower().startswith("chr") else text
    return rest.upper() in {"M", "MT"}


def _escape_info(value: str) -> str:
    return (
        str(value)
        .replace("\\", "\\\\")
        .replace(";", "%3B")
        .replace("=", "%3D")
        .replace(" ", "_")
        .replace(",", "%2C")
        .replace("\t", "_")
    )


def load_mitimpact(path: Path) -> dict[tuple[int, str, str], dict[str, str]]:
    """Index MitImpact rows by (pos, ref, alt)."""
    table: dict[tuple[int, str, str], dict[str, str]] = {}
    with path.open("r", encoding="utf-8", errors="replace", newline="") as handle:
        reader = csv.DictReader(handle)
        required = {"Start", "Ref", "Alt"}
        if reader.fieldnames is None or not required.issubset(set(reader.fieldnames)):
            raise ValueError(
                f"MitImpact table {path} missing Start/Ref/Alt columns; got {reader.fieldnames}"
            )
        for row in reader:
            try:
                pos = int(str(row["Start"]).strip())
            except (TypeError, ValueError):
                continue
            ref = str(row.get("Ref") or "").strip().upper()
            alt = str(row.get("Alt") or "").strip().upper()
            if not ref or not alt:
                continue
            values = {}
            for src, _dest in _FIELD_MAP:
                raw = row.get(src)
                if raw is None:
                    continue
                text = str(raw).strip()
                if text and text.lower() not in {"na", "nan", "."}:
                    values[src] = text
            if values:
                table[(pos, ref, alt)] = values
    return table


def _info_blob(values: dict[str, str]) -> str:
    parts = []
    for src, dest in _FIELD_MAP:
        if src not in values:
            continue
        parts.append(f"{dest}={_escape_info(values[src])}")
    return ";".join(parts)


def annotate(vcf: Path, out_vcf: Path, table: Path | None, status_path: Path) -> None:
    if table is None:
        status_path.write_text("not_assessed\n", encoding="utf-8")
        with _open_text(vcf) as src, _open_text(out_vcf, "wt") as dest:
            for line in src:
                dest.write(line if line.endswith("\n") else line + "\n")
        return

    index = load_mitimpact(table)
    status_path.write_text("called\n", encoding="utf-8")
    matched = 0
    with _open_text(vcf) as src, _open_text(out_vcf, "wt") as dest:
        header_done = False
        for line in src:
            if line.startswith("##"):
                dest.write(line if line.endswith("\n") else line + "\n")
                continue
            if line.startswith("#CHROM"):
                for header in _INFO_HEADERS:
                    dest.write(header + "\n")
                dest.write(
                    '##MitImpactAnnotation=MitImpact4 CSV join by POS/REF/ALT; '
                    f"rows_indexed={len(index)}\n"
                )
                dest.write(line if line.endswith("\n") else line + "\n")
                header_done = True
                continue
            if not header_done:
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 8:
                dest.write(line if line.endswith("\n") else line + "\n")
                continue
            chrom, pos_text, _vid, ref, alts = parts[0], parts[1], parts[2], parts[3], parts[4]
            if not _is_chr_m(chrom):
                dest.write(line if line.endswith("\n") else line + "\n")
                continue
            try:
                pos = int(pos_text)
            except ValueError:
                dest.write(line if line.endswith("\n") else line + "\n")
                continue
            alt0 = alts.split(",")[0].strip().upper()
            ref_u = ref.strip().upper()
            values = index.get((pos, ref_u, alt0))
            if values:
                matched += 1
                blob = _info_blob(values)
                info = parts[7]
                parts[7] = blob if info in {".", ""} else f"{info};{blob}"
            dest.write("\t".join(parts) + "\n")
    print(f"mitimpact_indexed={len(index)} matched={matched}", file=sys.stderr)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--vcf", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--table", type=Path, default=None)
    parser.add_argument("--status", type=Path, required=True)
    args = parser.parse_args(argv)
    annotate(args.vcf, args.out, args.table, args.status)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
