"""Parse a single-sample VEP VCF into the columns the scorer expects.

VEP ``--pick`` CSQ supplies consequence, symbol, and HGVS. Optional ``--custom``
INFO fields ``gnomAD`` (AF|AF_popmax|nhomalt) and ``ClinVar`` (CLNSIG) override
CSQ frequency and clinical significance when they are present. CADD, REVEL, and
SpliceAI are read from CSQ only when a plugin added those fields.
"""

from __future__ import annotations

import re
from pathlib import Path

import pandas as pd

from .normalize import add_variant_id, derive_zygosity, normalize_clinvar, variant_class_of
from .vcf_stream import _open_text, iter_vcf

_CSQ_FORMAT = re.compile(r"Format:\s*([^\"\s]+)")
_SPLICEAI = ("SpliceAI_pred_DS_AG", "SpliceAI_pred_DS_AL", "SpliceAI_pred_DS_DG", "SpliceAI_pred_DS_DL")


def csq_fields(header_text: str) -> list[str]:
    """Field names from the CSQ header description."""
    for line in header_text.splitlines():
        if line.startswith("##INFO=<ID=CSQ,"):
            match = _CSQ_FORMAT.search(line)
            if match:
                return match.group(1).split("|")
    return []


def _header_text(path: Path) -> str:
    lines = []
    with _open_text(path) as handle:
        for line in handle:
            if not line.startswith("#"):
                break
            lines.append(line.rstrip("\n"))
    return "\n".join(lines)


def _float(value: object) -> float | None:
    if value is None or value is True:
        return None
    text = str(value).split("&")[0].split(",")[0].strip()
    if text in {"", ".", "nan", "None"}:
        return None
    if "=" in text:
        text = text.split("=", 1)[1]
    try:
        return float(text)
    except ValueError:
        return None


def _pipe_fields(value: object) -> list[str]:
    if value is None or value is True:
        return []
    return str(value).split("|")


def _gnomad_from_info(info: dict) -> tuple[float | None, float | None]:
    """Custom VEP annotation ``gnomAD=AF|AF_popmax|nhomalt``, or lifted INFO flags."""
    af = _float(info.get("gnomAD_AF"))
    pop = _float(info.get("AF_popmax") or info.get("gnomAD_AF_popmax"))
    raw = info.get("gnomAD")
    if raw and raw is not True:
        parts = _pipe_fields(raw)
        if af is None and parts:
            af = _float(parts[0])
        if pop is None and len(parts) > 1:
            pop = _float(parts[1])
    return af, pop


def _clinvar_from_info(info: dict) -> str:
    raw = info.get("ClinVar")
    if not raw or raw is True:
        return ""
    return str(raw).split("|")[0]


def _csq_map(fields: list[str], value: str) -> dict[str, str]:
    parts = value.split("|")
    return {name: (parts[i] if i < len(parts) else "") for i, name in enumerate(fields)}


def _pick_csq(info: dict, fields: list[str], alt: str) -> dict[str, str]:
    raw = info.get("CSQ")
    if not raw or raw is True or not fields:
        return {}
    chosen = None
    for item in str(raw).split(","):
        mapped = _csq_map(fields, item)
        allele = mapped.get("Allele") or ""
        if allele and allele not in {alt, alt.split(",")[0]}:
            continue
        if mapped.get("CANONICAL") == "YES" or mapped.get("PICK") == "1":
            return mapped
        chosen = chosen or mapped
    return chosen or {}


def _max_present(mapped: dict[str, str], names: tuple[str, ...]) -> float | None:
    values = [_float(mapped.get(name)) for name in names]
    present = [value for value in values if value is not None]
    return max(present) if present else None


def _allele_balance(fmt: dict) -> float | None:
    ad = fmt.get("AD")
    if not ad or ad == ".":
        return None
    parts = []
    for item in str(ad).split(","):
        try:
            parts.append(float(item))
        except ValueError:
            return None
    if len(parts) < 2 or sum(parts) <= 0:
        return None
    return parts[1] / sum(parts)


def _carries(gt: str, alt_index: int) -> bool:
    alleles = [part for part in gt.replace("|", "/").split("/") if part not in {"", "."}]
    return str(alt_index) in alleles


def parse_vep(path: str | Path, sample_id: str) -> pd.DataFrame:
    """One row per non-reference allele the sample carries."""
    path = Path(path)
    fields = csq_fields(_header_text(path))
    rows: list[dict] = []
    wanted = str(sample_id)
    for record in iter_vcf(path, allow_full_scan=True):
        if record["filter"] not in {"PASS", "."}:
            continue
        alts = [alt for alt in str(record["alt"]).split(",") if alt and alt != "."]
        samples = record["samples"]
        if wanted in samples:
            chosen = {wanted: samples[wanted]}
        elif len(samples) == 1:
            chosen = {wanted: next(iter(samples.values()))}
        else:
            continue
        info = record["info"]
        info_af, info_pop = _gnomad_from_info(info)
        info_clinvar = _clinvar_from_info(info)
        for index, alt in enumerate(alts, start=1):
            csq = _pick_csq(info, fields, alt)
            gnomad_af = info_af if info_af is not None else _max_present(csq, ("gnomAD_AF", "gnomADe_AF", "gnomADg_AF", "MAX_AF"))
            popmax = info_pop if info_pop is not None else _max_present(
                csq, ("gnomAD_AF_popmax", "gnomADe_AF_popmax", "gnomADg_AF_popmax", "MAX_AF_POPS")
            )
            clinvar_raw = info_clinvar or csq.get("CLIN_SIG") or ""
            cadd = _float(csq.get("CADD_PHRED") or csq.get("CADD_RAW"))
            revel = _float(csq.get("REVEL"))
            spliceai = _max_present(csq, _SPLICEAI)
            gene = (csq.get("SYMBOL") or "").split("&")[0].strip() or "."
            for sample, fmt in chosen.items():
                gt = str(fmt.get("GT") or "./.")
                if not _carries(gt, index):
                    continue
                rows.append(
                    {
                        "sample_id": sample,
                        "family_id": sample,
                        "chrom": record["chrom"],
                        "pos": record["pos"],
                        "ref": record["ref"],
                        "alt": alt,
                        "gene": gene,
                        "consequence": csq.get("Consequence") or "sequence_variant",
                        "hgvsc": csq.get("HGVSc") or "",
                        "hgvsp": csq.get("HGVSp") or "",
                        "genotype": gt,
                        "gq": _float(fmt.get("GQ")),
                        "depth": _float(fmt.get("DP")),
                        "allele_balance": _allele_balance(fmt),
                        "gnomad_af": gnomad_af,
                        "gnomad_popmax_af": popmax,
                        "cadd": cadd,
                        "revel": revel,
                        "spliceai": spliceai,
                        "clinvar": clinvar_raw,
                        "variant_class": variant_class_of(record["ref"], alt),
                        "hit_source": "small_variant",
                    }
                )
    frame = pd.DataFrame(rows)
    if frame.empty:
        return frame
    frame = add_variant_id(frame)
    frame = derive_zygosity(frame)
    frame["clinvar_class"] = normalize_clinvar(frame["clinvar"])
    return frame
