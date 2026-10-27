"""Parse a single-sample VEP VCF into the columns the scorer expects.

VEP ``--pick`` CSQ supplies consequence, symbol, and HGVS. Optional ``--custom``
INFO fields from gnomAD genomes (``gnomADg_*``), gnomAD exomes (``gnomADe_*``),
legacy ``gnomAD`` / ``gnomAD_AF``, and ClinVar (``ClinVar_CLNSIG`` or ``ClinVar``)
override CSQ frequency and clinical significance when present. Unified
``gnomAD_AF`` / ``AF_popmax`` from ``LiftGnomadInfo`` are preferred when set.
CADD, REVEL, SpliceAI, and dbNSFP (Eigen / MutationTaster / fathmm-XF, plus CADD_phred
and REVEL_score from academic dbNSFP 5.4a) are read from CSQ when those plugins ran.
Cache-native ``--sift b --polyphen b`` fill SIFT and PolyPhen when enabled.
"""

from __future__ import annotations

import re
from pathlib import Path

import pandas as pd

from .normalize import add_variant_id, derive_zygosity, normalize_clinvar, variant_class_of
from .vcf_stream import _open_text, iter_vcf

_CSQ_FORMAT = re.compile(r"Format:\s*([^\"\s]+)")
_PRED_SCORE = re.compile(r"^(?P<pred>[A-Za-z_]+)\((?P<score>[-+0-9.eE]+)\)$")
_SPLICEAI = ("SpliceAI_pred_DS_AG", "SpliceAI_pred_DS_AL", "SpliceAI_pred_DS_DG", "SpliceAI_pred_DS_DL")
_EIGEN_KEYS = ("Eigen-phred_coding", "Eigen-phred", "Eigen-PC-phred_coding", "Eigen-raw_coding", "Eigen-raw")
_SIFT_SCORE_KEYS = ("SIFT_score",)
_SIFT_PRED_KEYS = ("SIFT_pred",)
_POLYPHEN_SCORE_KEYS = ("Polyphen2_HDIV_score", "PolyPhen2_HDIV_score", "Polyphen2_HVAR_score")
_POLYPHEN_PRED_KEYS = ("Polyphen2_HDIV_pred", "PolyPhen2_HDIV_pred", "Polyphen2_HVAR_pred")
_MUTTASTER_SCORE_KEYS = ("MutationTaster_score",)
_MUTTASTER_PRED_KEYS = ("MutationTaster_pred",)
# Academic dbNSFP 5.4a uses fathmm-XF_*; older builds used classic FATHMM_*.
_FATHMM_SCORE_KEYS = ("fathmm-XF_coding_score", "fathmm-XF_coding_rankscore", "FATHMM_score")
_FATHMM_PRED_KEYS = ("fathmm-XF_coding_pred", "FATHMM_pred")
# Plugin CADD/REVEL use CADD_PHRED / REVEL; dbNSFP 5.4a uses CADD_phred / REVEL_score.
_CADD_KEYS = ("CADD_PHRED", "CADD_phred", "CADD_RAW", "CADD_raw")
_REVEL_KEYS = ("REVEL", "REVEL_score")


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


def _max_float(values: list[float | None]) -> float | None:
    present = [value for value in values if value is not None]
    return max(present) if present else None


def _gnomad_from_info(info: dict) -> tuple[float | None, float | None]:
    """Prefer lifted ``gnomAD_AF`` / ``AF_popmax``, else max across genomes and exomes.

    Accepts VEP ``--custom`` short names ``gnomADg`` / ``gnomADe`` (separate INFO
    keys), legacy pipe ``gnomAD=AF|AF_popmax|nhomalt``, and ``LiftGnomadInfo``
    unified keys.
    """
    raw = info.get("gnomAD")
    legacy_af = legacy_pop = None
    if raw and raw is not True:
        parts = _pipe_fields(raw)
        if parts:
            legacy_af = _float(parts[0])
        if len(parts) > 1:
            legacy_pop = _float(parts[1])

    af = _max_float(
        [
            _float(info.get("gnomAD_AF")),
            _float(info.get("gnomADg_AF")),
            _float(info.get("gnomADe_AF")),
            legacy_af,
        ]
    )
    pop = _max_float(
        [
            _float(info.get("AF_popmax")),
            _float(info.get("gnomAD_AF_popmax")),
            _float(info.get("gnomADg_AF_popmax")),
            _float(info.get("gnomADe_AF_popmax")),
            legacy_pop,
        ]
    )
    return af, pop


def _clinvar_from_info(info: dict) -> str:
    """CLNSIG from VEP ``--custom`` (``ClinVar_CLNSIG``) or legacy ``ClinVar`` pipe."""
    for key in ("ClinVar_CLNSIG", "ClinVar"):
        raw = info.get(key)
        if raw and raw is not True:
            return str(raw).split("|")[0]
    return ""


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


def _first_token(value: object) -> str:
    """First non-empty / non-dot token from VEP or dbNSFP multi-value fields."""
    if value is None or value is True:
        return ""
    for part in re.split(r"[&;,/]", str(value)):
        text = part.strip()
        if text and text not in {".", "NA", "na", "None"}:
            return text
    return ""


def _pred_and_score(value: object) -> tuple[str | None, float | None]:
    """Parse VEP ``prediction(score)`` strings such as ``deleterious(0.01)``."""
    text = _first_token(value)
    if not text:
        return None, None
    match = _PRED_SCORE.match(text)
    if match:
        return match.group("pred").lower(), _float(match.group("score"))
    # Prediction-only token (tolerated, probably_damaging, D, ...).
    if re.fullmatch(r"[A-Za-z_]+", text):
        return text.lower(), None
    return None, _float(text)


def _csq_first(mapped: dict[str, str], names: tuple[str, ...]) -> str:
    for name in names:
        token = _first_token(mapped.get(name))
        if token:
            return token
    return ""


def _csq_float(mapped: dict[str, str], names: tuple[str, ...]) -> float | None:
    for name in names:
        number = _float(_first_token(mapped.get(name)))
        if number is not None:
            return number
    return None


def _in_silico_from_csq(csq: dict[str, str]) -> dict[str, object]:
    """SIFT/PolyPhen (cache or dbNSFP) plus Eigen / MutationTaster / fathmm-XF."""
    sift_pred, sift_score = _pred_and_score(csq.get("SIFT"))
    if sift_score is None:
        sift_score = _csq_float(csq, _SIFT_SCORE_KEYS)
    if not sift_pred:
        sift_pred = _csq_first(csq, _SIFT_PRED_KEYS).lower() or None

    poly_pred, poly_score = _pred_and_score(csq.get("PolyPhen"))
    if poly_score is None:
        poly_score = _csq_float(csq, _POLYPHEN_SCORE_KEYS)
    if not poly_pred:
        poly_pred = _csq_first(csq, _POLYPHEN_PRED_KEYS).lower() or None

    mut_pred = _csq_first(csq, _MUTTASTER_PRED_KEYS).upper() or None
    # fathmm-XF_coding_pred is D/N; classic FATHMM_pred is also D/N (uppercase).
    fathmm_pred = _csq_first(csq, _FATHMM_PRED_KEYS).upper() or None
    if fathmm_pred and fathmm_pred.startswith("D"):
        fathmm_pred = "D"
    elif fathmm_pred and fathmm_pred.startswith("N"):
        fathmm_pred = "N"
    return {
        "sift": sift_score,
        "sift_pred": sift_pred,
        "polyphen": poly_score,
        "polyphen_pred": poly_pred,
        "eigen": _csq_float(csq, _EIGEN_KEYS),
        "muttaster": _csq_float(csq, _MUTTASTER_SCORE_KEYS),
        "muttaster_pred": mut_pred,
        "fathmm": _csq_float(csq, _FATHMM_SCORE_KEYS),
        "fathmm_pred": fathmm_pred,
    }


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
            cadd = _csq_float(csq, _CADD_KEYS)
            revel = _csq_float(csq, _REVEL_KEYS)
            spliceai = _max_present(csq, _SPLICEAI)
            in_silico = _in_silico_from_csq(csq)
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
                        **in_silico,
                    }
                )
    frame = pd.DataFrame(rows)
    if frame.empty:
        return frame
    frame = add_variant_id(frame)
    frame = derive_zygosity(frame)
    frame["clinvar_class"] = normalize_clinvar(frame["clinvar"])
    return frame
