"""Cross-class second hits for recessive panel genes.

SNV, indel, SV, and splice rows are already in one table. This step:

- marks an apparent homozygous SNV inside a heterozygous deletion as hemizygous
- uses phased genotypes (``|`` plus a shared ``PS``) when they are present
- leaves every other pair ``unphased``
- lists recessive genes that still have only one heterozygous hit

A shared ``PS`` tag on phased genotypes (``|``) is preferred over the unphased
compound-het label: disjoint alt haplotypes become ``AR_comphet_phased``
(``phase_source=PS``) and the same haplotype becomes ``AR_single_het``.
Statistical phasing runs only when a tool is on ``PATH`` and a reference panel
is configured; otherwise the status explains why it did not run.
"""

from __future__ import annotations

import shutil
from pathlib import Path

import pandas as pd

from .vcf_stream import chrom_bare, iter_vcf

SINGLE_HET_REASON = "single heterozygous hit in recessive gene, second allele not found"
PHASE_TOOLS = ("shapeit5", "shapeit4", "shapeit", "eagle", "beagle")
SOLVED_AR = {
    "AR_hom",
    "AR_comphet_candidate",
    "AR_comphet_phased",
    "AR_hemizygous",
    "AR_second_hit_splice",
}
FOLLOWUP_COLUMNS = [
    "sample_id", "gene", "moi", "variant_ids", "variant_classes", "zygosities",
    "phase_relation", "inheritance_status", "reason",
]


def _text(value: object) -> str:
    if value is None or (isinstance(value, float) and pd.isna(value)):
        return ""
    return str(value).strip()


def _ps_equal(left: object, right: object) -> bool:
    a, b = _text(left), _text(right)
    if a in {"", ".", "nan", "None", "<NA>"} or b in {"", ".", "nan", "None", "<NA>"}:
        return False
    return a == b


def phase_relation(gt1: object, ps1: object, gt2: object, ps2: object) -> str:
    """``trans`` / ``cis`` when both genotypes are phased in the same PS; otherwise ``unphased``."""

    def split(gt: object) -> tuple[list[str], bool]:
        text = _text(gt)
        if "|" in text:
            return text.split("|"), True
        return text.replace("|", "/").split("/"), False

    a1, phased1 = split(gt1)
    a2, phased2 = split(gt2)
    if not (phased1 and phased2):
        return "unphased"

    def alt_haps(alleles: list[str]) -> set[int]:
        return {i for i, allele in enumerate(alleles) if allele not in {"0", ".", ""}}

    h1, h2 = alt_haps(a1), alt_haps(a2)
    if not h1 or not h2 or not _ps_equal(ps1, ps2):
        return "unphased"
    if h1.isdisjoint(h2):
        return "trans"
    return "cis"


def statistical_phasing_status(enabled: bool = False, reference_panel: str | None = None) -> dict:
    """Report whether cohort statistical phasing can run. Does not invent phase."""
    found = [tool for tool in PHASE_TOOLS if shutil.which(tool)]
    if not enabled:
        return {"ran": False, "tool": found[0] if found else None, "reason": "disabled"}
    if not found:
        return {"ran": False, "tool": None, "reason": "tool_not_installed"}
    if not reference_panel:
        return {"ran": False, "tool": found[0], "reason": "reference_panel_not_configured"}
    return {"ran": False, "tool": found[0], "reason": "reference_panel_configured_not_executed", "reference_panel": reference_panel}


def statistical_phasing_command(tool: str, vcf: str, reference: str, region: str) -> list[str]:
    return [tool, "--input", vcf, "--reference", reference, "--region", region]


def _is_deletion(row: pd.Series) -> bool:
    if _text(row.get("variant_class")) != "SV":
        return False
    if _text(row.get("svtype")) == "DEL" or _text(row.get("dosage_event")) == "deletion":
        return True
    cn, ecn = row.get("copy_number"), row.get("expected_cn")
    try:
        return float(cn) < float(ecn)
    except (TypeError, ValueError):
        return False


def _is_het_deletion(row: pd.Series) -> bool:
    if not _is_deletion(row):
        return False
    if _text(row.get("zygosity")) == "het":
        return True
    try:
        return float(row.get("copy_number")) == float(row.get("expected_cn")) - 1
    except (TypeError, ValueError):
        return False


def _is_hom_deletion(row: pd.Series) -> bool:
    if not _is_deletion(row):
        return False
    if _text(row.get("zygosity")) == "hom":
        return True
    try:
        return float(row.get("copy_number")) == 0
    except (TypeError, ValueError):
        return False


def _looks_hom_alt(row: pd.Series) -> bool:
    if _text(row.get("variant_class")) == "SV":
        return False
    gt = _text(row.get("genotype"))
    if gt in {"1/1", "1|1"}:
        return True
    return _text(row.get("zygosity")) == "hom"


def _interval(row: pd.Series) -> tuple[int, int] | None:
    try:
        start = int(row["pos"])
    except (TypeError, ValueError, KeyError):
        return None
    try:
        end = int(row["end"]) if "end" in row and pd.notna(row["end"]) else start
    except (TypeError, ValueError):
        end = start
    return start, max(start, end)


def _overlaps(snv: pd.Series, sv: pd.Series) -> bool:
    snv_iv = _interval(snv)
    sv_iv = _interval(sv)
    if snv_iv is None or sv_iv is None:
        return False
    pos = snv_iv[0]
    return sv_iv[0] <= pos <= sv_iv[1] and chrom_bare(_text(snv.get("chrom"))) == chrom_bare(_text(sv.get("chrom")))


def _empty_followup() -> pd.DataFrame:
    return pd.DataFrame(columns=FOLLOWUP_COLUMNS)


def refine(df: pd.DataFrame, cfg: dict | None = None) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Update AR statuses for hemizygosity and phase, and build the follow-up list."""
    cfg = cfg or {}
    phase_cfg = (cfg.get("statistical_phasing") or {}) if isinstance(cfg, dict) else {}
    phasing = statistical_phasing_status(
        enabled=bool(phase_cfg.get("enabled", False)),
        reference_panel=phase_cfg.get("reference_panel"),
    )
    if df.empty:
        out = df.copy()
        out["phase_relation"] = pd.Series(dtype=object)
        out["phase_source"] = pd.Series(dtype=object)
        out["statistical_phasing"] = phasing["reason"]
        return out, _empty_followup()

    df = df.copy()
    df["phase_relation"] = "unphased"
    df["phase_source"] = ""
    df["statistical_phasing"] = phasing["reason"]
    if "moi" not in df.columns or "inheritance_status" not in df.columns:
        return df, _empty_followup()

    ar = df["moi"].astype(str) == "AR"
    for (_, _), grp in df[ar].groupby(["sample_id", "gene"], sort=False):
        idx = list(grp.index)
        _apply_hemizygosity(df, idx)
        _apply_phase(df, idx)

    follow_rows = []
    for (sample, gene), grp in df[ar].groupby(["sample_id", "gene"], sort=False):
        statuses = set(grp["inheritance_status"].astype(str))
        if statuses & SOLVED_AR:
            continue
        if "AR_single_het" not in statuses and "review_conflict" not in statuses:
            continue
        reason = (
            "homozygous SNV overlaps a homozygous deletion; review"
            if "review_conflict" in statuses
            else SINGLE_HET_REASON
        )
        follow_rows.append(
            {
                "sample_id": sample,
                "gene": gene,
                "moi": "AR",
                "variant_ids": ";".join(grp["variant_id"].astype(str)) if "variant_id" in grp.columns else "",
                "variant_classes": ";".join(grp.get("variant_class", pd.Series(dtype=str)).astype(str)),
                "zygosities": ";".join(grp["zygosity"].astype(str)) if "zygosity" in grp.columns else "",
                "phase_relation": ";".join(sorted(set(grp["phase_relation"].astype(str)))),
                "inheritance_status": ";".join(sorted(statuses)),
                "reason": reason,
            }
        )
    follow = pd.DataFrame(follow_rows) if follow_rows else _empty_followup()
    return df, follow


def _apply_hemizygosity(df: pd.DataFrame, idx: list) -> None:
    rows = df.loc[idx]
    snvs = [i for i in idx if _text(df.at[i, "variant_class"]) != "SV"] if "variant_class" in df.columns else list(idx)
    dels = [i for i in idx if _is_deletion(df.loc[i])]
    for snv_i in snvs:
        snv = df.loc[snv_i]
        if not _looks_hom_alt(snv):
            continue
        for del_i in dels:
            deletion = df.loc[del_i]
            if not _overlaps(snv, deletion):
                continue
            if _is_het_deletion(deletion):
                df.at[snv_i, "zygosity"] = "hemi"
                df.at[snv_i, "inheritance_status"] = "AR_hemizygous"
                df.at[del_i, "inheritance_status"] = "AR_hemizygous"
                df.at[snv_i, "phase_relation"] = "hemizygous_overlap"
                df.at[del_i, "phase_relation"] = "hemizygous_overlap"
            elif _is_hom_deletion(deletion):
                df.at[snv_i, "inheritance_status"] = "review_conflict"
                df.at[del_i, "inheritance_status"] = "review_conflict"
                df.at[snv_i, "phase_relation"] = "conflict_homozygous_deletion"
                df.at[del_i, "phase_relation"] = "conflict_homozygous_deletion"


def _qualifying(df: pd.DataFrame, index: int) -> bool:
    return _text(df.at[index, "inheritance_status"]) not in {"", "not_qualifying", "nan"}


def _apply_phase(df: pd.DataFrame, idx: list) -> None:
    hets = []
    for i in idx:
        if _text(df.at[i, "phase_relation"]) in {"hemizygous_overlap", "conflict_homozygous_deletion"}:
            continue
        if _text(df.at[i, "zygosity"]) != "het":
            continue
        if not _qualifying(df, i):
            continue
        hets.append(i)
    if len(hets) < 2:
        return
    relations = []
    for left in range(len(hets)):
        for right in range(left + 1, len(hets)):
            a, b = hets[left], hets[right]
            relations.append(
                phase_relation(
                    df.at[a, "genotype"] if "genotype" in df.columns else "",
                    df.at[a, "phase_set"] if "phase_set" in df.columns else "",
                    df.at[b, "genotype"] if "genotype" in df.columns else "",
                    df.at[b, "phase_set"] if "phase_set" in df.columns else "",
                )
            )
    if any(rel == "trans" for rel in relations):
        for i in hets:
            df.at[i, "inheritance_status"] = "AR_comphet_phased"
            df.at[i, "phase_relation"] = "trans"
            df.at[i, "phase_source"] = "PS"
        return
    if relations and all(rel == "cis" for rel in relations):
        for i in hets:
            df.at[i, "inheritance_status"] = "AR_single_het"
            df.at[i, "phase_relation"] = "cis"
            df.at[i, "phase_source"] = "PS"
        return
    for i in hets:
        if _text(df.at[i, "phase_relation"]) == "unphased":
            df.at[i, "phase_relation"] = "unphased"


def attach_phase_from_vcf(df: pd.DataFrame, vcf: str | Path) -> pd.DataFrame:
    """Replace genotypes at candidate positions with GT/PS from a region query.

    Only rows that already exist are updated. The joint VCF is not scanned.
    """
    if df.empty or not {"chrom", "pos", "sample_id", "ref", "alt"}.issubset(df.columns):
        return df
    vcf = Path(vcf)
    positions = (
        df[["chrom", "pos"]]
        .dropna()
        .drop_duplicates()
    )
    if positions.empty:
        return df
    import tempfile

    with tempfile.TemporaryDirectory() as tmp:
        bed = Path(tmp) / "sites.bed"
        with bed.open("w", encoding="utf-8") as fh:
            for chrom, pos in positions.itertuples(index=False):
                try:
                    coordinate = int(pos)
                except (TypeError, ValueError):
                    continue
                bare = chrom_bare(str(chrom))
                fh.write(f"chr{bare}\t{coordinate - 1}\t{coordinate}\n")
        wanted = {
            (str(row.sample_id), chrom_bare(str(row.chrom)), int(row.pos), str(row.ref), str(row.alt))
            for row in df.itertuples(index=False)
            if pd.notna(row.pos)
        }
        updates: dict[tuple, dict] = {}
        for record in iter_vcf(vcf, regions=bed):
            for sample, fmt in record["samples"].items():
                for alt in str(record["alt"]).split(","):
                    key = (str(sample), chrom_bare(record["chrom"]), int(record["pos"]), str(record["ref"]), alt)
                    if key in wanted:
                        updates[key] = {"genotype": fmt.get("GT", "./."), "phase_set": fmt.get("PS", "")}
    if not updates:
        return df
    df = df.copy()
    if "phase_set" not in df.columns:
        df["phase_set"] = pd.NA
    for idx, row in df.iterrows():
        try:
            key = (str(row["sample_id"]), chrom_bare(str(row["chrom"])), int(row["pos"]), str(row["ref"]), str(row["alt"]))
        except (TypeError, ValueError):
            continue
        if key in updates:
            df.at[idx, "genotype"] = updates[key]["genotype"]
            df.at[idx, "phase_set"] = updates[key]["phase_set"]
    return df
