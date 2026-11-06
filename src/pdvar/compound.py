"""Cross-class second hits for recessive panel genes.

SNV, indel, SV, and splice rows are already in one table. This step:

- marks an apparent homozygous SNV inside a heterozygous deletion as hemizygous
- pairs a LoF / exon-disrupting SV with a coding (or ClinVar P/LP) allele as
  ``AR_comphet_candidate``, and demotes same-gene Introme-only intron "solves"
- lists recessive genes that still have only one heterozygous hit
"""

from __future__ import annotations

import pandas as pd

from .scoring import NULL_CONSEQUENCES
from .vcf_stream import chrom_bare

SINGLE_HET_REASON = "single heterozygous hit in recessive gene, second allele not found"
SOLVED_AR = {
    "AR_hom",
    "AR_comphet_candidate",
    "AR_hemizygous",
    "AR_second_hit_splice",
}
# Intronic Introme support kept visible but ranked below coding / SV solves.
INTRONIC_SUPPORT_STATUS = "AR_intronic_candidate"
CODING_CONSEQUENCES = NULL_CONSEQUENCES | {
    "missense_variant",
    "inframe_insertion",
    "inframe_deletion",
    "protein_altering_variant",
    "stop_lost",
    "splice_region_variant",
}
# Gene-level demotion: require true coding/null/P-LP/LoF-SV, not splice_region+intron noise.
STRONG_CODING_CONSEQUENCES = NULL_CONSEQUENCES | {
    "missense_variant",
    "inframe_insertion",
    "inframe_deletion",
    "protein_altering_variant",
    "stop_lost",
}
WEAK_NONCODING_TERMS = {
    "intron_variant",
    "non_coding_transcript_variant",
    "splice_region_variant",
    "splice_polypyrimidine_tract_variant",
    "donor_canonical",
    "acceptor_canonical",
    "donor_region",
    "acceptor_region",
    "3_prime_UTR_variant",
    "5_prime_UTR_variant",
}
FOLLOWUP_COLUMNS = [
    "sample_id", "gene", "moi", "variant_ids", "variant_classes", "zygosities",
    "inheritance_status", "reason",
]


def _text(value: object) -> str:
    if value is None or (isinstance(value, float) and pd.isna(value)):
        return ""
    return str(value).strip()


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


def _consequence_terms(row: pd.Series) -> set[str]:
    return {t.strip() for t in _text(row.get("consequence")).replace(",", "&").split("&") if t.strip()}


def _is_lof_sv(row: pd.Series) -> bool:
    """Exon-disrupting / LoF SV that can act as an AR second hit (DEL or DUP)."""
    if _text(row.get("variant_class")) != "SV" and _text(row.get("hit_source")) != "sv":
        return False
    if _text(row.get("dosage_mechanism")).lower() == "lof":
        return True
    if _text(row.get("scored_dosage")) in {"True", "true", "1"} or row.get("scored_dosage") is True:
        if _text(row.get("svtype")) in {"DEL", "DUP", "CNV"} or _text(row.get("dosage_event")) in {
            "deletion",
            "duplication",
            "triplication_or_higher",
        }:
            return True
    return bool(_consequence_terms(row) & {"transcript_ablation"})


def _is_coding_support(row: pd.Series) -> bool:
    """Coding / ClinVar P/LP small-variant support (not intron-only noise)."""
    if _is_lof_sv(row):
        return False
    if _text(row.get("variant_class")) == "SV":
        return False
    if _text(row.get("clinvar_class")) == "P/LP":
        return True
    terms = _consequence_terms(row)
    if terms & CODING_CONSEQUENCES:
        return True
    return False


def _is_strong_allele(row: pd.Series) -> bool:
    """True coding/null, ClinVar P/LP, or LoF SV — not intron/splice noise."""
    if _is_lof_sv(row):
        return True
    if _text(row.get("clinvar_class")) == "P/LP":
        return True
    if _text(row.get("variant_class")) == "SV":
        return False
    terms = _consequence_terms(row)
    return bool(terms & STRONG_CODING_CONSEQUENCES)


def _is_weak_intronic(row: pd.Series) -> bool:
    """Deep-intronic / splice-region noise without strong coding evidence."""
    if _is_strong_allele(row):
        return False
    terms = _consequence_terms(row)
    if not terms:
        return False
    if terms & STRONG_CODING_CONSEQUENCES:
        return False
    if "intron_variant" in terms or "non_coding_transcript_variant" in terms:
        return True
    if terms <= WEAK_NONCODING_TERMS:
        return True
    if terms & {
        "donor_canonical",
        "acceptor_canonical",
        "donor_region",
        "acceptor_region",
    }:
        return True
    return False


def _locked_status(status: str) -> bool:
    return status in {"AR_hemizygous", "review_conflict"}


def _pair_sv_with_coding(df: pd.DataFrame, idx: list) -> None:
    """Promote LoF SV + coding allele to comp-het; demote same-gene weak intron solves."""
    rescue_statuses = {"AR_single_het", "not_qualifying", ""}
    lof_svs = [i for i in idx if _is_lof_sv(df.loc[i]) and _qualifying(df, i)]
    # Rescue LoF SVs left not_qualifying (tier floor) so they can still pair.
    lof_svs += [
        i
        for i in idx
        if i not in lof_svs
        and _is_lof_sv(df.loc[i])
        and _text(df.at[i, "inheritance_status"]) in rescue_statuses
    ]
    coding = [i for i in idx if _is_coding_support(df.loc[i]) and _qualifying(df, i)]
    # Also rescue coding rows left as AR_single_het / not yet paired.
    coding += [
        i
        for i in idx
        if i not in coding
        and _is_coding_support(df.loc[i])
        and _text(df.at[i, "inheritance_status"]) in rescue_statuses
    ]
    if not lof_svs or not coding:
        return
    partners = sorted(set(lof_svs + coding))
    for i in partners:
        status = _text(df.at[i, "inheritance_status"])
        if _locked_status(status):
            continue
        df.at[i, "inheritance_status"] = "AR_comphet_candidate"
    for i in idx:
        if i in partners:
            continue
        if not _is_weak_intronic(df.loc[i]):
            continue
        status = _text(df.at[i, "inheritance_status"])
        if status in SOLVED_AR or status == "AR_single_het":
            df.at[i, "inheritance_status"] = INTRONIC_SUPPORT_STATUS


def _demote_intron_only_ar(df: pd.DataFrame, idx: list) -> None:
    """Demote intron/Introme AR labels that should not outrank coding evidence.

    - No strong allele in the gene: demote intron solves / single-hets to support-only
      (avoid intron-only AR_hom "solves").
    - Strong coding/LoF allele present: demote weak intron rows labelled ``AR_hom``
      so an Introme homozygous intron cannot bury coding compound-het partners.
      If that demotion leaves a lone ``AR_comphet_candidate``, drop it to
      ``AR_single_het`` (the demoted hom was the phantom second hit).
    """
    has_strong = any(_is_strong_allele(df.loc[i]) for i in idx)
    demoted_hom = False
    for i in idx:
        status = _text(df.at[i, "inheritance_status"])
        if _locked_status(status):
            continue
        if not _is_weak_intronic(df.loc[i]):
            continue
        if has_strong:
            if status == "AR_hom":
                df.at[i, "inheritance_status"] = INTRONIC_SUPPORT_STATUS
                demoted_hom = True
            continue
        if status not in SOLVED_AR and status != "AR_single_het":
            continue
        df.at[i, "inheritance_status"] = INTRONIC_SUPPORT_STATUS
    if has_strong and demoted_hom:
        comphets = [i for i in idx if _text(df.at[i, "inheritance_status"]) == "AR_comphet_candidate"]
        other_solved = [
            i
            for i in idx
            if _text(df.at[i, "inheritance_status"]) in SOLVED_AR
            and _text(df.at[i, "inheritance_status"]) != "AR_comphet_candidate"
        ]
        if len(comphets) == 1 and not other_solved:
            df.at[comphets[0], "inheritance_status"] = "AR_single_het"


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


def refine(df: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Update AR statuses for hemizygosity and SV second hits, and build the follow-up list."""
    if df.empty or "moi" not in df.columns or "inheritance_status" not in df.columns:
        return df.copy(), _empty_followup()

    df = df.copy()
    ar = df["moi"].astype(str) == "AR"
    for (_, _), grp in df[ar].groupby(["sample_id", "gene"], sort=False):
        idx = list(grp.index)
        _apply_hemizygosity(df, idx)
        _pair_sv_with_coding(df, idx)
        _demote_intron_only_ar(df, idx)

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
            elif _is_hom_deletion(deletion):
                df.at[snv_i, "inheritance_status"] = "review_conflict"
                df.at[del_i, "inheritance_status"] = "review_conflict"


def _qualifying(df: pd.DataFrame, index: int) -> bool:
    return _text(df.at[index, "inheritance_status"]) not in {"", "not_qualifying", "nan"}
