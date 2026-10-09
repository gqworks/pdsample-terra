"""Column mapping and variant-key harmonisation."""

from __future__ import annotations

import ast
import re

import pandas as pd

NUMERIC_COLUMNS = (
    "pos", "depth", "gq", "allele_balance", "cohort_af", "gnomad_af", "gnomad_popmax_af",
    "gnomad_hom", "cadd", "revel", "spliceai", "clinvar_stars", "introme_score",
    "age_onset", "age_sampling",
)


def apply_column_map(df: pd.DataFrame, mapping: dict[str, str]) -> pd.DataFrame:
    """Rename source columns to canonical names; unmapped source columns are kept as-is."""
    rename = {src: canon for canon, src in mapping.items() if src in df.columns}
    out = df.rename(columns=rename)
    for col in NUMERIC_COLUMNS:
        if col in out.columns:
            out[col] = pd.to_numeric(out[col], errors="coerce")
    return out


def normalize_chrom(chrom: pd.Series) -> pd.Series:
    return chrom.astype(str).str.replace(r"^chr", "", regex=True).str.upper().replace({"MT": "M"})


def add_variant_id(df: pd.DataFrame) -> pd.DataFrame:
    """Add ``variant_id`` as chrom-pos-ref-alt (chrom without 'chr' prefix)."""
    df = df.copy()
    df["chrom"] = normalize_chrom(df["chrom"])
    df["pos"] = df["pos"].astype("Int64")
    df["variant_id"] = (
        df["chrom"] + "-" + df["pos"].astype(str) + "-" + df["ref"].astype(str) + "-" + df["alt"].astype(str)
    )
    return df


def derive_zygosity(df: pd.DataFrame) -> pd.DataFrame:
    """Fill ``zygosity`` (het / hom / hemi) from ``genotype`` when not supplied."""
    df = df.copy()
    if "zygosity" in df.columns and df["zygosity"].notna().any():
        df["zygosity"] = df["zygosity"].astype(str).str.lower()
        return df
    if "genotype" not in df.columns:
        df["zygosity"] = pd.NA
        return df

    def _zyg(gt: object) -> object:
        if pd.isna(gt):
            return pd.NA
        alleles = str(gt).replace("|", "/").split("/")
        if len(alleles) == 1:
            return "hemi" if alleles[0] not in ("0", ".") else "ref"
        alt = [a for a in alleles if a not in ("0", ".")]
        if not alt:
            return "ref"
        return "hom" if len(alt) == len(alleles) and len(set(alt)) == 1 else "het"

    df["zygosity"] = df["genotype"].map(_zyg)
    return df


def normalize_clinvar(sig: pd.Series) -> pd.Series:
    """Collapse ClinVar significance strings to P/LP, VUS, conflicting, B/LB or none."""
    s = sig.fillna("").astype(str).str.lower().str.replace("_", " ")

    def _cat(v: str) -> str:
        if not v:
            return "none"
        if "conflicting" in v:
            return "conflicting"
        if "pathogenic" in v and "benign" not in v:
            return "P/LP"
        if "benign" in v and "pathogenic" not in v:
            return "B/LB"
        if "uncertain" in v:
            return "VUS"
        return "other"

    return s.map(_cat)


def apply_gene_aliases(genes: pd.Series, aliases: dict[str, str] | None) -> pd.Series:
    """Upper-case gene symbols and map legacy symbols (e.g. GBA -> GBA1)."""
    genes = genes.astype(str).str.strip().str.upper()
    if aliases:
        genes = genes.replace({k.upper(): v.upper() for k, v in aliases.items()})
    return genes


def variant_class_of(ref: object, alt: object) -> str:
    """SNV, INDEL, or SV from allele strings."""
    r, a = str(ref), str(alt)
    if r.startswith("<") or a.startswith("<"):
        return "SV"
    if len(r) == len(a) == 1:
        return "SNV"
    return "INDEL"


INTROME_NUMERIC = (
    "introme_score", "introme_qual", "introme_pangolin_gain", "introme_pangolin_loss",
    "introme_mmsplice_pathogenicity", "introme_spip_score",
    "introme_spliceai_ds_ag", "introme_spliceai_ds_al", "introme_spliceai_ds_dg", "introme_spliceai_ds_dl",
)
INTROME_SPLICEAI_DS = INTROME_NUMERIC[-4:]
INTROME_ANNOTATIONS = (
    "introme_rsid", "introme_qual", "introme_score", "introme_gene_region", "introme_intronic",
    "introme_spliceai_max", "introme_pangolin_gain", "introme_pangolin_loss",
    "introme_mmsplice_pathogenicity", "introme_spip_score",
)


def _literal(value: object, default: object) -> object:
    """Parse Python-literal strings such as ``"('C',)"`` or ``"[('SAMPLE_1', 1)]"``."""
    if value is None or (isinstance(value, float) and pd.isna(value)) or str(value).strip() in ("", "."):
        return default
    try:
        return ast.literal_eval(str(value))
    except (ValueError, SyntaxError):
        return default


def _as_list(value: object) -> list:
    return list(value) if isinstance(value, (list, tuple)) else [value]


def _gene_symbols(value: object) -> list[str]:
    """Accept a legacy Python list, a bare symbol, or a comma/semicolon list.

    ``.`` and blank cells become an empty list. A raw Introme v2 TSV omits ``gene``.
    """
    if value is None or (isinstance(value, float) and pd.isna(value)):
        return []
    text = str(value).strip()
    if text in ("", ".", "nan", "None", "[]"):
        return []
    parsed = _literal(value, None)
    if isinstance(parsed, (list, tuple)):
        symbols = [str(item) for item in parsed]
    elif isinstance(parsed, str):
        symbols = re.split(r"[;,]", parsed)
    else:
        symbols = re.split(r"[;,]", text)
    cleaned = []
    for item in symbols:
        symbol = str(item).strip().strip("'\"")
        if symbol and symbol not in (".", "nan", "None"):
            cleaned.append(symbol)
    return cleaned


def _carrier_pairs(value: object) -> list:
    """Legacy ``[('SAMPLE_A', 1)]`` pairs. Missing or non-list cells yield no carriers."""
    parsed = _literal(value, [])
    if not isinstance(parsed, (list, tuple)):
        return []
    pairs = []
    for item in parsed:
        if isinstance(item, (list, tuple)) and len(item) >= 2:
            pairs.append(item)
    return pairs


def normalize_introme(
    df: pd.DataFrame, mapping: dict[str, str], aliases: dict[str, str] | None = None
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Parse the Introme export into a variant-level lookup and a per-sample carrier table.

    Returns ``(variants, carriers)``:
      variants: one row per ``variant_id`` with Introme scores, ``introme_genes`` and
                carrier counts (``introme_n_carriers``, ``introme_n_hom``, ``introme_cohort_ac``)
      carriers: one row per variant x gene x sample with ``allele_count`` and ``zygosity``
    Rows repeated for the same variant (split carrier lists) are merged.
    """
    df = apply_column_map(df, mapping)
    if "gene" not in df.columns and "INFO:gene" in df.columns:
        df["gene"] = df["INFO:gene"]
    df["alt"] = df["alt"].map(lambda v: str(_as_list(_literal(v, v))[0]))
    df = add_variant_id(df)
    for col in INTROME_NUMERIC:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")
    ds = [c for c in INTROME_SPLICEAI_DS if c in df.columns]
    if ds:
        df["introme_spliceai_max"] = df[ds].max(axis=1)

    if "gene" not in df.columns:
        df["genes"] = [[] for _ in range(len(df))]
    else:
        df["genes"] = df["gene"].map(
            lambda v: list(apply_gene_aliases(pd.Series(_gene_symbols(v), dtype=str), aliases))
        )
    if "carriers" not in df.columns:
        df["carriers"] = [[] for _ in range(len(df))]
    else:
        df["carriers"] = df["carriers"].map(_carrier_pairs)

    carriers = df[["variant_id", "genes", "carriers"]].explode("carriers").dropna(subset=["carriers"])
    carriers["sample_id"] = carriers["carriers"].map(lambda c: str(c[0]))
    carriers["allele_count"] = carriers["carriers"].map(lambda c: int(c[1]))
    carriers["zygosity"] = carriers["allele_count"].map({1: "het", 2: "hom"}).fillna("het")
    carriers = (
        carriers.explode("genes")
        .rename(columns={"genes": "gene"})
        [["variant_id", "gene", "sample_id", "allele_count", "zygosity"]]
        .drop_duplicates(["variant_id", "gene", "sample_id"])
        .reset_index(drop=True)
    )

    per_sample = carriers.drop_duplicates(["variant_id", "sample_id"])
    counts = per_sample.groupby("variant_id").agg(
        introme_n_carriers=("sample_id", "nunique"),
        introme_n_hom=("allele_count", lambda s: int((s == 2).sum())),
        introme_cohort_ac=("allele_count", "sum"),
    )
    keep = [c for c in INTROME_ANNOTATIONS if c in df.columns]
    variants = df.drop_duplicates("variant_id")[["variant_id", "genes", *keep]].copy()
    variants["introme_genes"] = variants.pop("genes").map(";".join)
    variants = variants.merge(counts, on="variant_id", how="left")
    for col in ("introme_n_carriers", "introme_n_hom", "introme_cohort_ac"):
        if col not in variants.columns:
            variants[col] = 0
        else:
            variants[col] = variants[col].fillna(0)
    return variants.reset_index(drop=True), carriers


def normalize_phenotype(df: pd.DataFrame, mapping: dict[str, str]) -> pd.DataFrame:
    df = apply_column_map(df, mapping)
    df["sample_id"] = df["sample_id"].astype(str)
    return df


def normalize_panel(df: pd.DataFrame, mapping: dict[str, str], aliases: dict[str, str] | None = None) -> pd.DataFrame:
    df = apply_column_map(df, mapping)
    df["gene"] = apply_gene_aliases(df["gene"], aliases)
    if "moi" in df.columns:
        df["moi"] = df["moi"].fillna("unknown").astype(str).str.strip().str.upper()
    if "risk_allele_exception" in df.columns:
        df["risk_allele_exception"] = (
            df["risk_allele_exception"].astype(str).str.strip().str.lower().isin(["true", "1", "yes"])
        )
    return df.drop_duplicates("gene")
