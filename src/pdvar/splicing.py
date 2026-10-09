"""Merge Introme splice predictions and fold them into the evidence score."""

from __future__ import annotations

import pandas as pd

from .scoring import assign_tier


def classify_introme(score: pd.Series, spcfg: dict) -> pd.Series:
    cls = pd.Series("none", index=score.index, dtype="object")
    cls[score >= spcfg["introme_moderate"]] = "moderate"
    cls[score >= spcfg["introme_high"]] = "high"
    return cls


def merge_introme(variants: pd.DataFrame, introme: pd.DataFrame) -> pd.DataFrame:
    cols = [c for c in introme.columns if c != "gene"]
    return variants.merge(introme[cols], on="variant_id", how="left")


def _lookup_consequence(region: pd.Series, intronic: pd.Series) -> pd.Series:
    region = region.fillna(".").astype(str)
    fallback = intronic.fillna("").astype(str).map({"intronic": "intron_variant", "exonic": "exonic_introme"})
    return region.where(region != ".", fallback.fillna("introme_variant"))


def find_second_hits(
    df: pd.DataFrame,
    carriers: pd.DataFrame,
    introme: pd.DataFrame,
    panel: pd.DataFrame,
    cfg: dict,
) -> pd.DataFrame:
    """Look up Introme splice variants as second alleles in the same sample.

    For every sample x AR-gene pair that already has a qualifying first hit,
    return that sample's Introme variants in the same gene that are not already
    in ``df`` and pass ``lookup_min_score``. When ``gnomad_af`` is present,
    alleles above ``lookup_max_gnomad_af`` are dropped. A missing AF is kept.
    There is no cohort allele frequency.
    Returned rows have ``hit_source = "introme_lookup"`` and only count as
    second hits.
    """
    spcfg, scfg, icfg = cfg["splicing"], cfg["scoring"], cfg["inheritance"]
    ar_genes = set(panel.loc[panel["moi"] == "AR", "gene"])
    is_first = df["gene"].isin(ar_genes) & (df["tier"] <= icfg["ar_first_hit_min_tier"])
    if "hit_source" in df.columns:
        is_first &= df["hit_source"] != "introme_lookup"
    pairs = df.loc[is_first, ["sample_id", "gene"]].astype(str).drop_duplicates()
    if pairs.empty or carriers.empty:
        return pd.DataFrame()

    cand = carriers.merge(pairs, on=["sample_id", "gene"])
    present = df[["sample_id", "variant_id"]].astype(str).drop_duplicates()
    cand = cand.merge(present, on=["sample_id", "variant_id"], how="left", indicator=True)
    cand = cand[cand["_merge"] == "left_only"].drop(columns="_merge")
    cand = cand.merge(introme, on="variant_id", how="left")
    max_af = spcfg.get("lookup_max_gnomad_af", 0.01)
    if "gnomad_af" in cand.columns and max_af is not None:
        af = pd.to_numeric(cand["gnomad_af"], errors="coerce")
        cand = cand[af.isna() | (af <= float(max_af))]
    cand = cand[cand["introme_score"] >= spcfg["lookup_min_score"]].reset_index(drop=True)
    if cand.empty:
        return pd.DataFrame()

    parts = cand["variant_id"].str.split("-", n=3, expand=True)
    cand["chrom"], cand["ref"], cand["alt"] = parts[0], parts[2], parts[3]
    cand["pos"] = parts[1].astype("Int64")
    cand["consequence"] = _lookup_consequence(
        cand.get("introme_gene_region", pd.Series(index=cand.index, dtype=object)),
        cand.get("introme_intronic", pd.Series(index=cand.index, dtype=object)),
    )
    cand["introme_class"] = classify_introme(cand["introme_score"], spcfg)
    high = cand["introme_class"] == "high"
    cand["evidence_score"] = high.astype(int) * scfg["points"]["introme_high"]
    cand["evidence_reasons"] = high.map({True: "introme_lookup;introme_high", False: "introme_lookup"})
    cand["tier"] = assign_tier(cand["evidence_score"], scfg["tiers"])
    cand["splice_candidate"] = True
    cand["hit_source"] = "introme_lookup"
    return cand.drop(columns=["allele_count"])


def apply_introme_evidence(df: pd.DataFrame, spcfg: dict, scfg: dict) -> pd.DataFrame:
    """Add ``introme_class`` and add Introme points to the score for relevant consequences.

    Introme points are only added when SpliceAI has not already contributed, to avoid
    double-counting the same splice evidence.
    """
    df = df.copy()
    score = df["introme_score"] if "introme_score" in df.columns else pd.Series(float("nan"), index=df.index)
    df["introme_class"] = classify_introme(score, spcfg)

    consider = set(spcfg["consider_consequences"])
    terms = df["consequence"].fillna("").astype(str).str.split(r"[&,]")
    relevant = terms.map(lambda ts: any(t.strip() in consider for t in ts)).astype(bool)
    already_splice = df["evidence_reasons"].fillna("").str.contains("spliceai").astype(bool)
    add = (df["introme_class"] == "high") & relevant & ~already_splice

    df["splice_candidate"] = (df["introme_class"] != "none") | already_splice
    df.loc[add, "evidence_score"] += scfg["points"]["introme_high"]
    df.loc[add, "evidence_reasons"] = (df.loc[add, "evidence_reasons"].fillna("") + ";introme_high").str.strip(";")
    df["tier"] = assign_tier(df["evidence_score"], scfg["tiers"])
    return df


def ensure_hit_source(df: pd.DataFrame, default: str = "small_variant") -> pd.DataFrame:
    """Fill a missing ``hit_source`` without overwriting SV or callset rows."""
    df = df.copy()
    if "hit_source" not in df.columns:
        df["hit_source"] = default
        return df
    blank = df["hit_source"].isna() | df["hit_source"].astype(str).str.strip().isin(["", "nan", "None", "<NA>"])
    df.loc[blank, "hit_source"] = default
    return df


def merge_splicevault(variants: pd.DataFrame, table: pd.DataFrame) -> pd.DataFrame:
    """Attach the highest-proportion SpliceVault outcome per variant. No extra points."""
    if table.empty or variants.empty:
        return variants
    table = table.copy()
    for column in ("chrom", "pos", "ref", "alt"):
        if column not in table.columns:
            return variants
    table["chrom"] = table["chrom"].astype(str).str.replace(r"^chr", "", regex=True)
    table["variant_id"] = (
        table["chrom"].astype(str)
        + "-"
        + table["pos"].astype(int).astype(str)
        + "-"
        + table["ref"].astype(str)
        + "-"
        + table["alt"].astype(str)
    )
    if "proportion" in table.columns:
        table = table.sort_values("proportion", ascending=False)
    table = table.drop_duplicates("variant_id")
    rename = {
        "outcome": "splicevault_outcome",
        "event": "splicevault_event",
        "proportion": "splicevault_proportion",
        "n_reads": "splicevault_n_reads",
        "transcript_id": "splicevault_transcript",
        "citation": "splicevault_citation",
        "exon_chrom": "splicevault_exon_chrom",
        "exon_start": "splicevault_exon_start",
        "exon_end": "splicevault_exon_end",
    }
    keep = ["variant_id", *[src for src in rename if src in table.columns]]
    extra = table[keep].rename(columns=rename)
    return variants.merge(extra, on="variant_id", how="left")


def summarize_splice_scores(df: pd.DataFrame) -> pd.DataFrame:
    """Record Pangolin and AlphaMissense when the columns already exist. No new points."""
    df = df.copy()
    parts = []
    for column in ("introme_pangolin_gain", "introme_pangolin_loss"):
        if column in df.columns:
            parts.append(pd.to_numeric(df[column], errors="coerce"))
    if parts:
        df["pangolin_max"] = pd.concat(parts, axis=1).max(axis=1)
    else:
        df["pangolin_max"] = pd.NA
    if "alphamissense" not in df.columns:
        df["alphamissense"] = pd.NA
    return df
