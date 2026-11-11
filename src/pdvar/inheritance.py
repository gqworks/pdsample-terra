"""Mode-of-inheritance logic, including AR second-hit / compound-het detection."""

from __future__ import annotations

import pandas as pd


PANEL_ANNOTATIONS = {
    "moi": "moi",
    "moi_notes": "moi_notes",
    "confidence": "panel_confidence",
    "risk_allele_exception": "risk_allele_exception",
    "mechanism": "mechanism",
    "mechanism_notes": "mechanism_notes",
}

# AD genes where splice LoF does not match the Mendelian disease mechanism.
GOF_LIKE_MECHANISMS = frozenset({"gof", "dn"})

# Weak noncoding terms that Introme may elevate but that do not create a GoF/DN protein.
_WEAK_INTRONIC_TERMS = frozenset(
    {
        "intron_variant",
        "non_coding_transcript_variant",
        "splice_region_variant",
        "splice_polypyrimidine_tract_variant",
        "splice_donor_5th_base_variant",
        "donor_canonical",
        "acceptor_canonical",
        "donor_region",
        "acceptor_region",
    }
)
_CODINGISH_TERMS = frozenset(
    {
        "missense_variant",
        "frameshift_variant",
        "stop_gained",
        "stop_lost",
        "start_lost",
        "splice_acceptor_variant",
        "splice_donor_variant",
        "transcript_ablation",
        "inframe_insertion",
        "inframe_deletion",
        "protein_altering_variant",
        "structural_variant",
        "copy_number_variation",
        "copy_number_gain",
        "copy_number_loss",
    }
)

AD_GOF_INTRONIC_STATUS = "AD_gof_intronic"


def annotate_moi(df: pd.DataFrame, panel: pd.DataFrame) -> pd.DataFrame:
    """Add panel MOI, notes, confidence, mechanism and risk-allele flag; genes not on the panel get ``NOT_IN_PANEL``."""
    df = df.copy()
    genes = df["gene"].astype(str).str.upper()
    info = panel.set_index("gene")
    for src, dest in PANEL_ANNOTATIONS.items():
        if src in info.columns:
            df[dest] = genes.map(info[src])
    df["moi"] = df["moi"].fillna("NOT_IN_PANEL") if "moi" in df.columns else "NOT_IN_PANEL"
    if "mechanism" in df.columns:
        df["mechanism"] = df["mechanism"].fillna("").astype(str).str.strip().str.lower()
    return df


def _consequence_term_sets(series: pd.Series) -> list[set[str]]:
    out: list[set[str]] = []
    for raw in series.fillna("").astype(str):
        terms: set[str] = set()
        for part in raw.replace(",", "&").split("&"):
            t = part.strip()
            if t:
                terms.add(t)
        out.append(terms)
    return out


def clinvar_vus_mask(df: pd.DataFrame) -> pd.Series:
    """True where ``clinvar_class`` is ClinVar VUS (coding or non-coding)."""
    if "clinvar_class" not in df.columns:
        return pd.Series(False, index=df.index, dtype=bool)
    return df["clinvar_class"].fillna("").astype(str).eq("VUS")


def clinvar_vus_coding_mask(df: pd.DataFrame) -> pd.Series:
    """True for ClinVar VUS alleles with a codingish consequence (missense/null/inframe/…)."""
    vus = clinvar_vus_mask(df)
    if not vus.any():
        return vus
    cons = df["consequence"] if "consequence" in df.columns else pd.Series("", index=df.index)
    term_sets = _consequence_term_sets(cons)
    coding = pd.Series(
        [bool(terms & _CODINGISH_TERMS) for terms in term_sets],
        index=df.index,
        dtype=bool,
    )
    return vus & coding


def weak_intronic_mask(df: pd.DataFrame) -> pd.Series:
    """True for intron / splice-region noise without coding or ClinVar P/LP support."""
    clin = (
        df["clinvar_class"].fillna("").astype(str)
        if "clinvar_class" in df.columns
        else pd.Series("", index=df.index, dtype=str)
    )
    cons = df["consequence"] if "consequence" in df.columns else pd.Series("", index=df.index)
    term_sets = _consequence_term_sets(cons)
    weak = []
    for terms, cv in zip(term_sets, clin.tolist()):
        if cv == "P/LP":
            weak.append(False)
            continue
        if terms & _CODINGISH_TERMS:
            weak.append(False)
            continue
        weak.append(bool(terms & _WEAK_INTRONIC_TERMS) or (bool(terms) and terms <= _WEAK_INTRONIC_TERMS))
    return pd.Series(weak, index=df.index, dtype=bool)


def demote_gof_ad_intronic(df: pd.DataFrame, icfg: dict) -> pd.DataFrame:
    """Demote weak intron AD hits in GoF / DN genes to ``AD_gof_intronic``.

    Splice-disrupting alleles remove or truncate product; they do not create the
    toxic / hyperactive / dominant-negative proteins that drive these AD diseases.
    """
    if not icfg.get("demote_gof_ad_intronic", True):
        return df
    if "mechanism" not in df.columns or "inheritance_status" not in df.columns:
        return df
    mechanisms = {
        str(m).strip().lower()
        for m in (icfg.get("gof_intronic_mechanisms") or sorted(GOF_LIKE_MECHANISMS))
    } or set(GOF_LIKE_MECHANISMS)
    df = df.copy()
    gof_like = df["mechanism"].fillna("").astype(str).str.lower().isin(mechanisms)
    ad = df["inheritance_status"].fillna("").astype(str).eq("AD_candidate")
    demote = ad & gof_like & weak_intronic_mask(df)
    df.loc[demote, "inheritance_status"] = AD_GOF_INTRONIC_STATUS
    return df


def _qualifies(df: pd.DataFrame, min_tier: int) -> pd.Series:
    return df["tier"] <= min_tier


def evaluate_inheritance(df: pd.DataFrame, icfg: dict) -> pd.DataFrame:
    """Assign ``inheritance_status`` per variant, evaluated within each sample x gene.

    Statuses:
      AD_candidate          qualifying variant in an AD gene
      AD_gof_intronic       weak intron/Introme hit in a GoF/DN AD gene (demoted)
      AR_hom                qualifying homozygous variant in an AR gene
                            (other same-gene alleles are still labelled; a hom
                            does not leave coding partners as not_qualifying)
      AR_comphet_candidate  >=2 qualifying variants in an AR gene (phase unknown)
      AR_second_hit_splice  the Introme-supported hit of an AR comp-het (partner is AR_comphet_candidate)
      AR_single_het         only one qualifying het in an AR gene (look for missed 2nd hit / CNV)
      XL_hemi               qualifying hemizygous variant in an X-linked gene
      risk_carrier          qualifying variant in a risk-factor gene (e.g. GBA1)
      risk_biallelic        homozygous or >=2 qualifying variants in a risk gene (e.g. GBA1 -> Gaucher)
      MOI_unknown_candidate qualifying variant in a panel gene with inheritance "unknown"
      not_qualifying        below the tier threshold for its MOI, or not on the panel

    When ``qualify_clinvar_vus`` is on (default), ClinVar VUS alleles also qualify
    regardless of tier (coding and non-coding). Ranking prefers P/LP, then coding
    VUS, then non-coding VUS.
    """
    df = df.copy()
    df["inheritance_status"] = "not_qualifying"
    df["ar_hit_count"] = 0

    ad_q = _qualifies(df, icfg["ad_min_tier"])
    ar_first = _qualifies(df, icfg["ar_first_hit_min_tier"])
    ar_second = _qualifies(df, icfg["ar_second_hit_min_tier"])
    if icfg.get("allow_introme_second_hit", True) and "introme_class" in df.columns:
        ar_second |= df["introme_class"] == "high"
    # ClinVar VUS (coding or non-coding) qualify below the usual tier floor.
    if icfg.get("qualify_clinvar_vus", True):
        vus_q = clinvar_vus_mask(df)
        ad_q = ad_q | vus_q
        ar_first = ar_first | vus_q
        ar_second = ar_second | vus_q

    df.loc[(df["moi"] == "AD") & ad_q, "inheritance_status"] = "AD_candidate"
    df.loc[(df["moi"] == "RISK") & ad_q, "inheritance_status"] = "risk_carrier"
    df.loc[(df["moi"] == "XL") & ad_q & (df["zygosity"] == "hemi"), "inheritance_status"] = "XL_hemi"
    df.loc[(df["moi"] == "XL") & ad_q & (df["zygosity"] != "hemi"), "inheritance_status"] = "AD_candidate"
    unknown_q = _qualifies(df, icfg.get("unknown_moi_min_tier", icfg["ad_min_tier"]))
    if icfg.get("qualify_clinvar_vus", True):
        unknown_q = unknown_q | clinvar_vus_mask(df)
    df.loc[(df["moi"] == "UNKNOWN") & unknown_q, "inheritance_status"] = "MOI_unknown_candidate"

    for (_, _), grp in df[(df["moi"] == "RISK") & ad_q].groupby(["sample_id", "gene"]):
        if len(grp) >= 2 or (grp["zygosity"] == "hom").any():
            df.loc[grp.index, "inheritance_status"] = "risk_biallelic"

    is_ar = df["moi"] == "AR"
    ar_df = df[is_ar & (ar_first | ar_second)]
    for (_, _), grp in ar_df.groupby(["sample_id", "gene"]):
        idx = grp.index
        first = ar_first.loc[idx].to_numpy()
        second = ar_second.loc[idx].to_numpy()
        hom = first & (grp["zygosity"] == "hom").to_numpy()
        n_hits = int(second.sum())
        df.loc[idx, "ar_hit_count"] = n_hits
        if hom.any():
            df.loc[idx[hom], "inheritance_status"] = "AR_hom"
        # Non-hom alleles still need labels. A qualifying AR_hom (including an
        # Introme intron hom) must not leave coding partners as not_qualifying.
        other = ~hom
        if not other.any():
            continue
        first_o = first & other
        second_o = second & other
        # Do not let a same-gene hom (often a weak Introme intron) inflate the
        # non-hom hit count into a false AR_comphet_candidate.
        n_hits_other = int(second_o.sum())
        introme = grp["evidence_reasons"].fillna("").str.contains("introme_high").to_numpy()
        if first_o.any() and n_hits_other >= 2:
            splice_hit = second_o & (~first_o | introme)
            df.loc[idx[second_o], "inheritance_status"] = "AR_comphet_candidate"
            if splice_hit.any() and (second_o & ~splice_hit).any():
                df.loc[idx[splice_hit], "inheritance_status"] = "AR_second_hit_splice"
        elif first_o.any():
            df.loc[idx[first_o], "inheritance_status"] = "AR_single_het"

    return demote_gof_ad_intronic(df, icfg)
