"""mtDNA MitImpact evidence and status promotion."""

import pandas as pd

from pdvar import mtdna
from pdvar.config import load_config


def test_apply_mtdna_evidence_adds_mitomap_apogee_mitotip():
    cfg = load_config()
    df = pd.DataFrame(
        [
            {
                "hit_source": "mtdna",
                "evidence_score": 2,
                "evidence_reasons": "absent_gnomad",
                "tier": 4,
                "mitomap_status": "Cfrm-[P]",
                "apogee2_class": "Pathogenic",
                "mitotip_quartile": "",
                "mlc_score": 0.72,
                "clinvar_class": "P/LP",
                "gnomad_af": 0.0002,
            },
            {
                "hit_source": "small_variant",
                "evidence_score": 10,
                "evidence_reasons": "clinvar_plp",
                "tier": 1,
                "mitomap_status": "",
                "apogee2_class": "",
                "mitotip_quartile": "Q1",
                "clinvar_class": "P/LP",
                "gnomad_af": 0,
            },
        ]
    )
    out = mtdna.apply_mtdna_evidence(df, cfg["scoring"])
    mito = out[out["hit_source"] == "mtdna"].iloc[0]
    assert "mitomap_pathogenic" in mito["evidence_reasons"]
    assert "apogee_pathogenic" in mito["evidence_reasons"]
    assert "mlc_high" in mito["evidence_reasons"]
    assert mito["evidence_score"] == 2 + 8 + 3 + 1
    assert mito["tier"] == 1
    nuclear = out[out["hit_source"] == "small_variant"].iloc[0]
    assert nuclear["evidence_score"] == 10
    assert "mitotip_likely_pathogenic" not in nuclear["evidence_reasons"]


def test_promote_mtdna_status_pathogenic_candidate():
    df = pd.DataFrame(
        [
            {
                "hit_source": "mtdna",
                "inheritance_status": "mtdna_variant",
                "clinvar_class": "P/LP",
                "mitomap_status": "",
                "tier": 4,
                "gnomad_af": 0.5,
                "gnomad_af_het": 0.0,
            },
            {
                "hit_source": "mtdna",
                "inheritance_status": "mtdna_variant",
                "clinvar_class": "none",
                "mitomap_status": "",
                "tier": 4,
                "gnomad_af": 0.5,
                "gnomad_af_het": 0.0,
            },
            {
                "hit_source": "mtdna",
                "inheritance_status": "mtdna_variant",
                "clinvar_class": "none",
                "mitomap_status": "Cfrm-[LP]",
                "tier": 4,
                "gnomad_af": 0.2,
                "gnomad_af_het": 0.0,
            },
            {
                # Common het-only allele: AF_hom=0 must not count as rare.
                "hit_source": "mtdna",
                "inheritance_status": "mtdna_variant",
                "clinvar_class": "none",
                "mitomap_status": "",
                "tier": 2,
                "gnomad_af": 0.0,
                "gnomad_af_het": 0.05,
            },
        ]
    )
    out = mtdna.promote_mtdna_status(df, max_gnomad_af=0.01)
    assert out.loc[0, "inheritance_status"] == "mtdna_pathogenic_candidate"
    assert out.loc[1, "inheritance_status"] == "mtdna_variant"
    assert out.loc[2, "inheritance_status"] == "mtdna_pathogenic_candidate"
    assert out.loc[3, "inheritance_status"] == "mtdna_variant"


def test_must_keep_only_pathogenic_mtdna():
    from pdvar.report import _must_keep_mask

    df = pd.DataFrame(
        [
            {"inheritance_status": "mtdna_pathogenic_candidate", "clinvar_class": "none", "hit_source": "mtdna"},
            {"inheritance_status": "mtdna_variant", "clinvar_class": "none", "hit_source": "mtdna"},
        ]
    )
    keep = _must_keep_mask(df)
    assert bool(keep.iloc[0]) is True
    assert bool(keep.iloc[1]) is False
