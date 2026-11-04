"""GBA1-PD Browser annotation and score bumps for Gauchian SNVs."""

from pathlib import Path

import pandas as pd

from pdvar import gba1_pd
from pdvar.config import load_config


def _browser(tmp_path: Path) -> Path:
    path = tmp_path / "gba1_pd_browser.tsv"
    path.write_text(
        "\n".join(
            [
                "cdna\tprotein\tvariant_name\trsid\tpd_class\tgd_acmg\tpos_hg38\tallele_frequency\treferences\tcomments",
                "c.681T>G\tp.(Asn227Lys)\tp.N188K\t.\tsevere\tLikely pathogenic\t155238214\t8e-6\t1\t",
                "c.667T>C\tp.(Trp223Arg)\tp.W184R\t.\tsevere\tPathogenic\t155238228\t2e-5\t1\t",
                "c.1093G>A\tp.(Glu365Lys)\tp.E326K\trs.\trisk\tBenign\t155236376\t0.01\t1\t",
                "c.1226A>G\tp.(Asn409Ser)\tp.N370S\trs.\tmild\tLikely pathogenic\t155235843\t0.003\t1\t",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    return path


def test_annotate_and_score_gauchian_only(tmp_path):
    browser = _browser(tmp_path)
    cfg = load_config()
    df = pd.DataFrame(
        [
            {
                "gene": "GBA1",
                "hit_source": "gauchian",
                "hgvsc": "c.681T>G",
                "evidence_score": 10,
                "evidence_reasons": "clinvar_plp",
                "tier": 1,
            },
            {
                "gene": "GBA1",
                "hit_source": "small_variant",
                "hgvsc": "ENST00000368373.8:c.681T>G",
                "evidence_score": 14,
                "evidence_reasons": "clinvar_plp",
                "tier": 1,
            },
            {
                "gene": "GBA1",
                "hit_source": "gauchian",
                "hgvsc": "c.1093G>A",
                "evidence_score": 10,
                "evidence_reasons": "clinvar_plp",
                "tier": 1,
            },
        ]
    )
    out = gba1_pd.annotate_and_score(df, cfg, browser_path=browser)
    gauchian_severe = out[(out["hit_source"] == "gauchian") & (out["hgvsc"] == "c.681T>G")].iloc[0]
    assert gauchian_severe["gba1_pd_class"] == "severe"
    assert gauchian_severe["gba1_gd_acmg"] == "Likely pathogenic"
    assert gauchian_severe["evidence_score"] == 14  # 10 + gba1_pd_severe(4)
    assert "gba1_pd_severe" in gauchian_severe["evidence_reasons"]

    vep = out[out["hit_source"] == "small_variant"].iloc[0]
    assert pd.isna(vep["gba1_pd_class"]) or str(vep["gba1_pd_class"]) in {"", "nan", "<NA>"}
    assert vep["evidence_score"] == 14
    assert "gba1_pd_severe" not in str(vep["evidence_reasons"])

    risk = out[out["hgvsc"] == "c.1093G>A"].iloc[0]
    assert risk["gba1_pd_class"] == "risk"
    assert risk["evidence_score"] == 11  # 10 + gba1_pd_risk(1)
    assert "gba1_pd_risk" in risk["evidence_reasons"]


def test_default_browser_path_resolves_shipped_resource():
    cfg = load_config()
    path = gba1_pd.default_browser_path(cfg)
    assert path is not None
    assert path.is_file()
    table = gba1_pd.load_browser(path)
    assert "c.667T>C" in set(table["cdna"])
    assert table.loc[table["cdna"] == "c.667T>C", "pd_class"].iloc[0] == "severe"
