import pandas as pd

from pdvar import io, splicing
from pdvar.normalize import normalize_introme


def _parsed(cfg):
    raw = io.load_raw(cfg, "introme")
    return normalize_introme(raw, cfg["columns"]["introme"], cfg.get("gene_aliases"))


def test_parse_introme_variants(cfg):
    variants, _ = _parsed(cfg)
    assert variants["variant_id"].is_unique
    row = variants.set_index("variant_id").loc["6-161500000-TA-T"]
    assert row["introme_genes"] == "ENSG00000286805;PRKN"
    assert row["introme_n_carriers"] == 2  # duplicate rows merged
    common = variants.set_index("variant_id").loc["1-20670000-C-T"]
    assert common["introme_cohort_ac"] == 23
    assert common["introme_n_hom"] == 11
    assert variants.set_index("variant_id").loc["1-20660000-A-G", "introme_spliceai_max"] == 0.4


def test_parse_introme_carriers(cfg):
    _, carriers = _parsed(cfg)
    multi = carriers[carriers["variant_id"] == "6-161500000-TA-T"]
    assert set(multi["gene"]) == {"ENSG00000286805", "PRKN"}
    assert set(multi["sample_id"]) == {"P005", "X1"}
    hom = carriers[(carriers["sample_id"] == "S10")]
    assert hom["zygosity"].tolist() == ["hom"]


def test_merge_introme_scores_onto_callset_alleles_only(cfg):
    """Introme annotates existing callset rows; it does not invent missing alleles."""
    variants, _carriers = _parsed(cfg)
    scored = pd.DataFrame(
        {
            "sample_id": ["P003"],
            "gene": ["PINK1"],
            "variant_id": ["1-20660000-A-G"],
            "consequence": "intron_variant",
            "evidence_score": 0,
            "evidence_reasons": "",
        }
    )
    merged = splicing.merge_introme(scored, variants)
    out = splicing.apply_introme_evidence(merged, cfg["splicing"], cfg["scoring"])
    assert len(out) == 1
    assert out.loc[0, "variant_id"] == "1-20660000-A-G"
    assert out.loc[0, "introme_class"] == "high"
