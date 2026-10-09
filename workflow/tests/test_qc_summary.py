"""Single-sample QC summary."""

import json
from pathlib import Path

from qc_summary import build_summary, main, parse_somalier_sex


def _write(path: Path, text: str) -> Path:
    path.write_text(text, encoding="utf-8")
    return path


def test_provided_sex_wins_and_freemix_fails(tmp_path):
    selfsm = _write(tmp_path / "selfSM", "# comment\n#SEQ_ID\tFREEMIX\nSYN\t0.05\n")
    somalier = _write(tmp_path / "samples.tsv", "sample_id\tpredicted_sex\tx_n\tx_het\nSYN\tfemale\t100\t40\n")
    mosdepth = _write(tmp_path / "mosdepth.tsv", "chrom\tlength\tbases\tmean\ntotal_region\t1000\t30000\t30\n")
    reference = _write(tmp_path / "ref.json", json.dumps({"ok": True}))
    summary, sex = build_summary(
        sample_id="SYNTHETIC_SAMPLE",
        selfsm=selfsm,
        somalier=somalier,
        mosdepth=mosdepth,
        reference=reference,
        provided_sex="male",
        max_freemix=0.02,
        min_mean_coverage=20,
        hard_fail_contamination=True,
        hard_fail_sex=False,
        hard_fail_coverage=False,
    )
    assert sex == "male"
    assert summary["sex"]["source"] == "provided"
    assert summary["pass"] is False
    assert summary["coverage"]["mean"] == 30
    assert summary["relatedness"] == "not_run"


def test_x_het_rate_infers_female(tmp_path):
    somalier = _write(tmp_path / "samples.tsv", "sample_id\tx_n\tx_het\nSYN\t100\t40\n")
    assert parse_somalier_sex(somalier) == "female"
    out = tmp_path / "qc.json"
    sex_out = tmp_path / "sex.txt"
    code = main(
        [
            "--sample-id",
            "SYNTHETIC_SAMPLE",
            "--selfsm",
            str(_write(tmp_path / "selfSM", "#SEQ_ID\tFREEMIX\nSYN\t0.01\n")),
            "--somalier",
            str(somalier),
            "--mosdepth",
            str(_write(tmp_path / "m.tsv", "chrom\tmean\ntotal\t25\n")),
            "--reference",
            str(_write(tmp_path / "ref.json", "{}")),
            "--provided-sex",
            "unspecified",
            "--out",
            str(out),
            "--sex-out",
            str(sex_out),
        ]
    )
    assert code == 0
    assert sex_out.read_text(encoding="utf-8") == "female\n"
    assert json.loads(out.read_text(encoding="utf-8"))["sex"]["value"] == "female"
