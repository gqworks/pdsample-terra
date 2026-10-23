"""Regression for EnsembleInfer header column lookup under bash pipefail."""

from __future__ import annotations

import gzip
import subprocess
import textwrap
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]


def _write_predictions(path: Path, rows: int) -> None:
    with gzip.open(path, "wt", encoding="utf-8") as handle:
        handle.write("CHROM\tPOS\tREF\tALT\tintrome_score\n")
        for i in range(rows):
            handle.write(f"chr1\t{i}\tA\tG\t0.{i % 10}\n")


def test_zcat_head_under_pipefail_exits_141_on_large_tsv(tmp_path: Path) -> None:
    tsv = tmp_path / "introme.predictions.tsv.gz"
    _write_predictions(tsv, 200_000)
    script = textwrap.dedent(
        f"""\
        set -euo pipefail
        cd "{tmp_path}"
        introme_col=$(
            zcat introme.predictions.tsv.gz \\
            | head -1 \\
            | tr '\\t' '\\n' \\
            | awk '$0 == "introme_score" {{print NR; exit}}'
        )
        printf '%s\\n' "deadbeef" > introme.commit.txt
        echo "col=$introme_col"
        """
    )
    result = subprocess.run(["bash", "-c", script], capture_output=True, text=True)
    assert result.returncode == 141
    assert not (tmp_path / "introme.commit.txt").exists()


def test_python_header_lookup_reaches_commit_txt(tmp_path: Path) -> None:
    tsv = tmp_path / "introme.predictions.tsv.gz"
    _write_predictions(tsv, 200_000)
    # Same logic as EnsembleInfer after the fix.
    script = textwrap.dedent(
        f"""\
        set -euo pipefail
        cd "{tmp_path}"
        introme_col="$(python3 - <<'PY'
import gzip
import sys

with gzip.open("introme.predictions.tsv.gz", "rt", encoding="utf-8", newline="") as handle:
    header = handle.readline().rstrip("\\n").split("\\t")
try:
    print(header.index("introme_score") + 1)
except ValueError:
    sys.stderr.write("introme_score column missing from ensemble TSV header: %r\\n" % (header,))
    sys.exit(1)
PY
)"
        printf '%s\\n' "deadbeef" > introme.commit.txt
        echo "col=$introme_col"
        """
    )
    result = subprocess.run(["bash", "-c", script], capture_output=True, text=True, check=True)
    assert "col=5" in result.stdout
    assert (tmp_path / "introme.commit.txt").read_text(encoding="utf-8") == "deadbeef\n"
