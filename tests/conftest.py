from pathlib import Path

import pytest

from pdvar.config import load_config

FIXTURES = Path(__file__).parent / "fixtures"


@pytest.fixture
def cfg(tmp_path):
    """Project config pointed at the synthetic fixtures, writing outputs to tmp_path."""
    c = load_config()
    c["paths"].update(
        raw_dir=str(FIXTURES),
        interim_dir=str(tmp_path / "interim"),
        processed_dir=str(tmp_path / "processed"),
        results_dir=str(tmp_path / "results"),
    )
    c["inputs"] = {
        "introme": str(FIXTURES / "introme" / "*.csv"),
        "phenotype": str(FIXTURES / "phenotype" / "*.csv"),
        "panel": str(FIXTURES / "panel" / "*.csv"),
    }
    return c
