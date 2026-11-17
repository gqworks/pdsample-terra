from pathlib import Path

import pytest

from pdvar.config import load_config

FIXTURES = Path(__file__).parent / "fixtures"


@pytest.fixture
def cfg():
    """Project config with Introme fixtures for local table-load tests."""
    c = load_config()
    c["inputs"] = {
        **c.get("inputs", {}),
        "introme": str(FIXTURES / "introme" / "*.csv"),
    }
    return c
