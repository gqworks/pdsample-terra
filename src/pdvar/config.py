"""Load and validate YAML configuration."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import yaml

PROJECT_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_CONFIG = PROJECT_ROOT / "config" / "config.yaml"

REQUIRED_SECTIONS = ("paths", "inputs", "filters", "scoring", "splicing", "inheritance", "phenotype", "report")


def load_yaml(path: str | Path) -> dict[str, Any]:
    with open(path, encoding="utf-8") as fh:
        return yaml.safe_load(fh) or {}


def resolve(path: str | Path) -> Path:
    """Resolve a path relative to the project root."""
    p = Path(path)
    return p if p.is_absolute() else PROJECT_ROOT / p


def load_config(
    path: str | Path | None = None,
    columns_path: str | Path | None = None,
    phenotypes_path: str | Path | None = None,
) -> dict[str, Any]:
    """Load the main config plus the column mapping and gene phenotype profiles.

    Relative paths inside the config are resolved from the repository root when
    the file lives in ``config/``. Explicit ``columns_path`` and
    ``phenotypes_path`` override those lookups, which is how the Assess task
    points at Terra-localized files.
    """
    cfg_path = Path(path) if path else DEFAULT_CONFIG
    if not cfg_path.is_file():
        raise FileNotFoundError(cfg_path)
    cfg = load_yaml(cfg_path)
    missing = [section for section in REQUIRED_SECTIONS if section not in cfg]
    if missing:
        raise ValueError(f"Config missing sections: {missing}")
    root = cfg_path.resolve().parent.parent if cfg_path.resolve().parent.name == "config" else cfg_path.resolve().parent
    cfg["columns"] = load_yaml(_locate(cfg["paths"]["columns"], root, columns_path))
    cfg["gene_phenotypes"] = load_yaml(_locate(cfg["paths"]["gene_phenotypes"], root, phenotypes_path))
    cfg["_root"] = str(root)
    return cfg


def _locate(configured: str, root: Path, override: str | Path | None) -> Path:
    if override:
        path = Path(override)
        if not path.is_file():
            raise FileNotFoundError(path)
        return path
    raw = Path(str(configured))
    candidates = [raw] if raw.is_absolute() else [root / raw, Path.cwd() / raw]
    for candidate in candidates:
        if candidate.is_file():
            return candidate
    raise FileNotFoundError(f"Could not find {configured} from {root}")
