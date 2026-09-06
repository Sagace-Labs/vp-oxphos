"""Mitochondrial membrane-potential disruption from a SMILES string.

Collapsing the inner-membrane proton gradient starves the hepatocyte of ATP and
bursts reactive oxygen species. It is the shared endpoint of several distinct
mitochondrial liabilities and a principal route to drug-induced liver injury.

    from vp_oxphos import predict
    predict(["CC(=O)Oc1ccccc1C(=O)O"])      # -> DataFrame[oxphos_disrupt]
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import pandas as pd

from vp_core.registry import Version, VersionedPathway
from vp_oxphos.target import TARGET, TARGETS, Endpoint, all_names
from vp_oxphos.target import get as get_target

__version__ = "1.0.0"

PATHWAY = "oxphos"
VERSIONS_DIR = Path(__file__).resolve().parent / "versions"

# One probability column from an XGBoost model on a vp-core fingerprint is
# exactly what the registry's default predictor produces, so no override.
_pathway = VersionedPathway(PATHWAY, VERSIONS_DIR)

__all__ = [
    "PATHWAY",
    "TARGET",
    "TARGETS",
    "VERSIONS_DIR",
    "Endpoint",
    "Version",
    "__version__",
    "all_names",
    "current_version",
    "get",
    "get_target",
    "predict",
    "signature",
    "versions",
]


def predict(smiles: list[str], *, version: str | None = None) -> pd.DataFrame:
    """P(membrane-potential disruptor) for each SMILES, using ``version``."""
    return _pathway.predict(smiles, version=version)


def versions() -> list[str]:
    """Released version names, oldest first."""
    return _pathway.versions()


def current_version() -> str:
    """The newest released version."""
    return _pathway.current()


def get(version: str | None = None) -> Version:
    """Load a version and its validated manifest."""
    return _pathway.get(version)


def signature(version: str | None = None) -> dict[str, Any]:
    """The output contract of a version."""
    return get(version).manifest["signature"]
