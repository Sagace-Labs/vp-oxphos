"""Mitochondrial membrane-potential disruption from a SMILES string.

Collapsing the inner-membrane proton gradient starves the hepatocyte of ATP and
bursts reactive oxygen species.

    from vp_oxphos import predict
    predict(["CC(=O)Oc1ccccc1C(=O)O"])   # -> DataFrame[oxphos_disrupt, oxphos_cytotox]
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from vp_core.registry import Version, VersionedPathway
from vp_oxphos.target import CYTOTOX, TARGET, TARGETS, Endpoint, all_names
from vp_oxphos.target import get as get_target

__version__ = "2.1.0"

PATHWAY = "oxphos"
VERSIONS_DIR = Path(__file__).resolve().parent / "versions"


def _predict_values(model: Any, smiles: list[str], version: Version) -> np.ndarray:
    """Columns for ``version``, in the order its signature declares them.

    Signature 1 stores one model; signature 2 stores one per endpoint under
    its output name.
    """
    from rdkit import Chem, RDLogger

    from vp_core import xgb
    from vp_oxphos import model as oxphos_model

    RDLogger.DisableLog("rdApp.*")
    X = oxphos_model.featurize(smiles, str(version.features))
    if isinstance(model, dict):
        values = np.column_stack(
            [xgb.predict_proba(model[name], X) for name in version.output_names]
        )
    else:
        values = xgb.predict_proba(model, X).reshape(-1, 1)

    # An unparseable input is a declared NaN, not an error.
    unparseable = [Chem.MolFromSmiles(s) is None for s in smiles]
    values = values.astype(np.float32)
    values[np.asarray(unparseable)] = np.nan
    return values


_pathway = VersionedPathway(PATHWAY, VERSIONS_DIR, predict_fn=_predict_values)

__all__ = [
    "CYTOTOX",
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
    """Score each SMILES with ``version`` (default: newest)."""
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
