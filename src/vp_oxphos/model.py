"""The OXPHOS model: binary XGBoost on descriptors and structural alerts.

Hyperparameters, feature choice and featuriser."""

from __future__ import annotations

from typing import Any

import numpy as np

from vp_core import fingerprints, xgb
from vp_oxphos import alerts

__all__ = ["FEATURES", "HYPERPARAMS", "featurize", "fit", "predict"]

FEATURES = "rdkit_desc+alerts"

HYPERPARAMS: dict[str, Any] = {
    "n_estimators": 2000,
    "learning_rate": 0.05,
    "max_depth": 6,
    "min_child_weight": 1.0,
    "subsample": 0.8,
    "colsample_bytree": 0.8,
    "reg_lambda": 2.0,
    "reg_alpha": 0.05,
    "gamma": 0.1,
    "early_stopping_rounds": 40,
}


def featurize(smiles: list[str], kind: str = FEATURES) -> np.ndarray:
    """Feature matrix for ``kind``. Names vp-core does not know are built here."""
    if kind == "rdkit_desc+alerts":
        return np.hstack(
            [fingerprints.featurize(smiles, "rdkit_desc"), alerts.counts(smiles)]
        ).astype(np.float32)
    return fingerprints.featurize(smiles, kind)


def fit(X_train, y_train, X_val, y_val, *, seed: int = 0) -> Any:
    """Fit one endpoint's binary model. The validation fold stops boosting."""
    return xgb.fit_binary(
        X_train,
        y_train,
        X_val,
        y_val,
        params=dict(HYPERPARAMS),
        seed=seed,
    )


def predict(model: Any, smiles: list[str]) -> np.ndarray:
    """Positive-class probability for arbitrary SMILES."""
    return xgb.predict_proba(model, featurize(smiles))
