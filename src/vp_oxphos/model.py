"""The OXPHOS model: binary XGBoost on a Morgan, MACCS and descriptor stack.

Hyperparameters and feature choice.

The label is the assay's own Active/Inactive call, so the fit is binary rather
than censored: the qHTS potency is reported for too few compounds to carry an
interval. Balanced positive weighting comes from ``vp_core.xgb``, since a
quarter of the library is active.
"""

from __future__ import annotations

from typing import Any

import numpy as np

from vp_core import fingerprints, xgb

__all__ = ["FEATURES", "HYPERPARAMS", "fit", "predict"]

FEATURES = "combo3"

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


def fit(X_train, y_train, X_val, y_val, *, seed: int = 0) -> Any:
    """Fit one binary model. The validation fold stops boosting."""
    return xgb.fit_binary(
        X_train,
        y_train,
        X_val,
        y_val,
        params=dict(HYPERPARAMS),
        seed=seed,
    )


def predict(model: Any, smiles: list[str]) -> np.ndarray:
    """P(membrane-potential disruptor) for arbitrary SMILES."""
    X = fingerprints.featurize(smiles, FEATURES)
    return xgb.predict_proba(model, X)
