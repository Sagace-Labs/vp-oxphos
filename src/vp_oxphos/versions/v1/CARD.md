# oxphos v1

_Generated from `manifest.toml` and `metrics.json`. Do not edit._

Released 2026-09-01 · signature 1

**Why this version.** Initial release. Binary XGBoost on Morgan, MACCS and RDKit descriptors over the Tox21 membrane-potential screen, giving the panel its mitochondrial axis alongside the cholestatic one.

## Outputs

| column | dtype | range | meaning |
|---|---|---|---|
| `oxphos_disrupt` | float32 | 0.0–1.0 | P(collapses the mitochondrial membrane potential in the Tox21 qHTS reporter). The assay does not separate the mechanisms that reach it |

Missing values: NaN when RDKit cannot parse the input SMILES

## Performance

Protocol `scaffold-shuffle-5seed@1` — Bemis-Murcko scaffold split with scaffold groups permuted by seed, so distinct seeds give distinct test sets. Five seeds; report mean and standard deviation over the held-out test folds.

Evaluated 2026-09-01 on n_train=3694, n_val=268, n_test=963.

| metric | mean | std | per seed |
|---|---|---|---|
| auc_roc | 0.887 | 0.017 | 0.866, 0.884, 0.895, 0.915, 0.873 |
| auprc | 0.784 | 0.032 | 0.738, 0.776, 0.812, 0.828, 0.767 |
| mcc | 0.590 | 0.045 | 0.556, 0.548, 0.622, 0.662, 0.561 |
| brier | 0.120 | 0.007 | 0.128, 0.117, 0.123, 0.108, 0.125 |

> Comparable only with metrics carrying the same protocol id.

## Data

Tox21 mitochondrial membrane potential qHTS — PubChem BioAssay AID 720635 (qHTS assay for small molecule disruptors of the mitochondrial membrane potential), rows called Active or Inactive, one row per compound labelled by majority call across its assay records. Retrieved 2026-09-01, licensed public-domain, redistributed here.

`4925` compounds, positive rate `0.263`, table SHA-256 `f86374d71b3db436…`

Regenerate and check for upstream drift with `python -m vp_oxphos.data fetch --verify`.

## Model

xgboost-binary on `combo3` features. Shipped weights: full dataset minus a 10% scaffold carve used for early stopping

`weights.joblib` SHA-256 `1a0581ef95308c72…`

## Provenance

Environment: python 3.11.11, rdkit 2026.03.5, xgboost 3.2.0.

Reproducibility is to this dataset hash and this environment, not bit-exact: the sources are live endpoints and RDKit descriptor values move between releases.
