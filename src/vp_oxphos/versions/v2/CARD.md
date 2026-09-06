# oxphos v2

_Generated from `manifest.toml` and `metrics.json`. Do not edit._

Released 2026-09-03 · signature 2

**Why this version.** adds the viability counter-screen as a second endpoint and output; drops the structural fingerprint for the descriptor panel

## Outputs

| column | dtype | range | meaning |
|---|---|---|---|
| `oxphos_disrupt` | float32 | 0.0–1.0 | P(collapses the mitochondrial membrane potential in the Tox21 qHTS reporter). The assay does not separate the mechanisms that reach it |
| `oxphos_cytotox` | float32 | 0.0–1.0 | P(reduces viability in the counter-screen over the same library). A compound scoring high on both readouts lost membrane potential in a cell that was also dying |

Missing values: NaN when RDKit cannot parse the input SMILES

## Performance

Protocol `scaffold-shuffle-5seed@1` — Bemis-Murcko scaffold split with scaffold groups permuted by seed, so distinct seeds give distinct test sets. Five seeds; report mean and standard deviation over the held-out test folds.

Evaluated 2026-09-03 on n_train=3694, n_val=268, n_test=963.

### `oxphos_disrupt`

| metric | mean | std | per seed |
|---|---|---|---|
| auc_roc | 0.883 | 0.018 | 0.866, 0.874, 0.886, 0.917, 0.871 |
| auprc | 0.776 | 0.040 | 0.724, 0.764, 0.798, 0.839, 0.754 |
| mcc | 0.588 | 0.037 | 0.561, 0.565, 0.613, 0.648, 0.553 |
| brier | 0.120 | 0.008 | 0.124, 0.118, 0.131, 0.107, 0.121 |

### `oxphos_cytotox`

Measured on n_train=3449, n_val=233, n_test=897.

| metric | mean | std | per seed |
|---|---|---|---|
| auc_roc | 0.872 | 0.019 | 0.862, 0.889, 0.893, 0.875, 0.841 |
| auprc | 0.370 | 0.057 | 0.363, 0.380, 0.471, 0.333, 0.305 |
| mcc | 0.331 | 0.087 | 0.308, 0.268, 0.486, 0.354, 0.239 |
| brier | 0.064 | 0.017 | 0.059, 0.051, 0.047, 0.069, 0.096 |

> Comparable only with metrics carrying the same protocol id.

## Data

Tox21 mitochondrial membrane potential qHTS — PubChem BioAssay AID 720635 (qHTS assay for small molecule disruptors of the mitochondrial membrane potential) and AID 720634 (qHTS assay for small molecule disruptors of the mitochondrial membrane potential - cell viability), rows called Active or Inactive, one row per compound labelled by majority call across its assay records. Retrieved 2026-09-01, licensed public-domain, redistributed here.

`4925` compounds, positive rate `0.263`, table SHA-256 `95f4fa7fb544025d…`

Regenerate and check for upstream drift with `python -m vp_oxphos.data fetch --verify`.

## Model

xgboost-binary on `rdkit_desc` features. Shipped weights: one model per output, each on every compound its endpoint labels minus a 10% scaffold carve used for early stopping

`weights.joblib` SHA-256 `90ef550c0ae64674…`

## Provenance

Environment: python 3.11.11, rdkit 2026.03.5, xgboost 3.2.0.

Reproducibility is to this dataset hash and this environment, not bit-exact: the sources are live endpoints and RDKit descriptor values move between releases.
