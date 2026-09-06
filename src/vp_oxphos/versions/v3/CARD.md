# oxphos v3

_Generated from `manifest.toml` and `metrics.json`. Do not edit._

Released 2026-09-03 · signature 2 · supersedes v2

**Why this version.** adds published mitochondrial structural alerts to the descriptor panel

## Outputs

| column | dtype | range | meaning |
|---|---|---|---|
| `oxphos_disrupt` | float32 | 0.0–1.0 | P(collapses the mitochondrial membrane potential in the Tox21 qHTS reporter). The assay does not separate the mechanisms that reach it |
| `oxphos_cytotox` | float32 | 0.0–1.0 | P(reduces viability in the counter-screen over the same library). A compound scoring high on both readouts lost membrane potential in a cell that was also dying |

Missing values: NaN when RDKit cannot parse the input SMILES

## Performance — `scaffold-shuffle-5seed@1`

Protocol `scaffold-shuffle-5seed@1` — Bemis-Murcko scaffold split with scaffold groups permuted by seed, so distinct seeds give distinct test sets. Five seeds; report mean and standard deviation over the held-out test folds.

Evaluated 2026-09-03 on n_train=3694, n_val=268, n_test=963.

### `oxphos_disrupt`

| metric | mean | std | per seed |
|---|---|---|---|
| auc_roc | 0.888 | 0.015 | 0.868, 0.881, 0.894, 0.914, 0.883 |
| auprc | 0.788 | 0.033 | 0.737, 0.775, 0.819, 0.828, 0.779 |
| mcc | 0.588 | 0.043 | 0.521, 0.569, 0.625, 0.644, 0.583 |
| brier | 0.117 | 0.006 | 0.124, 0.116, 0.125, 0.109, 0.113 |

### `oxphos_cytotox`

Measured on n_train=3449, n_val=233, n_test=897.

| metric | mean | std | per seed |
|---|---|---|---|
| auc_roc | 0.880 | 0.010 | 0.865, 0.888, 0.888, 0.885, 0.872 |
| auprc | 0.393 | 0.055 | 0.343, 0.407, 0.493, 0.352, 0.371 |
| mcc | 0.312 | 0.054 | 0.285, 0.234, 0.386, 0.358, 0.299 |
| brier | 0.060 | 0.007 | 0.061, 0.051, 0.053, 0.067, 0.070 |

## Performance — `scaffold-balanced-5seed@1`

Protocol `scaffold-balanced-5seed@1` — Bemis-Murcko scaffold split with scaffold groups permuted by seed and each group placed in the fold it overfills least, so a group larger than a fold's capacity settles in train instead of starving that fold. Same fold fractions, seeds and metrics as scaffold-shuffle-5seed@1; only the packing differs. Five seeds; report mean and standard deviation over the held-out test folds.

Evaluated 2026-09-19 on n_train=3790, n_val=454, n_test=681.

### `oxphos_disrupt`

| metric | mean | std | per seed |
|---|---|---|---|
| auc_roc | 0.906 | 0.005 | 0.911, 0.907, 0.897, 0.903, 0.910 |
| auprc | 0.825 | 0.022 | 0.856, 0.823, 0.840, 0.794, 0.810 |
| mcc | 0.635 | 0.028 | 0.666, 0.625, 0.599, 0.616, 0.671 |
| brier | 0.118 | 0.007 | 0.114, 0.120, 0.128, 0.121, 0.106 |

### `oxphos_cytotox`

Measured on n_train=3553, n_val=416, n_test=610.

| metric | mean | std | per seed |
|---|---|---|---|
| auc_roc | 0.856 | 0.023 | 0.837, 0.878, 0.821, 0.875, 0.868 |
| auprc | 0.348 | 0.039 | 0.394, 0.378, 0.305, 0.365, 0.300 |
| mcc | 0.319 | 0.043 | 0.340, 0.256, 0.287, 0.380, 0.331 |
| brier | 0.066 | 0.005 | 0.070, 0.063, 0.070, 0.070, 0.058 |

> Comparable only with metrics carrying the same protocol id.

## Data

Tox21 mitochondrial membrane potential qHTS — PubChem BioAssay AID 720635 (qHTS assay for small molecule disruptors of the mitochondrial membrane potential) and AID 720634 (qHTS assay for small molecule disruptors of the mitochondrial membrane potential - cell viability), rows called Active or Inactive, one row per compound labelled by majority call across its assay records. Retrieved 2026-09-01, licensed public-domain, redistributed here.

`4925` compounds, positive rate `0.263`, table SHA-256 `95f4fa7fb544025d…`

Regenerate and check for upstream drift with `python -m vp_oxphos.data fetch --verify`.

## Model

xgboost-binary on `rdkit_desc+alerts` features. Shipped weights: one model per output, each on every compound its endpoint labels minus a 10% scaffold carve used for early stopping

`weights.joblib` SHA-256 `20ff8b1b1096b632…`

## Provenance

Environment: python 3.11.11, rdkit 2026.03.5, xgboost 3.2.0.

Reproducibility is to this dataset hash and this environment, not bit-exact: the sources are live endpoints and RDKit descriptor values move between releases.
