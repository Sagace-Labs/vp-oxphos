# vp-oxphos

Predicts collapse of the mitochondrial membrane potential from a SMILES string
— a molecular initiating event for drug-induced liver injury, where a
hepatocyte that cannot hold its inner-membrane proton gradient loses ATP
production and floods with reactive oxygen species.

## Install

    pip install vp-oxphos

## Use

    from vp_oxphos import predict
    predict(["CC(=O)Oc1ccccc1C(=O)O"])      # -> DataFrame[oxphos_disrupt]

Returns one row per input and one column per declared output: a probability
that the compound disrupts the membrane potential. Several distinct
mechanisms reach that endpoint and the assay does not separate them, so read
the score as a mitochondrial liability rather than as one mode of action.
Unparseable SMILES come back as NaN. Pin a version with
`predict(smiles, version="v1")`; list what is available with `versions()`.

## Current version

**v1**, signature 1, measured under protocol `scaffold-shuffle-5seed@1`. The
full record — metrics per seed, dataset hash, environment — is in
[`src/vp_oxphos/versions/v1/CARD.md`](src/vp_oxphos/versions/v1/CARD.md).

## Data

PubChem BioAssay AID 720635, the Tox21 qHTS screen for disruptors of the
mitochondrial membrane potential, reduced to one row per compound labelled by
the majority call across its assay records, retrieved 2026-09-01 and
redistributed here as a United States government work in the public domain.
Rebuild it and check for upstream drift with
`python -m vp_oxphos.data fetch --verify`; see [`data/README.md`](data/README.md)
for the expected layout.

## Retrain

    python -m vp_oxphos.train --version v2 --reason "why this version exists"
    python -m vp_oxphos.evaluate --version v2

`train` fits the deployment model on the whole dataset and writes a new version
directory; `evaluate` refits per seed under the protocol and records what those
held-out models scored. Reproducibility is to the recorded dataset hash and
environment, which can change — the source is a live endpoint and RDKit descriptor
values move between releases.

## Licence

Code is Apache-2.0 ([`LICENSE`](LICENSE)). The bundled dataset is in the public
domain ([`LICENSE-DATA`](LICENSE-DATA)).

## Cite

See [`CITATION.cff`](CITATION.cff).
