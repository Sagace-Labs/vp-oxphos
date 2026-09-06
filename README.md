# vp-oxphos

Predicts collapse of the mitochondrial membrane potential from a SMILES string
— a molecular initiating event for drug-induced liver injury, where a
hepatocyte that cannot hold its inner-membrane proton gradient loses ATP
production and floods with reactive oxygen species.

## Install

    pip install vp-oxphos

## Use

    from vp_oxphos import predict
    predict(["CC(=O)Oc1ccccc1C(=O)O"])   # -> DataFrame[oxphos_disrupt, oxphos_cytotox]

Returns one row per input and one column per declared output. `oxphos_disrupt`
is the probability of collapsing the membrane potential; `oxphos_cytotox` is the
probability of reducing viability in the counter-screen over the same library. A
compound scoring high on both lost membrane potential in a cell that was also
dying. Several distinct mechanisms reach the endpoint and the assay does not
separate them. Unparseable SMILES come back as NaN. Pin a version with `predict(smiles, version="v3")`; list what is
available with `versions()`.

## Current version

**v3**, signature 2, measured under protocol `scaffold-shuffle-5seed@1`. The
full record — metrics per output and per seed, dataset hash, environment — is in
[`src/vp_oxphos/versions/v3/CARD.md`](src/vp_oxphos/versions/v3/CARD.md).

## Data

PubChem BioAssay AID 720635, the Tox21 qHTS screen for disruptors of the
mitochondrial membrane potential, and AID 720634, its cell-viability
counter-screen over the same library. Both are reduced to one row per compound
labelled by the majority call across its assay records, retrieved 2026-09-01 and
redistributed here as a United States government work in the public domain.
Rebuild and check for upstream drift with
`python -m vp_oxphos.data fetch --verify`; see [`data/README.md`](data/README.md)
for the expected layout.

## Retrain

    python -m vp_oxphos.train --version v4 --reason "why this version exists"
    python -m vp_oxphos.evaluate --version v4

`train` fits one deployment model per output on the whole dataset and writes a
new version directory; `evaluate` refits per seed under the protocol and records
what those held-out models scored. Reproducibility is to the recorded dataset
hash and environment, which can change.

## Licence

Code is Apache-2.0 ([`LICENSE`](LICENSE)). The bundled dataset is in the public
domain ([`LICENSE-DATA`](LICENSE-DATA)).

## Cite

See [`CITATION.cff`](CITATION.cff).
