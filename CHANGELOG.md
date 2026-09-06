# Changelog

Released versions are immutable. A correction to a released version is a new
patch version. This file records *why* each version exists; what
it contains is in its own `manifest.toml` and `CARD.md`.

## v2

Adds the viability counter-screen as a second endpoint. The two readouts are
confounded by construction — a compound that kills a cell collapses its membrane
potential as a consequence — and a quarter of the compounds labelled positive on
the membrane potential are also active in the counter-screen. One output cannot
say which of the two a compound is doing, so v2 predicts both and leaves the
comparison to the reader.

The counter-screen calls 4579 of the table's compounds. The rest keep a null in
`cytotox`, which is an absence of evidence rather than a negative, and each
endpoint is fit and scored only on the compounds it labels.

Features drop the Morgan fingerprint and MACCS keys for the RDKit descriptor
panel alone: 217 columns against 2432, for a held-out difference smaller than
the spread across seeds. This endpoint is reached by partitioning into a
membrane and carrying a proton back out, which bulk properties state directly
and a substructure key can only approximate.

Signature 2: `oxphos_cytotox` joins `oxphos_disrupt`. Adding a column is
breaking for anyone consuming positionally.

## v1

Initial release.

The label is the assay's own Active/Inactive call rather than a threshold on
the reported potency, because the qHTS curve reports a potency for a minority
of records. The potency travels in the table as provenance, so a later censored
fit can use it without re-fetching and moving the dataset hash.

An Inconclusive call is dropped rather than counted as a negative: a curve the
assay declined to call either way is an absence of evidence.

Signature 1: a single column `oxphos_disrupt`.
