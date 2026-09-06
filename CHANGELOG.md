# Changelog

Released versions are immutable. A correction to a released version is a new
patch version. This file records *why* each version exists; what
it contains is in its own `manifest.toml` and `CARD.md`.

## v1

Initial release.

The label is the assay's own Active/Inactive call rather than a threshold on
the reported potency, because the qHTS curve reports a potency for a minority
of records. The potency travels in the table as provenance, so a later censored
fit can use it without re-fetching and moving the dataset hash.

An Inconclusive call is dropped rather than counted as a negative: a curve the
assay declined to call either way is an absence of evidence.

Signature 1: a single column `oxphos_disrupt`.
