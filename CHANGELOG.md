# Changelog

Released versions are immutable. A correction to a released version is a new
patch version. This file records *why* each version exists.

## v3

Adds seventeen published structural alerts for mitochondrial toxicity as
substructure counts beside the descriptor panel.

## v2

Adds the viability counter-screen as a second endpoint. A compound that kills a
cell collapses its membrane potential as a consequence. The counter-screen
calls 4579 of the table's compounds; the rest carry a null, which records that
no call was made rather than a negative, and each endpoint is fit on the
compounds it labels.

Features drop the Morgan fingerprint and MACCS keys for the RDKit descriptor
panel.

Signature 2: `oxphos_cytotox` joins `oxphos_disrupt`.

## v1

Initial release.

Signature 1: a single column `oxphos_disrupt`.
