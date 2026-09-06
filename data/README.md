# OXPHOS data

    data/
      oxphos_tox21.parquet         standardised table (producing the hash)
      example/
        oxphos_example.parquet     small stratified fixture, used by the tests

`oxphos_tox21.parquet` holds one row per compound with the three contract
columns — `inchikey`, `smiles` (standardised), `label` — plus `potency_um`,
`n_calls` and `active_frac` as provenance. Only the contract columns are
hashed, so a provenance column may be added without moving the dataset hash.

## Origin and processing

PubChem BioAssay AID 720635, the Tox21 qHTS screen for disruptors of the
mitochondrial membrane potential, retrieved through the public PUG-REST
concise assay endpoint. Rows the assay called Active or Inactive are kept and
rows it called Inconclusive are dropped. Compound identifiers are resolved to
isomeric SMILES through the compound property endpoint; structures are
standardised (normalise, largest fragment, neutralise); records are collapsed
to one row per InChIKey by majority call, and a tie is dropped.

A standardised structure that RDKit will not read back is dropped too, in order
to avoid featurization outputting all-zeros.

`potency_um` is the median of the potencies reported across a compound's
records and is NaN when none reported one, which is the usual case. `n_calls`
is how many assay records collapsed into the row and `active_frac` is the share
of them that were Active, so the majority vote stays auditable.

## Rebuilding it

    python -m vp_oxphos.data fetch --verify

This downloads, re-parses and re-hashes, then compares against the hash the
current version recorded. A mismatch in the hash can indicate upstream changes
to the source data or the processing pipeline.

## Licence

The bundled table is a United States government work in the public domain; see
`../LICENSE-DATA`. That file covers these files only, not the package code or
the trained weights.
