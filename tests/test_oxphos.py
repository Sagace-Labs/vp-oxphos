"""OXPHOS package tests.

Everything here runs on the committed fixture and the shipped weights, so a
fresh clone with no network and no licensed data can still verify the package.
"""

from __future__ import annotations

import numpy as np
import pytest

import vp_oxphos
from vp_oxphos import contract, data

ASPIRIN = "CC(=O)Oc1ccccc1C(=O)O"
CAFFEINE = "Cn1cnc2c1c(=O)n(C)c(=O)n2C"
NONSENSE = "not-a-molecule"


def test_versions_are_discoverable_and_ordered():
    versions = vp_oxphos.versions()
    assert versions, "no released versions found"
    assert vp_oxphos.current_version() == versions[-1]


def test_declared_signature_matches_the_shipped_manifest():
    shipped = vp_oxphos.signature()
    assert shipped["version"] == contract.SIGNATURE_VERSION
    assert [o["name"] for o in shipped["outputs"]] == contract.column_names()


def test_predictions_carry_the_declared_column_names():
    frame = vp_oxphos.predict([ASPIRIN, CAFFEINE])
    assert list(frame.columns) == contract.column_names()
    assert len(frame) == 2


def test_probabilities_are_in_range():
    values = vp_oxphos.predict([ASPIRIN, CAFFEINE])["oxphos_disrupt"].to_numpy()
    assert np.all((values >= 0.0) & (values <= 1.0))


def test_unparseable_input_becomes_nan_rather_than_raising():
    frame = vp_oxphos.predict([ASPIRIN, NONSENSE])
    assert np.isfinite(frame["oxphos_disrupt"].iloc[0])
    assert np.isnan(frame["oxphos_disrupt"].iloc[1])


def test_a_version_can_be_pinned():
    """Every release stays loadable under its own contract, not the newest one."""
    for name in vp_oxphos.versions():
        frame = vp_oxphos.predict([ASPIRIN], version=name)
        declared = [o["name"] for o in vp_oxphos.signature(name)["outputs"]]
        assert list(frame.columns) == declared

    assert vp_oxphos.predict([ASPIRIN]).equals(
        vp_oxphos.predict([ASPIRIN], version=vp_oxphos.current_version())
    )


def test_a_single_string_is_rejected():
    with pytest.raises(TypeError, match="list of SMILES"):
        vp_oxphos.predict(ASPIRIN)  # type: ignore[arg-type]


def test_unknown_version_names_the_available_ones():
    with pytest.raises(FileNotFoundError, match="v1"):
        vp_oxphos.predict([ASPIRIN], version="v999")


def test_fixture_satisfies_the_dataset_contract():
    from vp_core import dataset as core_dataset

    fixture = data.example()
    assert core_dataset.validate_table(fixture) == []
    assert set(fixture["label"]) == {0, 1}


def test_majority_vote_drops_a_compound_the_assay_called_both_ways():
    """A tie carries no clean label, so it must not be broken arbitrarily."""
    import pandas as pd

    records = pd.DataFrame(
        {
            "cid": [1, 2, 3, 3],
            "outcome": ["Active", "Inactive", "Active", "Inactive"],
            "potency_um": [4.0, float("nan"), 8.0, float("nan")],
        }
    )
    structures = pd.DataFrame(
        {"cid": [1, 2, 3], "smiles_raw": [ASPIRIN, CAFFEINE, "CCO"]}
    )
    table = data._to_compounds(records, structures)

    assert len(table) == 2
    assert set(table["label"]) == {0, 1}
    assert "CCO" not in set(table["smiles"])


def test_target_registry_is_uniform():
    assert vp_oxphos.all_names() == ["MMP"]
    assert vp_oxphos.get_target("mmp").pubchem_aid == 720635
    with pytest.raises(KeyError):
        vp_oxphos.get_target("nope")


def test_evaluation_harness_runs_on_the_fixture():
    """A smoke run must never be recordable as a released result."""
    from vp_oxphos.evaluate import evaluate_version

    record = evaluate_version(vp_oxphos.current_version(), use_example=True, write=False)
    assert record["dataset"] == "example fixture"
    with pytest.raises(ValueError, match="refusing to record"):
        evaluate_version(vp_oxphos.current_version(), use_example=True, write=True)
