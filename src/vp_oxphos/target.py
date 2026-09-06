"""The endpoint this pathway predicts.

Frozen so the data source is auditable. Unlike a single-protein pathway there
is no one target here: collapsing the mitochondrial membrane potential is the
shared consequence of electron-transport-chain inhibition, protonophoric
uncoupling, ATP-synthase inhibition and permeability-transition opening. The
assay is therefore the definition of the label, and it is recorded verbatim.

The assay was verified live against PubChem on 2026-09-01:

    AID 720635   qHTS assay for small molecule disruptors of the mitochondrial
                 membrane potential   Confirmatory   10486 substances
"""

from __future__ import annotations

from dataclasses import dataclass

__all__ = ["TARGET", "TARGETS", "Endpoint", "all_names", "get"]


@dataclass(frozen=True)
class Endpoint:
    """One molecular initiating event, identified by the assay that reads it out."""

    name: str
    pathway: str
    mie: str
    pubchem_aid: int
    assay_name: str


TARGET = Endpoint(
    name="MMP",
    pathway="OXPHOS / mitochondrial integrity",
    mie=(
        "loss of the inner-membrane proton gradient through electron-transport-chain "
        "inhibition, protonophoric uncoupling or permeability-transition opening"
    ),
    pubchem_aid=720635,
    assay_name=(
        "qHTS assay for small molecule disruptors of the mitochondrial membrane potential"
    ),
)

TARGETS: dict[str, Endpoint] = {"MMP": TARGET}


def get(name: str) -> Endpoint:
    try:
        return TARGETS[name.upper()]
    except KeyError:
        raise KeyError(f"unknown target {name!r}; known: {sorted(TARGETS)}") from None


def all_names() -> list[str]:
    return sorted(TARGETS)
