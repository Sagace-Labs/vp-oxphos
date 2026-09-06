"""The endpoints this pathway predicts.

Frozen so the data sources are auditable.

Both assays were verified live against PubChem on 2026-09-01, and screen the
same library:

    AID 720635   qHTS assay for small molecule disruptors of the mitochondrial
                 membrane potential                       10486 substances
    AID 720634   the same, cell viability                 10486 substances
"""

from __future__ import annotations

from dataclasses import dataclass

__all__ = ["CYTOTOX", "TARGET", "TARGETS", "Endpoint", "all_names", "get"]


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

CYTOTOX = Endpoint(
    name="VIABILITY",
    pathway="OXPHOS / mitochondrial integrity",
    mie="loss of cell viability, which collapses the membrane potential as a consequence",
    pubchem_aid=720634,
    assay_name=(
        "qHTS assay for small molecule disruptors of the mitochondrial membrane "
        "potential - cell viability"
    ),
)

TARGETS: dict[str, Endpoint] = {"MMP": TARGET, "VIABILITY": CYTOTOX}


def get(name: str) -> Endpoint:
    try:
        return TARGETS[name.upper()]
    except KeyError:
        raise KeyError(f"unknown target {name!r}; known: {sorted(TARGETS)}") from None


def all_names() -> list[str]:
    return sorted(TARGETS)
