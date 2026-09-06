"""Structural alerts for mitochondrial toxicity."""

from __future__ import annotations

from functools import lru_cache

import numpy as np

__all__ = ["ALERTS", "column_names", "counts", "width"]

ALERTS: tuple[tuple[str, str], ...] = (
    ("c1cc2C(=O)c3c(C(=O)c2c(O)c1)c(O)ccc3", "danthron"),
    ("CNc1ccc(C(F))cc1", "4-(fluoromethyl)-N-methylaniline"),
    ("C(=O)(c1cc(O)ccc1)c1ccc(cc1)", "4-hydroxybenzophenone"),
    ("O=C1c2ccccc2C(=O)c2ccccc21", "anthracene-9,10-dione"),
    ("c1(Cl)cccc(c1O)Cl", "2,6-dichlorophenol"),
    ("C1=CNC=CC1", "dihydropyridine"),
    ("[N+](=O)([O-])c1ccc(c(c1)C)O", "2-methyl-4-nitrophenol"),
    ("c1(Cl)c(Cl)ccc(c1)O", "3,4-dichlorophenol"),
    ("c1(N)ccc(cc1)N=N", "4-diazenylaniline"),
    ("O=C1C=CCO1", "2,5-dihydrofuran-2-one"),
    ("c1ccc2cc3ccccc3cc2c1", "anthracene"),
    ("O=C(Nc1ccccc1)c1ccccc1", "N-phenylbenzamide"),
    ("CCOC(=O)c1ccc(O)cc1", "ethyl 4-hydroxybenzoate"),
    ("C1=CCOC1", "2,5-dihydrofuran"),
    ("C(=O)(Nc1ccc(N)cc1)", "N-(4-aminophenyl)formamide"),
    ("ClC(Cl)Cc1ccccc1", "(2,2-dichloroethyl)benzene"),
    ("c1c(cccc1)C(=C)c1ccccc1", "(1-phenylethenyl)benzene"),
)


def column_names() -> list[str]:
    return [name for _, name in ALERTS]


def width() -> int:
    return len(ALERTS)


@lru_cache(maxsize=1)
def _patterns():
    from rdkit import Chem, RDLogger

    RDLogger.DisableLog("rdApp.*")
    compiled = [Chem.MolFromSmarts(smarts) for smarts, _ in ALERTS]
    missing = [name for patt, (_, name) in zip(compiled, ALERTS, strict=True) if patt is None]
    if missing:
        raise ValueError(f"unparseable alert SMARTS: {missing}")
    return compiled


def counts(smiles: list[str]) -> np.ndarray:
    """Match count per alert, one row per input. Unparseable inputs are zero."""
    from rdkit import Chem

    patterns = _patterns()
    X = np.zeros((len(smiles), len(patterns)), dtype=np.float32)
    for i, smi in enumerate(smiles):
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            continue
        for j, patt in enumerate(patterns):
            X[i, j] = len(mol.GetSubstructMatches(patt))
    return X
