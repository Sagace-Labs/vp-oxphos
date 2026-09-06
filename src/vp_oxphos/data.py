"""The OXPHOS dataset: obtain, standardise, verify.

The shipped table is the *parsed* result — one row per compound, keyed on the
standardised InChIKey — and that is what the manifest hash covers. ``fetch``
rebuilds it from PubChem and compares hashes, so its real job is drift
detection: a differing hash means the upstream assay record changed, which is a
scientific event worth a new version.

Pipeline: pull the concise BioAssay table for AID 720635, keep the rows the
assay called Active or Inactive, resolve each compound identifier to a SMILES
string through the compound property endpoint, standardise, then collapse to
one row per compound by majority vote across its assay records.

Run as a command:

    python -m vp_oxphos.data fetch --verify
    python -m vp_oxphos.data verify
"""

from __future__ import annotations

import argparse
import io
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

from vp_oxphos.target import TARGET

__all__ = [
    "DATA_DIR",
    "EXAMPLE_PATH",
    "TABLE_PATH",
    "build_example",
    "example",
    "fetch",
    "load",
    "verify",
]

_PACKAGE_ROOT = Path(__file__).resolve().parents[2]
DATA_DIR = _PACKAGE_ROOT / "data"
TABLE_PATH = DATA_DIR / "oxphos_tox21.parquet"
EXAMPLE_PATH = DATA_DIR / "example" / "oxphos_example.parquet"

PUG_BASE = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
_MIN_INTERVAL = 0.25  # PubChem asks for no more than five requests a second
_CID_BATCH = 200  # the property endpoint carries its identifiers in the URL

# The two calls that carry a label. A qHTS curve the assay could not call
# either way is not a negative, so "Inconclusive" is dropped rather than zeroed.
_ACTIVE = "Active"
_INACTIVE = "Inactive"


def _missing(path: Path) -> str:
    return (
        f"{path} is not present. The dataset is not included in the wheel; run "
        "`python -m vp_oxphos.data fetch` to rebuild it from PubChem."
    )


def load() -> pd.DataFrame:
    """The full standardised dataset. Raises with guidance when absent."""
    from vp_core import dataset

    if not TABLE_PATH.exists():
        raise FileNotFoundError(_missing(TABLE_PATH))
    return dataset.read_table(TABLE_PATH)


def example() -> pd.DataFrame:
    """The committed test fixture — small, stratified, always available."""
    from vp_core import dataset

    if not EXAMPLE_PATH.exists():
        raise FileNotFoundError(_missing(EXAMPLE_PATH))
    return dataset.read_table(EXAMPLE_PATH)


# ---------------------------------------------------------------------------
# Rebuild from source
# ---------------------------------------------------------------------------


def _assay_records() -> pd.DataFrame:
    """The concise BioAssay table for the endpoint's AID, as PubChem serves it."""
    import requests

    url = f"{PUG_BASE}/assay/aid/{TARGET.pubchem_aid}/concise/CSV"
    response = requests.get(url, headers={"Accept": "text/csv"}, timeout=300)
    response.raise_for_status()
    raw = pd.read_csv(io.BytesIO(response.content), dtype=str, low_memory=False)

    for column in ("CID", "Activity Outcome"):
        if column not in raw.columns:
            raise RuntimeError(
                f"AID {TARGET.pubchem_aid} returned no {column!r} column; the concise "
                f"format changed upstream. Columns: {raw.columns.tolist()}"
            )

    df = pd.DataFrame(
        {
            "cid": pd.to_numeric(raw["CID"], errors="coerce"),
            "outcome": raw["Activity Outcome"].astype(str).str.strip(),
            "potency_um": pd.to_numeric(raw.get("Activity Value [uM]"), errors="coerce"),
        }
    )
    df = df[df["outcome"].isin([_ACTIVE, _INACTIVE])].dropna(subset=["cid"])
    df["cid"] = df["cid"].astype(int)
    print(
        f"  {len(df)} labelled records "
        f"({int((df['outcome'] == _ACTIVE).sum())} active)",
        file=sys.stderr,
    )
    return df.reset_index(drop=True)


def _smiles_for(cids: list[int]) -> pd.DataFrame:
    """``(cid, smiles_raw)`` from the compound property endpoint, in batches.

    Isomeric SMILES with the connectivity form as a fallback, since PubChem
    renamed both properties and older mirrors still answer under the old names.
    """
    import requests

    session = requests.Session()
    frames: list[pd.DataFrame] = []
    unique = sorted(set(cids))
    for start in range(0, len(unique), _CID_BATCH):
        batch = unique[start : start + _CID_BATCH]
        ids = ",".join(str(c) for c in batch)
        url = f"{PUG_BASE}/compound/cid/{ids}/property/SMILES,ConnectivitySMILES/CSV"
        response = session.get(url, headers={"Accept": "text/csv"}, timeout=120)
        response.raise_for_status()
        page = pd.read_csv(io.BytesIO(response.content))
        primary = "SMILES" if "SMILES" in page.columns else "IsomericSMILES"
        fallback = (
            "ConnectivitySMILES"
            if "ConnectivitySMILES" in page.columns
            else "CanonicalSMILES"
        )
        frames.append(
            pd.DataFrame(
                {
                    "cid": pd.to_numeric(page["CID"], errors="coerce"),
                    "smiles_raw": page[primary].fillna(page.get(fallback)),
                }
            )
        )
        print(
            f"  resolved {min(start + _CID_BATCH, len(unique))} / {len(unique)}",
            file=sys.stderr,
        )
        time.sleep(_MIN_INTERVAL)

    if not frames:
        raise RuntimeError("PubChem returned no compound structures at all.")
    out = pd.concat(frames, ignore_index=True).dropna(subset=["cid", "smiles_raw"])
    out["cid"] = out["cid"].astype(int)
    return out.drop_duplicates(subset=["cid"]).reset_index(drop=True)


def _to_compounds(records: pd.DataFrame, structures: pd.DataFrame) -> pd.DataFrame:
    """Collapse assay records to one row per standardised compound.

    Salt and charge variants of the same compound carry separate identifiers
    upstream, so several records can land on one InChIKey. The label is the
    majority call across them; a tie is dropped rather than broken arbitrarily,
    because a compound the assay called both ways carries no clean label.
    """
    from rdkit import Chem, RDLogger

    from vp_core.standardise import standardise_many

    RDLogger.DisableLog("rdApp.*")

    df = records.merge(structures, on="cid", how="inner")
    smiles, keys = standardise_many(df["smiles_raw"].astype(str).tolist())
    df["smiles"] = smiles
    df["inchikey"] = keys
    df = df.dropna(subset=["smiles", "inchikey"])

    # Standardisation can emit a SMILES that RDKit will not read back, and a
    # featuriser turns one of those into an all-zero row rather than an error.
    # Requiring the round trip keeps a silently empty compound out of the table.
    df = df[[Chem.MolFromSmiles(s) is not None for s in df["smiles"]]]

    df["active"] = (df["outcome"] == _ACTIVE).astype(int)

    rows: list[dict] = []
    for inchikey, group in df.groupby("inchikey"):
        active_frac = float(group["active"].mean())
        if active_frac == 0.5:
            continue
        reported = group["potency_um"].dropna()
        rows.append(
            {
                "inchikey": inchikey,
                "smiles": group["smiles"].iloc[0],
                "label": int(active_frac > 0.5),
                # Provenance. The qHTS potency is reported for a minority of
                # records, which is why the label is the assay's own call and
                # not a threshold on this column.
                "potency_um": float(np.median(reported)) if len(reported) else float("nan"),
                "n_calls": len(group),
                "active_frac": round(active_frac, 6),
            }
        )
    return pd.DataFrame(rows).sort_values("inchikey").reset_index(drop=True)


def fetch(*, write: bool = True) -> pd.DataFrame:
    """Rebuild the dataset from PubChem. Returns the standardised table."""
    from vp_core import dataset

    print(f"fetching AID {TARGET.pubchem_aid} from PubChem", file=sys.stderr)
    records = _assay_records()
    structures = _smiles_for(records["cid"].tolist())
    table = _to_compounds(records, structures)

    problems = dataset.validate_table(table)
    if problems:
        raise ValueError(f"rebuilt table is invalid: {'; '.join(problems)}")
    if write:
        dataset.write_table(table, TABLE_PATH)
    print(
        f"{len(table)} compounds, {int(table['label'].sum())} positive "
        f"({table['label'].mean():.1%}), sha256 {dataset.dataset_hash(table)}",
        file=sys.stderr,
    )
    return table


def verify(version: str | None = None) -> dict:
    """Compare the on-disk table against the hash a released version recorded."""
    import vp_oxphos
    from vp_core import dataset

    resolved = vp_oxphos.get(version)
    declared = resolved.manifest.get("dataset", {}).get("sha256")
    actual = dataset.dataset_hash(load())
    return {
        "version": resolved.name,
        "declared": declared,
        "actual": actual,
        "match": declared == actual,
    }


def build_example(n: int = 200, seed: int = 0) -> pd.DataFrame:
    """Regenerate the committed fixture from the full table."""
    from vp_core import dataset

    sample = dataset.stratified_example(load(), n=n, seed=seed)
    dataset.write_table(sample, EXAMPLE_PATH)
    return sample


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="python -m vp_oxphos.data")
    sub = parser.add_subparsers(dest="command", required=True)

    fetch_cmd = sub.add_parser("fetch", help="rebuild the dataset from PubChem")
    fetch_cmd.add_argument(
        "--verify",
        action="store_true",
        help="after fetching, compare the hash against the current version",
    )
    sub.add_parser("verify", help="check the on-disk table against the recorded hash")
    example_cmd = sub.add_parser("build-example", help="regenerate the test fixture")
    example_cmd.add_argument("--n", type=int, default=200)

    args = parser.parse_args(argv)

    if args.command == "fetch":
        fetch()
        if not args.verify:
            return 0
    if args.command == "build-example":
        sample = build_example(n=args.n)
        print(f"wrote {len(sample)} rows to {EXAMPLE_PATH}")
        return 0

    report = verify()
    if report["match"]:
        print(f"dataset matches {report['version']}: {report['actual']}")
        return 0
    print(
        f"DRIFT: {report['version']} recorded {report['declared']}\n"
        f"       the current table hashes to {report['actual']}\n"
        "The upstream source has changed. Record it as a new version rather than "
        "overwriting the released one.",
        file=sys.stderr,
    )
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
