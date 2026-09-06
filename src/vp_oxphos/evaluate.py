"""Measure a version under an evaluation protocol.

    python -m vp_oxphos.evaluate --version v1

Writes ``versions/<version>/metrics.json`` and regenerates ``CARD.md``. The
protocol named in the version's manifest supplies the split, the fold sizes,
the seed set and the metric list.

Each seed refits the model on its own training fold. These are therefore *not* the
shipped weights, which are fit on everything; see ``train``.

A seed that yields an empty fold raises rather than being skipped. Dropping a
seed silently would change the estimator without changing the protocol id.
"""

from __future__ import annotations

import argparse
import json
import sys
from datetime import date
from pathlib import Path

import numpy as np

from vp_oxphos import data as oxphos_data
from vp_oxphos import model as oxphos_model

__all__ = ["evaluate_version", "main"]


def evaluate_version(
    version: str, *, use_example: bool = False, write: bool = True
) -> dict:
    """Run the version's protocol and return the metrics record."""
    import vp_oxphos
    from vp_core import card, dataset, fingerprints, protocols
    from vp_core import metrics as metrics_mod

    resolved = vp_oxphos.get(version)
    protocol = protocols.get(resolved.protocol_id)

    table = oxphos_data.example() if use_example else oxphos_data.load()
    smiles = table["smiles"].tolist()
    y = table["label"].to_numpy(dtype=int)
    X = fingerprints.featurize(smiles, oxphos_model.FEATURES)

    per_seed: list[dict[str, float]] = []
    folds: list[dict[str, int]] = []
    for seed in protocol.seeds:
        train_idx, val_idx, test_idx = protocol.split_indices(smiles, seed)
        fitted = oxphos_model.fit(
            X[train_idx],
            y[train_idx],
            X[val_idx],
            y[val_idx],
            seed=seed,
        )
        proba = oxphos_model.predict(fitted, [smiles[i] for i in test_idx])
        scored = metrics_mod.binary_metrics(y[test_idx], proba)
        per_seed.append(scored)
        folds.append(
            {
                "seed": int(seed),
                "train": len(train_idx),
                "val": len(val_idx),
                "test": len(test_idx),
            }
        )
        print(
            f"  seed {seed}: test AUC {scored['auc_roc']:.4f} "
            f"(n_test={len(test_idx)})",
            file=sys.stderr,
        )

    record = {
        "pathway": "oxphos",
        "version": resolved.name,
        "protocol_id": resolved.protocol_id,
        "dataset_sha256": dataset.dataset_hash(table),
        "dataset": "example fixture" if use_example else "full",
        "evaluated": date.today().isoformat(),
        "n": {
            "total": len(table),
            "train": folds[0]["train"],
            "val": folds[0]["val"],
            "test": folds[0]["test"],
        },
        "folds": folds,
        "test": metrics_mod.aggregate(per_seed, protocol.metrics),
    }

    if write:
        if use_example:
            raise ValueError(
                "refusing to record fixture metrics as a released result — "
                "--example is for smoke-checking the harness only"
            )
        path = Path(resolved.directory) / "metrics.json"
        path.write_text(json.dumps(record, indent=2) + "\n", encoding="utf-8")
        card.write_card(resolved.directory)
        auc = record["test"]["auc_roc"]
        print(
            f"wrote {path}\n"
            f"  {resolved.protocol_id}: AUC {auc['mean']:.4f} +/- {auc['std']:.4f} "
            f"over {len(protocol.seeds)} seeds",
            file=sys.stderr,
        )
    return record


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="python -m vp_oxphos.evaluate")
    parser.add_argument("--version", default=None, help="default: the newest version")
    parser.add_argument(
        "--example",
        action="store_true",
        help="run on the committed fixture without writing (smoke check only)",
    )
    args = parser.parse_args(argv)

    import vp_oxphos

    version = args.version or vp_oxphos.current_version()
    record = evaluate_version(version, use_example=args.example, write=not args.example)
    if args.example:
        auc = record["test"]["auc_roc"]["mean"]
        print(f"fixture smoke AUC {auc if np.isfinite(auc) else float('nan'):.4f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
