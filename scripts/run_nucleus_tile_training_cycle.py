#!/usr/bin/env python3
"""Rebuild the nucleus tile training round and resume Cellpose training."""

from __future__ import annotations

import argparse
import os
import subprocess
import sys
from pathlib import Path

PROJECT = Path(__file__).resolve().parent.parent
DEFAULT_SOURCE_ROUND = PROJECT / "output" / "tile_training_round_v1"
DEFAULT_ROUND = PROJECT / "output" / "nucleus_tile_training_round_v1"
DEFAULT_ARTIFACT_ROOT = PROJECT / "output" / "runs" / "threshold_tuned_v1" / "nucleus_iod" / "artifacts"
DEFAULT_MODEL_NAME = "desmognathus_nucleus_tile_round1"
DEFAULT_CELLPOSE_PYTHON = os.environ.get("CELLPOSE_PYTHON", "")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-round", type=Path, default=DEFAULT_SOURCE_ROUND)
    parser.add_argument("--round-dir", type=Path, default=DEFAULT_ROUND)
    parser.add_argument("--artifact-root", type=Path, default=DEFAULT_ARTIFACT_ROOT)
    parser.add_argument(
        "--cellpose-python",
        default=DEFAULT_CELLPOSE_PYTHON,
        help="Python interpreter for Cellpose training helpers; defaults to current uv env when unset",
    )
    parser.add_argument(
        "--label-dir",
        type=Path,
        action="append",
        default=[],
        help="Optional directory containing corrected tile images plus sibling *_seg.npy files",
    )
    parser.add_argument(
        "--corrected-only",
        action="store_true",
        help="Use only manually corrected labels from --label-dir and skip unlabeled tiles",
    )
    parser.add_argument("--model-name", default=DEFAULT_MODEL_NAME)
    parser.add_argument("--base-model", default="nuclei", help="Used when no trained model exists yet")
    parser.add_argument("--min-masks", type=int, default=3)
    parser.add_argument("--train-batch-size", type=int, default=1)
    parser.add_argument("--bsize", type=int, default=256)
    parser.add_argument("--n-epochs", type=int, default=100)
    parser.add_argument("--save-every", type=int, default=10)
    parser.add_argument("--use-gpu", action="store_true")
    parser.add_argument("--dry-run", action="store_true")
    return parser.parse_args()


def choose_pretrained_model(round_dir: Path, model_name: str, base_model: str) -> str:
    model_path = round_dir / "train" / "models" / model_name
    if model_path.exists():
        return str(model_path)
    return base_model


def run(cmd: list[str], dry_run: bool) -> None:
    print("$", " ".join(cmd))
    if dry_run:
        return
    subprocess.run(cmd, check=True)


def main() -> None:
    args = parse_args()

    prepare_cmd = [
        sys.executable,
        str(PROJECT / "scripts" / "prepare_nucleus_tile_training_round.py"),
        "--source-round",
        str(args.source_round),
        "--output-dir",
        str(args.round_dir),
        "--artifact-root",
        str(args.artifact_root),
        "--min-masks",
        str(args.min_masks),
    ]
    for label_dir in args.label_dir:
        prepare_cmd.extend(["--label-dir", str(label_dir)])
    if args.corrected_only:
        prepare_cmd.append("--corrected-only")
    run(prepare_cmd, args.dry_run)

    pretrained_model = choose_pretrained_model(args.round_dir, args.model_name, args.base_model)
    train_cmd = [
        args.cellpose_python if args.cellpose_python else sys.executable,
        str(PROJECT / "scripts" / "train_cellpose_with_bsize.py"),
        "--dir",
        str(args.round_dir / "train"),
        "--test-dir",
        str(args.round_dir / "test"),
        "--mask-filter",
        "_masks",
        "--pretrained-model",
        pretrained_model,
        "--model-name-out",
        args.model_name,
        "--train-batch-size",
        str(args.train_batch_size),
        "--bsize",
        str(args.bsize),
        "--n-epochs",
        str(args.n_epochs),
        "--save-every",
        str(args.save_every),
        "--save-each",
        "--verbose",
    ]
    if args.use_gpu:
        train_cmd.append("--keep-bfloat16")
        train_cmd.append("--use-gpu")
    run(train_cmd, args.dry_run)


if __name__ == "__main__":
    main()
