#!/usr/bin/env python3
"""Rebuild the tile training round and resume Cellpose training."""

from __future__ import annotations

import argparse
import subprocess
from pathlib import Path

PROJECT = Path(__file__).resolve().parent.parent
DEFAULT_BUNDLE = PROJECT / "output" / "tile_annotation_bundle_v1"
DEFAULT_ROUND = PROJECT / "output" / "tile_training_round_v1"
DEFAULT_MODEL_NAME = "desmognathus_tile_round1"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tile-bundle", type=Path, default=DEFAULT_BUNDLE)
    parser.add_argument("--round-dir", type=Path, default=DEFAULT_ROUND)
    parser.add_argument(
        "--label-dir",
        type=Path,
        action="append",
        default=[],
        help="Optional directory containing corrected tile images plus sibling *_seg.npy files",
    )
    parser.add_argument("--model-name", default=DEFAULT_MODEL_NAME)
    parser.add_argument("--base-model", default="cpsam", help="Used when no trained model exists yet")
    parser.add_argument("--test-fraction", type=float, default=0.2)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--min-masks", type=int, default=3)
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
        "python3",
        str(PROJECT / "scripts" / "prepare_tile_training_round.py"),
        "--tile-bundle",
        str(args.tile_bundle),
        "--output-dir",
        str(args.round_dir),
        "--test-fraction",
        str(args.test_fraction),
        "--seed",
        str(args.seed),
        "--min-masks",
        str(args.min_masks),
    ]
    for label_dir in args.label_dir:
        prepare_cmd.extend(["--label-dir", str(label_dir)])
    run(prepare_cmd, args.dry_run)

    pretrained_model = choose_pretrained_model(args.round_dir, args.model_name, args.base_model)
    train_cmd = [
        "python3",
        "-m",
        "cellpose",
        "--train",
        "--dir",
        str(args.round_dir / "train"),
        "--test_dir",
        str(args.round_dir / "test"),
        "--mask_filter",
        "_masks",
        "--pretrained_model",
        pretrained_model,
        "--model_name_out",
        args.model_name,
        "--n_epochs",
        str(args.n_epochs),
        "--save_every",
        str(args.save_every),
        "--save_each",
        "--verbose",
    ]
    if args.use_gpu:
        train_cmd.append("--use_gpu")
    run(train_cmd, args.dry_run)


if __name__ == "__main__":
    main()
