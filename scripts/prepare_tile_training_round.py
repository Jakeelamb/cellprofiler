#!/usr/bin/env python3
"""Build a Cellpose training round from manually annotated tile bundles."""

from __future__ import annotations

import argparse
import json
import shutil
from pathlib import Path

import numpy as np
import pandas as pd

PROJECT = Path(__file__).resolve().parent.parent

from prepare_cellpose_training_round import labeled_u16, png_bytes, resolve_corrected_mask  # noqa: E402

DEFAULT_BUNDLE = PROJECT / "output" / "tile_annotation_bundle_v1"
DEFAULT_OUTPUT = PROJECT / "output" / "tile_training_round_v1"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tile-bundle", type=Path, default=DEFAULT_BUNDLE)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument(
        "--label-dir",
        type=Path,
        action="append",
        default=[],
        help="Optional directory containing corrected tile images plus sibling *_seg.npy files",
    )
    parser.add_argument("--test-fraction", type=float, default=0.2)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--min-masks", type=int, default=3)
    return parser.parse_args()


def require_exists(path: Path) -> None:
    if not path.exists():
        raise FileNotFoundError(path)


def ensure_clean_dir(path: Path) -> None:
    if path.exists():
        shutil.rmtree(path)
    path.mkdir(parents=True, exist_ok=True)


def clear_dir_except(path: Path, keep_names: set[str] | None = None) -> None:
    keep_names = keep_names or set()
    path.mkdir(parents=True, exist_ok=True)
    for child in path.iterdir():
        if child.name in keep_names:
            continue
        if child.is_dir():
            shutil.rmtree(child)
        else:
            child.unlink()


def find_annotation_source(image_path: Path, label_dirs: list[Path]) -> tuple[Path, str, np.ndarray] | None:
    candidates = []
    for label_dir in label_dirs:
        candidates.append((label_dir / image_path.name).resolve())
    candidates.append(image_path.resolve())

    for candidate in candidates:
        if not candidate.exists():
            continue
        mask_source_kind, corrected = resolve_corrected_mask(candidate, candidate)
        if corrected is None:
            continue
        corrected = labeled_u16(corrected)
        return candidate, mask_source_kind, corrected
    return None


def split_reviewed(df: pd.DataFrame, test_fraction: float, seed: int) -> pd.DataFrame:
    if df.empty:
        out = df.copy()
        out["split"] = pd.Series(dtype=object)
        return out

    rng = np.random.default_rng(seed)
    frame = df.sample(frac=1.0, random_state=seed).reset_index(drop=True)
    split_map: dict[int, str] = {}
    for _, grp in frame.groupby("species", sort=False):
        idx = list(grp.index)
        rng.shuffle(idx)
        n = len(idx)
        if n <= 1:
            test_n = 0
        else:
            test_n = max(1, int(round(n * test_fraction)))
            test_n = min(test_n, n - 1)
        test_idx = set(idx[:test_n])
        for row_idx in idx:
            split_map[row_idx] = "test" if row_idx in test_idx else "train"
    frame["split"] = frame.index.map(split_map)
    return frame


def load_annotated_tiles(bundle_dir: Path, min_masks: int, label_dirs: list[Path]) -> pd.DataFrame:
    require_exists(bundle_dir)
    manifest_path = bundle_dir / "manifest.csv"
    require_exists(manifest_path)
    manifest = pd.read_csv(manifest_path)
    manifest["image_path_abs"] = manifest["image_path"].map(lambda p: str((bundle_dir / p).resolve()))

    rows = []
    for _, row in manifest.iterrows():
        image_path = Path(str(row["image_path_abs"]))
        annotation = find_annotation_source(image_path, label_dirs)
        if annotation is None:
            continue
        annotation_image_path, mask_source_kind, corrected = annotation
        n_masks = int(corrected.max(initial=0))
        if n_masks < min_masks:
            continue
        rows.append(
            {
                **row.to_dict(),
                "image_path_abs": str(annotation_image_path),
                "source_bundle_image_path_abs": str(image_path),
                "mask_source_kind": mask_source_kind,
                "n_masks": n_masks,
            }
        )
    return pd.DataFrame(rows)


def copy_training_pairs(reviewed: pd.DataFrame, output_dir: Path) -> pd.DataFrame:
    rows = []
    for _, row in reviewed.iterrows():
        split_dir = output_dir / row["split"]
        split_dir.mkdir(parents=True, exist_ok=True)
        src_image = Path(str(row["image_path_abs"]))
        dst_image = split_dir / src_image.name
        dst_mask = split_dir / f"{src_image.stem}_masks.png"
        shutil.copy2(src_image, dst_image)
        _, corrected = resolve_corrected_mask(src_image, src_image)
        if corrected is None:
            raise ValueError(f"missing corrected mask for {src_image}")
        dst_mask.write_bytes(png_bytes(labeled_u16(corrected)))
        rows.append(
            {
                "split": row["split"],
                "bundle_index": int(row["bundle_index"]),
                "filename": row["filename"],
                "species": row["species"],
                "decision": row["decision"],
                "tile_name": row["tile_name"],
                "n_masks": int(row["n_masks"]),
                "mask_source_kind": row["mask_source_kind"],
                "annotation_image_path": str(src_image),
                "image_path": str(dst_image.relative_to(output_dir)),
                "mask_path": str(dst_mask.relative_to(output_dir)),
            }
        )
    return pd.DataFrame(rows)


def write_readme(output_dir: Path, summary: dict[str, int | str]) -> None:
    lines = [
        "# Tile Training Round",
        "",
        "This directory contains manually annotated raw tiles exported from the tile annotation bundle.",
        "",
        "## Contents",
        "",
        "- `train/`: training images and `_masks.png` labels",
        "- `test/`: held-out images and `_masks.png` labels",
        "- `training_manifest.csv`: per-tile training metadata",
        "",
        "## Suggested Training Command",
        "",
        "```bash",
        f"python3 -m cellpose --train --dir {output_dir / 'train'} --test_dir {output_dir / 'test'} --mask_filter _masks --pretrained_model cpsam --model_name_out desmognathus_tile_round1 --use_gpu",
        "```",
        "",
        "Drop `--use_gpu` if you want CPU training.",
        "",
        "## Summary",
        "",
    ]
    for key, value in sorted(summary.items()):
        lines.append(f"- `{key}`: `{value}`")
    (output_dir / "README.md").write_text("\n".join(lines) + "\n")


def main() -> None:
    args = parse_args()
    label_dirs = [path.resolve() for path in args.label_dir]
    for path in label_dirs:
        require_exists(path)
    reviewed = load_annotated_tiles(args.tile_bundle, args.min_masks, label_dirs)
    if reviewed.empty:
        raise SystemExit("No annotated tiles found with usable Cellpose labels")

    reviewed = split_reviewed(reviewed, args.test_fraction, args.seed)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    clear_dir_except(args.output_dir / "train", keep_names={"models"})
    clear_dir_except(args.output_dir / "test")
    training_manifest = copy_training_pairs(reviewed, args.output_dir)
    training_manifest.to_csv(args.output_dir / "training_manifest.csv", index=False)

    summary = {
        "tile_bundle": str(args.tile_bundle.resolve()),
        "label_dirs": [str(path) for path in label_dirs],
        "n_annotated_tiles": int(len(reviewed)),
        "n_train": int((reviewed["split"] == "train").sum()),
        "n_test": int((reviewed["split"] == "test").sum()),
        "n_species": int(reviewed["species"].nunique()),
        "min_masks": int(args.min_masks),
        "median_masks_per_tile": int(reviewed["n_masks"].median()),
        "max_masks_per_tile": int(reviewed["n_masks"].max()),
    }
    (args.output_dir / "summary.json").write_text(json.dumps(summary, indent=2))
    write_readme(args.output_dir, summary)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
