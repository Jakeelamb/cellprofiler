#!/usr/bin/env python3
"""Build a nucleus tile training round from the existing tile split and nucleus artifacts."""

from __future__ import annotations

import argparse
import json
import shutil
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd
from roifile import ImagejRoi
from scipy import ndimage as ndi
from skimage.draw import polygon
import tifffile

PROJECT = Path(__file__).resolve().parent.parent

from prepare_cellpose_training_round import labeled_u16, png_bytes, resolve_corrected_mask  # noqa: E402

DEFAULT_SOURCE_ROUND = PROJECT / "output" / "tile_training_round_v1"
DEFAULT_OUTPUT = PROJECT / "output" / "nucleus_tile_training_round_v1"
DEFAULT_ARTIFACT_ROOT = PROJECT / "output" / "runs" / "threshold_tuned_v1" / "nucleus_iod" / "artifacts"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-round", type=Path, default=DEFAULT_SOURCE_ROUND)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--artifact-root", type=Path, default=DEFAULT_ARTIFACT_ROOT)
    parser.add_argument(
        "--label-dir",
        type=Path,
        action="append",
        default=[],
        help="Optional directory containing corrected tile images plus sibling *_seg.npy files",
    )
    parser.add_argument("--min-masks", type=int, default=3)
    parser.add_argument(
        "--corrected-only",
        action="store_true",
        help="Use only manually corrected labels from --label-dir and skip unlabeled tiles",
    )
    return parser.parse_args()


def require_exists(path: Path) -> None:
    if not path.exists():
        raise FileNotFoundError(path)


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


def artifact_dir(artifact_root: Path, filename: str) -> Path:
    return artifact_root / Path(str(filename)).stem


def roi_zip_path(artifact_root: Path, filename: str, tile_name: str) -> Path:
    return artifact_dir(artifact_root, filename) / "macro_artifacts" / f"{Path(tile_name).stem}__rois.zip"


def binary_mask_path(artifact_root: Path, filename: str, tile_name: str) -> Path:
    return artifact_dir(artifact_root, filename) / "macro_artifacts" / f"{Path(tile_name).stem}__mask.tif"


def render_roi_labels(path: Path, shape: tuple[int, int]) -> np.ndarray:
    labels = np.zeros(shape, dtype=np.uint16)
    next_label = 1
    with zipfile.ZipFile(path) as handle:
        for name in sorted(handle.namelist()):
            try:
                roi = ImagejRoi.frombytes(handle.read(name))
            except Exception:
                continue
            coords = np.asarray(roi.coordinates(), dtype=float)
            if coords.size == 0:
                continue
            rr, cc = polygon(coords[:, 1], coords[:, 0], shape=shape)
            if len(rr) == 0:
                continue
            labels[rr, cc] = next_label
            next_label += 1
    return labels


def seed_labels(artifact_root: Path, filename: str, tile_name: str, shape: tuple[int, int]) -> tuple[str, np.ndarray]:
    roi_path = roi_zip_path(artifact_root, filename, tile_name)
    if roi_path.exists():
        labels = render_roi_labels(roi_path, shape)
        if labels.max(initial=0) > 0:
            return "roi_zip", labels
    mask_path = binary_mask_path(artifact_root, filename, tile_name)
    require_exists(mask_path)
    binary = tifffile.imread(mask_path) > 0
    labels, _ = ndi.label(binary)
    return "binary_components", labeled_u16(labels)


def corrected_label_source(image_path: Path, label_dirs: list[Path]) -> tuple[Path, str, np.ndarray] | None:
    for label_dir in label_dirs:
        candidate = (label_dir / image_path.name).resolve()
        if not candidate.exists():
            continue
        mask_source_kind, corrected = resolve_corrected_mask(candidate, candidate)
        if corrected is None:
            continue
        return candidate, mask_source_kind, labeled_u16(corrected)
    return None


def load_rows(source_round: Path) -> pd.DataFrame:
    manifest_path = source_round / "training_manifest.csv"
    require_exists(manifest_path)
    frame = pd.read_csv(manifest_path)
    required = {"split", "filename", "species", "tile_name", "image_path"}
    missing = required - set(frame.columns)
    if missing:
        raise ValueError(f"missing columns in {manifest_path}: {sorted(missing)}")
    frame["source_image_abs"] = frame["image_path"].map(lambda p: str((source_round / str(p)).resolve()))
    return frame


def copy_training_pairs(
    frame: pd.DataFrame,
    output_dir: Path,
    artifact_root: Path,
    label_dirs: list[Path],
    min_masks: int,
    corrected_only: bool,
) -> pd.DataFrame:
    rows = []
    for _, row in frame.iterrows():
        src_image = Path(str(row["source_image_abs"]))
        require_exists(src_image)
        corrected = corrected_label_source(src_image, label_dirs)
        if corrected is None:
            if corrected_only:
                continue
            image_arr = tifffile.imread(src_image) if src_image.suffix.lower() in {".tif", ".tiff"} else None
            if image_arr is None:
                from PIL import Image

                image_shape = np.asarray(Image.open(src_image)).shape[:2]
            else:
                image_shape = image_arr.shape[:2]
            mask_source_kind, labels = seed_labels(
                artifact_root,
                str(row["filename"]),
                str(row["tile_name"]),
                tuple(int(x) for x in image_shape),
            )
            annotation_path = binary_mask_path(artifact_root, str(row["filename"]), str(row["tile_name"]))
        else:
            annotation_path, mask_source_kind, labels = corrected
        labels = labeled_u16(labels)
        n_masks = int(labels.max(initial=0))
        if n_masks < min_masks:
            continue

        split_dir = output_dir / str(row["split"])
        split_dir.mkdir(parents=True, exist_ok=True)
        dst_image = split_dir / src_image.name
        dst_mask = split_dir / f"{src_image.stem}_masks.png"
        shutil.copy2(src_image, dst_image)
        dst_mask.write_bytes(png_bytes(labels))
        rows.append(
            {
                "split": str(row["split"]),
                "filename": str(row["filename"]),
                "species": str(row["species"]),
                "tile_name": str(row["tile_name"]),
                "n_masks": n_masks,
                "mask_source_kind": mask_source_kind,
                "annotation_path": str(annotation_path),
                "image_path": str(dst_image.relative_to(output_dir)),
                "mask_path": str(dst_mask.relative_to(output_dir)),
            }
        )
    return pd.DataFrame(rows)


def write_readme(output_dir: Path, summary: dict[str, object]) -> None:
    lines = [
        "# Nucleus Tile Training Round",
        "",
        "This directory mirrors the tile split from `tile_training_round_v1`, but with nucleus-instance masks.",
        "",
        "## Contents",
        "",
        "- `train/`: training tiles and `_masks.png` labels",
        "- `test/`: held-out tiles and `_masks.png` labels",
        "- `training_manifest.csv`: per-tile nucleus training metadata",
        "",
        "## Suggested Training Command",
        "",
        "```bash",
        f"uv run python scripts/run_nucleus_tile_training_cycle.py --use-gpu",
        "```",
        "",
        "Drop `--use-gpu` if you want CPU training.",
        "",
        "## Summary",
        "",
    ]
    for key, value in sorted(summary.items()):
        lines.append(f"- `{key}`: `{value}`")
    (output_dir / "README.md").write_text("\n".join(lines) + "\n")


def main() -> None:
    args = parse_args()
    require_exists(args.source_round)
    require_exists(args.artifact_root)
    label_dirs = [path.resolve() for path in args.label_dir]
    for path in label_dirs:
        require_exists(path)

    frame = load_rows(args.source_round)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    clear_dir_except(args.output_dir / "train", keep_names={"models"})
    clear_dir_except(args.output_dir / "test")
    training_manifest = copy_training_pairs(
        frame,
        args.output_dir,
        args.artifact_root,
        label_dirs,
        args.min_masks,
        args.corrected_only,
    )
    if training_manifest.empty:
        if args.corrected_only:
            raise SystemExit("No manually corrected nucleus-labeled tiles found with usable masks")
        raise SystemExit("No nucleus-labeled tiles found with usable masks")
    training_manifest.to_csv(args.output_dir / "training_manifest.csv", index=False)

    summary = {
        "source_round": str(args.source_round.resolve()),
        "artifact_root": str(args.artifact_root.resolve()),
        "label_dirs": [str(path) for path in label_dirs],
        "n_tiles": int(len(training_manifest)),
        "n_train": int((training_manifest["split"] == "train").sum()),
        "n_test": int((training_manifest["split"] == "test").sum()),
        "n_species": int(training_manifest["species"].nunique()),
        "n_seed_roi_zip": int((training_manifest["mask_source_kind"] == "roi_zip").sum()),
        "n_seed_binary_components": int((training_manifest["mask_source_kind"] == "binary_components").sum()),
        "n_corrected": int(training_manifest["mask_source_kind"].isin({"seg_npy", "cp_masks_png", "cp_masks_tif"}).sum()),
        "corrected_only": bool(args.corrected_only),
        "min_masks": int(args.min_masks),
        "median_masks_per_tile": int(training_manifest["n_masks"].median()),
        "max_masks_per_tile": int(training_manifest["n_masks"].max()),
    }
    (args.output_dir / "summary.json").write_text(json.dumps(summary, indent=2))
    write_readme(args.output_dir, summary)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
