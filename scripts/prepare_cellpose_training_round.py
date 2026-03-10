#!/usr/bin/env python3
"""Build Cellpose train/test folders from reviewed annotation-pack crops."""

from __future__ import annotations

import argparse
import json
import shutil
from pathlib import Path

import numpy as np
import pandas as pd
import tifffile

PROJECT = Path(__file__).resolve().parent.parent
DEFAULT_PACK = PROJECT / "output" / "annotation_pack_v1"
DEFAULT_OUTPUT = PROJECT / "output" / "cellpose_training_round_v1"
APPROVED_STATUSES = {"approved", "ready", "corrected"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--annotation-pack", type=Path, default=DEFAULT_PACK)
    parser.add_argument("--review-status", type=Path, default=None)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--test-fraction", type=float, default=0.2)
    parser.add_argument("--seed", type=int, default=42)
    return parser.parse_args()


def require_exists(path: Path) -> None:
    if not path.exists():
        raise FileNotFoundError(path)


def normalize_status(series: pd.Series) -> pd.Series:
    return series.fillna("").astype(str).str.strip().str.lower()


def png_bytes(arr: np.ndarray) -> bytes:
    from io import BytesIO
    from PIL import Image

    image = Image.fromarray(arr)
    handle = BytesIO()
    image.save(handle, format="PNG", compress_level=1)
    return handle.getvalue()


def labeled_u16(mask: np.ndarray) -> np.ndarray:
    arr = np.asarray(mask)
    if arr.ndim > 2:
        arr = np.squeeze(arr)
    if arr.dtype.kind in {"b", "u", "i"} and arr.max(initial=0) <= np.iinfo(np.uint16).max:
        return arr.astype(np.uint16, copy=False)
    arr = np.rint(arr).astype(np.int64, copy=False)
    arr[arr < 0] = 0
    if arr.max(initial=0) > np.iinfo(np.uint16).max:
        raise ValueError("mask labels exceed uint16 range")
    return arr.astype(np.uint16, copy=False)


def resolve_corrected_mask(image_path: Path, seed_mask_path: Path) -> tuple[str, np.ndarray | None]:
    base = image_path.with_suffix("")
    seg_npy = base.parent / f"{base.name}_seg.npy"
    cp_png = base.parent / f"{base.name}_cp_masks.png"
    cp_tif = base.parent / f"{base.name}_cp_masks.tif"
    if seg_npy.exists():
        dat = np.load(seg_npy, allow_pickle=True).item()
        masks = dat.get("masks")
        if masks is None:
            raise ValueError(f"{seg_npy} does not contain 'masks'")
        return "seg_npy", labeled_u16(masks)
    if cp_png.exists():
        return "cp_masks_png", labeled_u16(tifffile.imread(cp_png))
    if cp_tif.exists():
        return "cp_masks_tif", labeled_u16(tifffile.imread(cp_tif))
    return "seed_mask", None


def split_reviewed(df: pd.DataFrame, status_col: str, test_fraction: float, seed: int) -> pd.DataFrame:
    reviewed = df[normalize_status(df[status_col]).isin(APPROVED_STATUSES)].copy()
    if reviewed.empty:
        reviewed["split"] = []
        return reviewed

    rng = np.random.default_rng(seed)
    reviewed = reviewed.sample(frac=1.0, random_state=seed).reset_index(drop=True)
    splits = []
    for _, grp in reviewed.groupby("species", sort=False):
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
            splits.append((row_idx, "test" if row_idx in test_idx else "train"))
    split_map = {row_idx: split for row_idx, split in splits}
    reviewed["split"] = reviewed.index.map(split_map)
    return reviewed


def copy_pairs(reviewed: pd.DataFrame, family: str, output_dir: Path) -> pd.DataFrame:
    rows = []
    image_col = f"{family}_image_path"
    mask_col = f"{family}_mask_path"
    for _, row in reviewed.iterrows():
        split_dir = output_dir / family / row["split"]
        split_dir.mkdir(parents=True, exist_ok=True)
        src_image = Path(row[image_col])
        dst_image = split_dir / src_image.name
        dst_mask = split_dir / f"{src_image.stem}_masks.png"
        shutil.copy2(src_image, dst_image)
        mask_source_kind, corrected = resolve_corrected_mask(src_image, Path(row[mask_col]))
        if corrected is None:
            shutil.copy2(Path(row[mask_col]), dst_mask)
        else:
            dst_mask.write_bytes(png_bytes(corrected))
        rows.append(
            {
                "family": family,
                "split": row["split"],
                "sample_id": row["sample_id"],
                "species": row["species"],
                "annotation_bucket": row["annotation_bucket"],
                "image_path": str(dst_image.relative_to(output_dir)),
                "mask_path": str(dst_mask.relative_to(output_dir)),
                "mask_source_kind": mask_source_kind,
            }
        )
    return pd.DataFrame(rows)


def write_readme(output_dir: Path, summary: dict) -> None:
    lines = [
        "# Cellpose Training Round",
        "",
        "This directory contains only reviewed crops that were explicitly marked ready for training.",
        "",
        "## Review Statuses Used",
        "",
        "- Accepted statuses: `approved`, `ready`, `corrected`",
        "- Any blank or other status is excluded from the training split",
        "",
        "## Suggested Training Commands",
        "",
        "```bash",
        f"python3 -m cellpose --train --dir {output_dir / 'cellpose_cell' / 'train'} --test_dir {output_dir / 'cellpose_cell' / 'test'} --mask_filter _masks --pretrained_model cyto3 --model_name_out desmognathus_cell_round1",
        f"python3 -m cellpose --train --dir {output_dir / 'cellpose_nucleus' / 'train'} --test_dir {output_dir / 'cellpose_nucleus' / 'test'} --mask_filter _masks --pretrained_model nuclei --model_name_out desmognathus_nucleus_round1",
        "```",
        "",
        "Add `--use_gpu` if you want CUDA training.",
        "",
        "## Summary",
        "",
    ]
    for key, value in sorted(summary.items()):
        lines.append(f"- `{key}`: `{value}`")
    (output_dir / "README.md").write_text("\n".join(lines) + "\n")


def main() -> None:
    args = parse_args()
    require_exists(args.annotation_pack)
    manifest_path = args.annotation_pack / "manifest.csv"
    require_exists(manifest_path)
    review_status = args.review_status or (args.annotation_pack / "review_status_template.csv")
    require_exists(review_status)

    manifest = pd.read_csv(manifest_path)
    review = pd.read_csv(review_status)
    merged = manifest.merge(
        review[
            [
                "sample_id",
                "cell_review_status",
                "nucleus_review_status",
                "reviewer",
                "review_notes",
            ]
        ],
        on="sample_id",
        how="left",
        validate="one_to_one",
    )

    for col in [
        "image_path",
        "cell_seed_path",
        "nucleus_seed_path",
        "overlay_path",
        "cellpose_cell_image_path",
        "cellpose_cell_mask_path",
        "cellpose_nucleus_image_path",
        "cellpose_nucleus_mask_path",
    ]:
        merged[col] = merged[col].map(lambda p: str((args.annotation_pack / p).resolve()))

    cell = split_reviewed(merged, "cell_review_status", args.test_fraction, args.seed)
    nucleus = split_reviewed(merged, "nucleus_review_status", args.test_fraction, args.seed + 1)

    args.output_dir.mkdir(parents=True, exist_ok=True)
    cell_manifest = copy_pairs(cell, "cellpose_cell", args.output_dir) if len(cell) else pd.DataFrame()
    nucleus_manifest = copy_pairs(nucleus, "cellpose_nucleus", args.output_dir) if len(nucleus) else pd.DataFrame()
    training_manifest = pd.concat([cell_manifest, nucleus_manifest], ignore_index=True)
    training_manifest.to_csv(args.output_dir / "training_manifest.csv", index=False)

    summary = {
        "annotation_pack": str(args.annotation_pack),
        "review_status": str(review_status),
        "n_cell_reviewed_ready": int(len(cell)),
        "n_cell_train": int((cell["split"] == "train").sum()) if len(cell) else 0,
        "n_cell_test": int((cell["split"] == "test").sum()) if len(cell) else 0,
        "n_nucleus_reviewed_ready": int(len(nucleus)),
        "n_nucleus_train": int((nucleus["split"] == "train").sum()) if len(nucleus) else 0,
        "n_nucleus_test": int((nucleus["split"] == "test").sum()) if len(nucleus) else 0,
        "n_species_cell": int(cell["species"].nunique()) if len(cell) else 0,
        "n_species_nucleus": int(nucleus["species"].nunique()) if len(nucleus) else 0,
    }
    (args.output_dir / "summary.json").write_text(json.dumps(summary, indent=2))
    write_readme(args.output_dir, summary)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
