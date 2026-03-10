#!/usr/bin/env python3
"""Convert Cellpose tile predictions into linkage-ready cell measurements."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import tifffile
from PIL import Image

PROJECT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT / "scripts"))

from build_yolo_nucleus_dataset import load_rows, row_image_name, row_source_image_path  # noqa: E402
from nucleus_iod_python import load_background_cache  # noqa: E402
from run_yolo_tile_measurements import (  # noqa: E402
    CELL_COLUMNS,
    grayscale_uint8,
    measure_objects,
    write_csv,
)
from stage_yolo_nucleus_predictions import filter_mask_instances  # noqa: E402


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--prediction-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--background-cache", type=Path, default=None)
    parser.add_argument("--min-mask-area", type=int, default=16)
    parser.add_argument("--max-mask-area", type=int, default=0)
    parser.add_argument("--min-circularity", type=float, default=0.0)
    parser.add_argument("--min-solidity", type=float, default=0.0)
    parser.add_argument("--max-aspect-ratio", type=float, default=0.0)
    return parser.parse_args()


def require_exists(path: Path) -> None:
    if not path.exists():
        raise FileNotFoundError(path)


def clear_dir(path: Path) -> None:
    if path.exists():
        for child in path.iterdir():
            if child.is_dir():
                clear_dir(child)
                child.rmdir()
            else:
                child.unlink()
    else:
        path.mkdir(parents=True, exist_ok=True)


def cellpose_mask_path(prediction_dir: Path, image_name: str) -> Path:
    return prediction_dir / f"{Path(image_name).stem}_masks.png"


def main() -> None:
    args = parse_args()
    require_exists(args.manifest)
    require_exists(args.prediction_dir)

    rows = load_rows(args.manifest)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    clear_dir(args.output_dir)
    masks_root = args.output_dir / "masks"
    masks_root.mkdir(parents=True, exist_ok=True)
    bg_cache = load_background_cache(args.background_cache)

    all_measurements: list[dict[str, object]] = []
    tile_manifest_rows: list[dict[str, object]] = []
    total_raw_instances = 0
    total_filtered_instances = 0
    total_positive_pixels = 0
    missing_masks: list[str] = []

    for row in rows:
        image_name = row_image_name(args.manifest, row)
        raw_path = row_source_image_path(args.manifest, row)
        require_exists(raw_path)

        prediction_path = cellpose_mask_path(args.prediction_dir, image_name)
        if not prediction_path.exists():
            missing_masks.append(image_name)
            continue

        image = np.asarray(Image.open(raw_path))
        tile_gray = grayscale_uint8(image)
        raw_mask = np.asarray(Image.open(prediction_path))
        raw_mask = np.squeeze(raw_mask).astype(np.uint32, copy=False)
        filtered_mask = filter_mask_instances(
            raw_mask,
            min_area=args.min_mask_area,
            max_area=args.max_mask_area,
            min_circularity=args.min_circularity,
            min_solidity=args.min_solidity,
            max_aspect_ratio=args.max_aspect_ratio,
        ).astype(np.uint32, copy=False)

        total_raw_instances += int(raw_mask.max(initial=0))
        total_filtered_instances += int(filtered_mask.max(initial=0))
        total_positive_pixels += int((filtered_mask > 0).sum())

        filename = str(row.get("filename", "")).strip()
        tile_name = str(row.get("tile_name", "")).strip() or f"{Path(image_name).stem}.tiff"
        image_stem = Path(filename).stem if filename else Path(image_name).stem
        mask_dir = masks_root / image_stem
        mask_dir.mkdir(parents=True, exist_ok=True)
        raw_mask_path = mask_dir / f"{tile_name[:-5]}__raw_labels.tiff"
        mask_path = mask_dir / f"{tile_name[:-5]}__filtered_labels.tiff"
        tifffile.imwrite(str(raw_mask_path), raw_mask, compression="zlib")
        tifffile.imwrite(str(mask_path), filtered_mask, compression="zlib")

        tile_manifest_rows.append(
            {
                "tile_name": tile_name,
                "tile_y0": int(str(row.get("tile_y0", "")).strip() or 0),
                "tile_x0": int(str(row.get("tile_x0", "")).strip() or 0),
                "tile_h": int(filtered_mask.shape[0]),
                "tile_w": int(filtered_mask.shape[1]),
                "tile_score": float(str(row.get("tile_score", "")).strip() or 1.0),
                "tile_status": "processed",
                "mask_path": str(mask_path.resolve()),
                "raw_mask_path": str(raw_mask_path.resolve()),
            }
        )

        i_bg = float(max(1.0, bg_cache.get(filename, np.percentile(tile_gray, 95))))
        all_measurements.extend(
            measure_objects(
                mask=filtered_mask,
                tile_gray=tile_gray,
                i_bg=i_bg,
                row=row,
                source_image_path=raw_path,
                run_manifest_path=args.output_dir / "summary.json",
                object_kind="cell",
                mask_path=mask_path,
                tile_manifest_path=args.output_dir / "tile_manifest.csv",
            )
        )

    if tile_manifest_rows:
        write_csv(args.output_dir / "tile_manifest.csv", tile_manifest_rows, list(tile_manifest_rows[0].keys()))
    else:
        write_csv(
            args.output_dir / "tile_manifest.csv",
            [],
            ["tile_name", "tile_y0", "tile_x0", "tile_h", "tile_w", "tile_score", "tile_status", "mask_path", "raw_mask_path"],
        )
    write_csv(args.output_dir / "all_measurements.csv", all_measurements, CELL_COLUMNS)

    summary = {
        "manifest": str(args.manifest.resolve()),
        "prediction_dir": str(args.prediction_dir.resolve()),
        "output_dir": str(args.output_dir.resolve()),
        "n_tiles": len(rows),
        "n_tiles_with_predictions": len(rows) - len(missing_masks),
        "n_missing_masks": len(missing_masks),
        "missing_masks": missing_masks,
        "n_objects": len(all_measurements),
        "n_total_predicted_instances_raw": total_raw_instances,
        "n_total_predicted_instances_filtered": total_filtered_instances,
        "n_total_positive_pixels": total_positive_pixels,
        "min_mask_area": int(args.min_mask_area),
        "max_mask_area": int(args.max_mask_area),
        "min_circularity": float(args.min_circularity),
        "min_solidity": float(args.min_solidity),
        "max_aspect_ratio": float(args.max_aspect_ratio),
    }
    (args.output_dir / "summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
