#!/usr/bin/env python3
"""Run a YOLO segmentation model on tile images and emit linkage-ready masks and measurements."""

from __future__ import annotations

import argparse
import csv
import json
import re
import sys
from pathlib import Path

import numpy as np
import tifffile
from PIL import Image
from skimage.measure import regionprops

PROJECT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT / "scripts"))

from build_yolo_nucleus_dataset import load_rows, row_image_name, row_source_image_path  # noqa: E402
from nucleus_iod_python import PIXEL_AREA_UM2, load_background_cache  # noqa: E402
from stage_yolo_nucleus_predictions import choose_device, tile_predictions  # noqa: E402

CELL_COLUMNS = [
    "label",
    "area_px",
    "area_um2",
    "solidity",
    "circularity",
    "iod",
    "mean_od",
    "centroid_y",
    "centroid_x",
    "i_bg",
    "i_bg_source",
    "tile_name",
    "tile_y0",
    "tile_x0",
    "tile_h",
    "tile_w",
    "tile_score",
    "mask_path",
    "raw_mask_path",
    "mask_label_id",
    "tile_manifest_path",
    "source_image_path",
    "run_manifest_path",
    "overlay_path",
    "filename",
    "slide_id",
    "specimen_id",
    "species",
    "image_type",
]

NUCLEUS_COLUMNS = [
    "filename",
    "label",
    "area_px",
    "area_um2",
    "iod",
    "mean_od",
    "centroid_x",
    "centroid_y",
    "i_bg",
    "i_bg_source",
    "tile_name",
    "tile_y0",
    "tile_x0",
    "tile_height_px",
    "tile_width_px",
    "mask_path",
    "roi_zip_path",
    "tile_manifest_path",
    "source_image_path",
    "run_manifest_path",
    "raw_imagej_results_path",
    "slide_id",
    "specimen_id",
    "species",
    "image_type",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--model", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--object-kind", choices=["cell", "nucleus"], required=True)
    parser.add_argument("--image-type", default="brightfield", choices=["brightfield"])
    parser.add_argument("--background-cache", type=Path, default=None)
    parser.add_argument(
        "--allow-missing-background-cache",
        action="store_true",
        help="Permit per-tile 95th-percentile background fallback when the source image is absent from the cache.",
    )
    parser.add_argument("--patch-size", type=int, default=1024)
    parser.add_argument("--stride", type=int, default=1024)
    parser.add_argument("--imgsz", type=int, default=1024)
    parser.add_argument("--batch", type=int, default=8)
    parser.add_argument("--conf", type=float, default=0.10)
    parser.add_argument("--iou", type=float, default=0.50)
    parser.add_argument("--max-det", type=int, default=512)
    parser.add_argument("--min-mask-area", type=int, default=16)
    parser.add_argument("--max-mask-area", type=int, default=0)
    parser.add_argument("--min-circularity", type=float, default=0.0)
    parser.add_argument("--min-solidity", type=float, default=0.0)
    parser.add_argument("--max-aspect-ratio", type=float, default=0.0)
    parser.add_argument("--device", default=None)
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


def parse_tile_origin(tile_name: str) -> tuple[int, int]:
    match = re.search(r"tile_y(\d+)_x(\d+)\.tiff$", str(tile_name))
    if not match:
        return 0, 0
    return int(match.group(1)), int(match.group(2))


def slide_id_from_filename(filename: str) -> str:
    match = re.search(r"Process_(\d+)", filename)
    return match.group(1) if match else ""


def grayscale_uint8(image: np.ndarray) -> np.ndarray:
    arr = np.asarray(image)
    if arr.ndim == 3 and arr.shape[-1] in (3, 4):
        arr = arr[..., 0]
    elif arr.ndim == 3 and arr.shape[0] in (1, 3, 4):
        arr = arr[0]
    arr = np.squeeze(arr)
    if arr.dtype == np.uint16:
        return (arr >> 8).astype(np.uint8)
    if arr.dtype != np.uint8:
        return np.clip(arr, 0, 255).astype(np.uint8)
    return arr


def tile_background(
    tile_gray: np.ndarray,
    filename: str,
    cache: dict[str, int],
    *,
    allow_missing_cache: bool,
) -> tuple[float, str]:
    cached = cache.get(filename)
    if cached is not None:
        return float(cached), "image_cache"
    if not allow_missing_cache:
        raise KeyError(
            f"Missing cached I_bg for {filename!r}. "
            "Provide a complete --background-cache or pass --allow-missing-background-cache."
        )
    return float(max(1.0, np.percentile(tile_gray, 95))), "tile_p95_fallback"


def measure_objects(
    *,
    mask: np.ndarray,
    tile_gray: np.ndarray,
    i_bg: float,
    i_bg_source: str,
    row: dict[str, str],
    source_image_path: Path,
    run_manifest_path: Path,
    object_kind: str,
    mask_path: Path,
    tile_manifest_path: Path,
) -> list[dict[str, object]]:
    props = regionprops(mask)
    filename = str(row.get("filename", "")).strip()
    tile_name = str(row.get("tile_name", "")).strip() or f"{Path(mask_path).stem}.tiff"
    tile_y0 = int(str(row.get("tile_y0", "")).strip() or parse_tile_origin(tile_name)[0])
    tile_x0 = int(str(row.get("tile_x0", "")).strip() or parse_tile_origin(tile_name)[1])
    tile_h, tile_w = mask.shape
    species = str(row.get("species", "")).strip()
    specimen_id = str(row.get("specimen_id", "")).strip()
    slide_id = str(row.get("slide_id", "")).strip() or slide_id_from_filename(filename)
    image_type = str(row.get("image_type", "")).strip() or "brightfield"
    tile_score = float(str(row.get("tile_score", "")).strip() or 1.0)

    measurements: list[dict[str, object]] = []
    for prop in props:
        area_px = int(prop.area)
        perimeter = float(prop.perimeter)
        circularity = 0.0 if perimeter <= 0 else float(4 * np.pi * area_px / (perimeter * perimeter))
        coords = prop.coords
        pixel_vals = np.clip(tile_gray[coords[:, 0], coords[:, 1]].astype(np.float64), 1, None)
        od_per_pixel = np.log10(max(i_bg, 1.0) / pixel_vals)
        iod = float(np.sum(od_per_pixel))
        mean_od = float(np.mean(od_per_pixel))
        centroid_y = float(prop.centroid[0] + tile_y0)
        centroid_x = float(prop.centroid[1] + tile_x0)

        if object_kind == "cell":
            measurements.append(
                {
                    "label": int(prop.label),
                    "area_px": area_px,
                    "area_um2": round(area_px * PIXEL_AREA_UM2, 4),
                    "solidity": round(float(prop.solidity), 6),
                    "circularity": round(circularity, 6),
                    "iod": round(iod, 6),
                    "mean_od": round(mean_od, 6),
                    "centroid_y": round(centroid_y, 2),
                    "centroid_x": round(centroid_x, 2),
                    "i_bg": round(i_bg, 4),
                    "i_bg_source": i_bg_source,
                    "tile_name": tile_name,
                    "tile_y0": tile_y0,
                    "tile_x0": tile_x0,
                    "tile_h": tile_h,
                    "tile_w": tile_w,
                    "tile_score": tile_score,
                    "mask_path": str(mask_path.resolve()),
                    "raw_mask_path": str(mask_path.resolve()),
                    "mask_label_id": int(prop.label),
                    "tile_manifest_path": str(tile_manifest_path.resolve()),
                    "source_image_path": str(source_image_path.resolve()),
                    "run_manifest_path": str(run_manifest_path.resolve()),
                    "overlay_path": "",
                    "filename": filename,
                    "slide_id": slide_id,
                    "specimen_id": specimen_id,
                    "species": species,
                    "image_type": image_type,
                }
            )
        else:
            measurements.append(
                {
                    "filename": filename,
                    "label": int(prop.label),
                    "area_px": area_px,
                    "area_um2": round(area_px * PIXEL_AREA_UM2, 4),
                    "iod": round(iod, 6),
                    "mean_od": round(mean_od, 6),
                    "centroid_x": round(centroid_x, 2),
                    "centroid_y": round(centroid_y, 2),
                    "i_bg": round(i_bg, 4),
                    "i_bg_source": i_bg_source,
                    "tile_name": tile_name,
                    "tile_y0": tile_y0,
                    "tile_x0": tile_x0,
                    "tile_height_px": tile_h,
                    "tile_width_px": tile_w,
                    "mask_path": str(mask_path.resolve()),
                    "roi_zip_path": "",
                    "tile_manifest_path": str(tile_manifest_path.resolve()),
                    "source_image_path": str(source_image_path.resolve()),
                    "run_manifest_path": str(run_manifest_path.resolve()),
                    "raw_imagej_results_path": "",
                    "slide_id": slide_id,
                    "specimen_id": specimen_id,
                    "species": species,
                    "image_type": image_type,
                }
            )
    return measurements


def write_csv(path: Path, rows: list[dict[str, object]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    args = parse_args()
    require_exists(args.manifest)
    require_exists(args.model)

    rows = load_rows(args.manifest)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    clear_dir(args.output_dir)
    masks_root = args.output_dir / "masks"
    masks_root.mkdir(parents=True, exist_ok=True)

    from ultralytics import YOLO

    device = choose_device(args.device)
    model = YOLO(str(args.model))
    bg_cache = load_background_cache(args.background_cache)
    require_cached_bg = args.object_kind == "nucleus" and not args.allow_missing_background_cache
    if require_cached_bg and not bg_cache:
        raise ValueError(
            "Nucleus IOD measurement requires a populated --background-cache unless "
            "--allow-missing-background-cache is set."
        )

    all_measurements: list[dict[str, object]] = []
    tile_manifest_rows: list[dict[str, object]] = []
    total_instances = 0
    total_positive_pixels = 0
    missing_bg_filenames: set[str] = set()
    bg_source_counts: dict[str, int] = {}

    for row in rows:
        image_name = row_image_name(args.manifest, row)
        raw_path = row_source_image_path(args.manifest, row)
        require_exists(raw_path)

        image = np.asarray(Image.open(raw_path))
        full_mask, n_instances = tile_predictions(
            model=model,
            image=image,
            patch_size=args.patch_size,
            stride=args.stride,
            imgsz=args.imgsz,
            batch=args.batch,
            conf=args.conf,
            iou=args.iou,
            max_det=args.max_det,
            device=device,
            min_mask_area=args.min_mask_area,
            max_mask_area=args.max_mask_area,
            min_circularity=args.min_circularity,
            min_solidity=args.min_solidity,
            max_aspect_ratio=args.max_aspect_ratio,
        )
        total_instances += n_instances
        total_positive_pixels += int((full_mask > 0).sum())

        filename = str(row.get("filename", "")).strip()
        tile_name = str(row.get("tile_name", "")).strip() or f"{Path(image_name).stem}.tiff"
        image_stem = Path(filename).stem if filename else Path(image_name).stem
        mask_dir = masks_root / image_stem
        mask_dir.mkdir(parents=True, exist_ok=True)
        mask_path = mask_dir / f"{tile_name[:-5]}__labels.tiff"
        tifffile.imwrite(str(mask_path), full_mask.astype(np.uint32), compression="zlib")

        tile_y0, tile_x0 = parse_tile_origin(tile_name)
        tile_manifest_rows.append(
            {
                "tile_name": tile_name,
                "tile_y0": int(str(row.get("tile_y0", "")).strip() or tile_y0),
                "tile_x0": int(str(row.get("tile_x0", "")).strip() or tile_x0),
                "tile_h": int(full_mask.shape[0]),
                "tile_w": int(full_mask.shape[1]),
                "tile_score": float(str(row.get("tile_score", "")).strip() or 1.0),
                "tile_status": "processed",
                "mask_path": str(mask_path.resolve()),
                "raw_mask_path": str(mask_path.resolve()),
            }
        )

        tile_gray = grayscale_uint8(image)
        try:
            i_bg, i_bg_source = tile_background(
                tile_gray,
                filename,
                bg_cache,
                allow_missing_cache=not require_cached_bg,
            )
        except KeyError:
            missing_bg_filenames.add(filename)
            raise
        bg_source_counts[i_bg_source] = bg_source_counts.get(i_bg_source, 0) + 1
        all_measurements.extend(
            measure_objects(
                mask=full_mask,
                tile_gray=tile_gray,
                i_bg=i_bg,
                i_bg_source=i_bg_source,
                row=row,
                source_image_path=raw_path,
                run_manifest_path=args.output_dir / "summary.json",
                object_kind=args.object_kind,
                mask_path=mask_path,
                tile_manifest_path=args.output_dir / "tile_manifest.csv",
            )
        )

    write_csv(args.output_dir / "tile_manifest.csv", tile_manifest_rows, list(tile_manifest_rows[0].keys()))
    csv_columns = CELL_COLUMNS if args.object_kind == "cell" else NUCLEUS_COLUMNS
    write_csv(args.output_dir / "all_measurements.csv", all_measurements, csv_columns)

    summary = {
        "manifest": str(args.manifest.resolve()),
        "model": str(args.model.resolve()),
        "output_dir": str(args.output_dir.resolve()),
        "object_kind": args.object_kind,
        "image_type": args.image_type,
        "background_cache": str(args.background_cache.resolve()) if args.background_cache else "",
        "allow_missing_background_cache": bool(args.allow_missing_background_cache),
        "bg_source_counts": bg_source_counts,
        "n_missing_bg_filenames": len(missing_bg_filenames),
        "device": device,
        "n_tiles": len(rows),
        "n_objects": len(all_measurements),
        "n_total_predicted_instances": int(total_instances),
        "n_total_positive_pixels": int(total_positive_pixels),
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
