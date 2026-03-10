#!/usr/bin/env python3
"""Fast rule-based brightfield cell segmentation.

This backend mirrors the existing nucleus-first brightfield logic without
Cellpose:
  1. Segment nuclei with the same threshold/blur/watershed workflow used by
     the Python nucleus IOD pipeline.
  2. Build a permissive cell foreground mask on the original brightfield tile.
  3. Grow one cell mask per nucleus by watershed within that foreground.
  4. Filter by realistic geometry and nucleus-in-cell constraints.

Outputs match the existing cell-size contract closely enough for the current
linkage/QC tooling:
  - per-image measurement CSVs
  - per-image or per-tile label masks
  - tile manifests
  - optional crops
  - optional debug overlays
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
import time
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import tifffile
import zarr
from scipy import ndimage
from skimage.filters import gaussian, threshold_otsu
from skimage.measure import regionprops
from skimage.morphology import closing, disk
from skimage.segmentation import find_boundaries, watershed

PROJECT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT / "scripts"))

from nucleus_iod_python import (  # noqa: E402
    EDGE_MARGIN,
    PIXEL_AREA_UM2,
    TILE_SIZE,
    _normalize_uint8,
    classify_tile,
    filter_nuclei,
    segment_nuclei,
)

CSV_COLUMNS = [
    "label",
    "seed_label",
    "area_px",
    "area_um2",
    "solidity",
    "circularity",
    "iod",
    "mean_od",
    "centroid_y",
    "centroid_x",
    "i_bg",
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
    "overlay_path",
    "nucleus_area_px_seed",
    "nucleus_area_um2_seed",
    "nc_ratio_seed",
    "distance_over_cell_radius_seed",
    "cell_extent",
    "cell_edge_touch",
    "cell_threshold",
    "nucleus_threshold",
    "source_image_path",
]

DEFAULT_CROP_SIZE = EDGE_MARGIN * 2


@dataclass
class RuleSettings:
    cell_blur_sigma: float = 1.2
    cell_threshold_correction: float = 0.95
    cell_closing_radius: int = 3
    cell_min_area_px: int = 2500
    cell_max_area_px: int = 26000
    cell_min_circularity: float = 0.35
    cell_min_solidity: float = 0.75
    cell_min_extent: float = 0.40
    nc_ratio_min: float = 0.08
    nc_ratio_max: float = 0.50
    max_distance_over_cell_radius: float = 0.75
    min_cytoplasm_area_px: int = 100
    max_expansion_px: int = 80
    min_foreground_object_px: int = 400


def load_settings(path: Path | None) -> RuleSettings:
    settings = RuleSettings()
    if path is None:
        return settings
    payload = json.loads(Path(path).read_text())
    for key, value in payload.items():
        if hasattr(settings, key):
            setattr(settings, key, value)
    return settings


def _artifact_paths(output_dir: Path, image_stem: str) -> dict[str, Path]:
    root = Path(output_dir)
    return {
        "mask_dir": root / "masks" / image_stem,
        "tile_manifest": root / "tile_manifests" / f"{image_stem}_tile_manifest.csv",
        "debug_dir": root / "debug",
    }


def _remove_small_binary(binary: np.ndarray, min_pixels: int) -> np.ndarray:
    if min_pixels <= 1:
        return np.asarray(binary, dtype=bool)
    labeled, count = ndimage.label(binary)
    if count == 0:
        return np.zeros_like(binary, dtype=bool)
    sizes = np.bincount(labeled.ravel())
    keep = sizes >= int(min_pixels)
    keep[0] = False
    return keep[labeled]


def _write_label_mask(label_image: np.ndarray, mask_path: Path) -> None:
    mask_path = Path(mask_path)
    mask_path.parent.mkdir(parents=True, exist_ok=True)
    tifffile.imwrite(str(mask_path), label_image.astype(np.uint32), compression="zlib")


def _write_tile_manifest(manifest_rows: list[dict], manifest_path: Path) -> None:
    manifest_path = Path(manifest_path)
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    if not manifest_rows:
        return
    fieldnames = [
        "tile_name",
        "tile_y0",
        "tile_x0",
        "tile_h",
        "tile_w",
        "tile_score",
        "tile_status",
        "mask_path",
        "raw_mask_path",
    ]
    with open(manifest_path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(manifest_rows)


def _save_csv(measurements: list[dict], output_dir: Path, name: str) -> None:
    if not measurements:
        return
    csv_path = Path(output_dir) / f"{name}_measurements.csv"
    with open(csv_path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=measurements[0].keys())
        writer.writeheader()
        writer.writerows(measurements)
    print(f"  {name}: saved {len(measurements)} measurements to {csv_path.name}")


def save_cell_crop(
    image: np.ndarray,
    centroid_y: float,
    centroid_x: float,
    label: int,
    crop_dir: Path,
    image_stem: str,
    crop_size: int = DEFAULT_CROP_SIZE,
) -> str:
    half = crop_size // 2
    cy = int(round(float(centroid_y)))
    cx = int(round(float(centroid_x)))
    h, w = image.shape[:2]

    sy0, sy1 = cy - half, cy + half
    sx0, sx1 = cx - half, cx + half

    if image.ndim == 2:
        pad_val = int(np.median(image))
        crop = np.full((crop_size, crop_size), pad_val, dtype=np.uint8)
    else:
        med = np.median(image.reshape(-1, image.shape[-1]), axis=0).astype(np.uint8)
        crop = np.tile(med, (crop_size, crop_size, 1))

    dy0 = max(0, -sy0)
    dx0 = max(0, -sx0)
    dy1 = crop_size - max(0, sy1 - h)
    dx1 = crop_size - max(0, sx1 - w)
    iy0, iy1 = max(0, sy0), min(h, sy1)
    ix0, ix1 = max(0, sx0), min(w, sx1)
    crop[dy0:dy1, dx0:dx1] = image[iy0:iy1, ix0:ix1]

    crop_dir.mkdir(parents=True, exist_ok=True)
    fname = f"{image_stem}_cell_{label:04d}.tiff"
    tifffile.imwrite(str(crop_dir / fname), crop)
    return fname


def _build_nucleus_markers(
    tile: np.ndarray,
    manual_threshold: int | None,
) -> tuple[np.ndarray, dict[int, object], int | None]:
    labeled = segment_nuclei(tile, manual_threshold=manual_threshold)
    props = filter_nuclei(labeled)
    interior = [
        p
        for p in props
        if EDGE_MARGIN < p.centroid[0] < tile.shape[0] - EDGE_MARGIN
        and EDGE_MARGIN < p.centroid[1] < tile.shape[1] - EDGE_MARGIN
    ]
    markers = np.zeros(tile.shape[:2], dtype=np.int32)
    nucleus_props: dict[int, object] = {}
    next_id = 1
    for p in interior:
        markers[labeled == p.label] = next_id
        nucleus_props[next_id] = p
        next_id += 1
    return markers, nucleus_props, manual_threshold


def _cell_foreground(tile: np.ndarray, markers: np.ndarray, settings: RuleSettings) -> tuple[np.ndarray, float]:
    inverted = 255.0 - tile.astype(np.float32)
    smoothed = gaussian(inverted, sigma=settings.cell_blur_sigma, preserve_range=True)
    try:
        otsu = float(threshold_otsu(smoothed))
    except ValueError:
        otsu = 255.0
    cell_threshold = float(np.clip(otsu * settings.cell_threshold_correction, 1.0, 254.0))
    foreground = smoothed > cell_threshold
    if settings.cell_closing_radius > 0:
        foreground = closing(foreground, footprint=disk(int(settings.cell_closing_radius)))
    foreground = ndimage.binary_fill_holes(foreground)
    foreground = _remove_small_binary(foreground, max(1, int(settings.min_foreground_object_px)))

    seed_dist = ndimage.distance_transform_edt(markers == 0)
    foreground &= seed_dist <= int(settings.max_expansion_px)
    foreground |= markers > 0
    return foreground, cell_threshold


def _grow_cells(tile: np.ndarray, markers: np.ndarray, foreground: np.ndarray, settings: RuleSettings) -> np.ndarray:
    surface = ndimage.gaussian_gradient_magnitude(tile.astype(np.float32), sigma=settings.cell_blur_sigma)
    return watershed(surface, markers=markers, mask=foreground)


def _filter_cells(
    raw_cells: np.ndarray,
    nucleus_props: dict[int, object],
    settings: RuleSettings,
) -> tuple[np.ndarray, list[dict], list[dict]]:
    filtered = np.zeros_like(raw_cells, dtype=np.int32)
    kept_rows: list[dict] = []
    rejected_rows: list[dict] = []
    next_id = 1

    for prop in regionprops(raw_cells):
        label_id = int(prop.label)
        nuc = nucleus_props.get(label_id)
        if nuc is None:
            continue
        area_px = int(prop.area)
        cytoplasm_area_px = int(area_px - nuc.area)
        circ = float(4.0 * math.pi * prop.area / (prop.perimeter ** 2)) if prop.perimeter > 0 else 0.0
        solidity = float(prop.solidity)
        extent = float(prop.extent)
        edge_touch = bool(
            prop.bbox[0] == 0
            or prop.bbox[1] == 0
            or prop.bbox[2] >= raw_cells.shape[0]
            or prop.bbox[3] >= raw_cells.shape[1]
        )
        nc_ratio = float(nuc.area / max(prop.area, 1))
        cell_radius = math.sqrt(float(prop.area) / math.pi) if prop.area > 0 else 0.0
        dist_over_r = (
            float(np.hypot(prop.centroid[0] - nuc.centroid[0], prop.centroid[1] - nuc.centroid[1]) / cell_radius)
            if cell_radius > 0
            else float("inf")
        )

        reason = ""
        if edge_touch:
            reason = "edge_touch"
        elif area_px < settings.cell_min_area_px or area_px > settings.cell_max_area_px:
            reason = "cell_area_out_of_range"
        elif circ < settings.cell_min_circularity:
            reason = "cell_low_circularity"
        elif solidity < settings.cell_min_solidity:
            reason = "cell_low_solidity"
        elif extent < settings.cell_min_extent:
            reason = "cell_low_extent"
        elif cytoplasm_area_px < settings.min_cytoplasm_area_px:
            reason = "cytoplasm_too_small"
        elif nc_ratio < settings.nc_ratio_min or nc_ratio > settings.nc_ratio_max:
            reason = "nc_ratio_out_of_range"
        elif dist_over_r > settings.max_distance_over_cell_radius:
            reason = "nucleus_off_center"

        row = {
            "label": label_id,
            "area_px": area_px,
            "area_um2": area_px * PIXEL_AREA_UM2,
            "solidity": solidity,
            "circularity": circ,
            "centroid_y": float(prop.centroid[0]),
            "centroid_x": float(prop.centroid[1]),
            "nucleus_area_px_seed": int(nuc.area),
            "nucleus_area_um2_seed": float(nuc.area * PIXEL_AREA_UM2),
            "nc_ratio_seed": nc_ratio,
            "distance_over_cell_radius_seed": dist_over_r,
            "cell_extent": extent,
            "cell_edge_touch": edge_touch,
            "reject_reason": reason,
        }

        if reason:
            rejected_rows.append(row)
            continue

        filtered[raw_cells == label_id] = next_id
        row["seed_label"] = label_id
        row["label"] = next_id
        kept_rows.append(row)
        next_id += 1

    return filtered, kept_rows, rejected_rows


def _measure_cells(tile: np.ndarray, kept_rows: list[dict], filtered_cells: np.ndarray) -> list[dict]:
    if not kept_rows:
        return []
    bg_mask = filtered_cells == 0
    i_bg = float(np.percentile(tile[bg_mask], 95)) if bg_mask.any() else 255.0
    rows_by_label = {int(row["label"]): dict(row) for row in kept_rows}
    for prop in regionprops(filtered_cells):
        row = rows_by_label.get(int(prop.label))
        if row is None:
            continue
        coords = prop.coords
        pixel_vals = np.clip(tile[coords[:, 0], coords[:, 1]].astype(np.float64), 1, None)
        od_per_pixel = np.log10(i_bg / pixel_vals)
        row["iod"] = float(np.sum(od_per_pixel))
        row["mean_od"] = float(np.mean(od_per_pixel))
        row["i_bg"] = i_bg
    return [rows_by_label[k] for k in sorted(rows_by_label)]


def process_tile(
    tile: np.ndarray,
    tile_y0: int,
    tile_x0: int,
    settings: RuleSettings,
    manual_threshold: int | None = None,
) -> tuple[list[dict], np.ndarray, np.ndarray]:
    markers, nucleus_props, _ = _build_nucleus_markers(tile, manual_threshold=manual_threshold)
    if not nucleus_props:
        return [], np.zeros(tile.shape[:2], dtype=np.int32), np.zeros(tile.shape[:2], dtype=np.int32)

    foreground, cell_threshold = _cell_foreground(tile, markers, settings)
    raw_cells = _grow_cells(tile, markers, foreground, settings)
    filtered_cells, kept_rows, _ = _filter_cells(raw_cells, nucleus_props, settings)
    measured = _measure_cells(tile, kept_rows, filtered_cells)
    for row in measured:
        row["centroid_y"] += tile_y0
        row["centroid_x"] += tile_x0
        row["tile_name"] = f"tile_y{tile_y0:06d}_x{tile_x0:06d}.tiff"
        row["tile_y0"] = tile_y0
        row["tile_x0"] = tile_x0
        row["tile_h"] = int(tile.shape[0])
        row["tile_w"] = int(tile.shape[1])
        row["cell_threshold"] = float(cell_threshold)
        row["nucleus_threshold"] = "" if manual_threshold is None else int(manual_threshold)
    return measured, raw_cells, filtered_cells


def _save_overlay(
    tile: np.ndarray,
    raw_cells: np.ndarray,
    filtered_cells: np.ndarray,
    debug_dir: Path,
    name: str,
) -> Path:
    from PIL import Image

    debug_dir = Path(debug_dir)
    debug_dir.mkdir(parents=True, exist_ok=True)
    rgb = np.stack([tile, tile, tile], axis=-1)
    raw_edge = find_boundaries(raw_cells, mode="outer")
    rgb[raw_edge] = [255, 200, 40]
    kept_edge = find_boundaries(filtered_cells, mode="outer")
    rgb[kept_edge] = [40, 255, 80]
    overlay_path = debug_dir / f"{name}_overlay.png"
    Image.fromarray(rgb).save(str(overlay_path))
    return overlay_path


def process_image(
    image_path: str | Path,
    output_dir: str | Path | None = None,
    debug: bool = False,
    crops: bool = True,
    crop_size: int = DEFAULT_CROP_SIZE,
    tile_filter: str = "auto",
    image_type: str = "brightfield",
    settings: RuleSettings | None = None,
    manual_threshold: int | None = None,
    max_cells: int = 0,
    max_tiles: int = 0,
) -> list[dict]:
    if image_type != "brightfield":
        raise ValueError("rule-based backend currently supports brightfield only")

    settings = settings or RuleSettings()
    image_path = Path(image_path)
    name = image_path.stem
    output_dir = Path(output_dir or PROJECT / "output" / "segmentation")
    output_dir.mkdir(parents=True, exist_ok=True)

    with tifffile.TiffFile(str(image_path)) as tif:
        if not tif.pages:
            raise ValueError(f"TIFF has no pages: {image_path}")
        page = tif.pages[0]
        img_h, img_w = page.shape[:2]

    needs_tiling = img_h > TILE_SIZE * 1.5 or img_w > TILE_SIZE * 1.5
    if needs_tiling:
        return _process_tiled(
            image_path,
            name,
            img_h,
            img_w,
            output_dir,
            debug=debug,
            crops=crops,
            crop_size=crop_size,
            tile_filter=tile_filter,
            settings=settings,
            manual_threshold=manual_threshold,
            max_cells=max_cells,
            max_tiles=max_tiles,
        )
    return _process_small(
        image_path,
        name,
        output_dir,
        debug=debug,
        crops=crops,
        crop_size=crop_size,
        settings=settings,
        manual_threshold=manual_threshold,
    )


def _process_small(
    image_path: Path,
    name: str,
    output_dir: Path,
    *,
    debug: bool,
    crops: bool,
    crop_size: int,
    settings: RuleSettings,
    manual_threshold: int | None,
) -> list[dict]:
    tile = _normalize_uint8(tifffile.imread(str(image_path)))
    print(f"  {name}: {tile.shape[1]}x{tile.shape[0]} (rule-based)")
    measured, raw_cells, filtered_cells = process_tile(tile, 0, 0, settings, manual_threshold=manual_threshold)
    raw_count = int(len(np.unique(raw_cells)) - 1)
    print(f"  {name}: {raw_count} raw -> {len(measured)} kept")

    artifacts = _artifact_paths(output_dir, name)
    raw_mask_path = artifacts["mask_dir"] / f"{name}_raw_labels.tiff"
    mask_path = artifacts["mask_dir"] / f"{name}_filtered_labels.tiff"
    _write_label_mask(raw_cells, raw_mask_path)
    _write_label_mask(filtered_cells, mask_path)

    tile_manifest_rows = [
        {
            "tile_name": f"{name}_full_image.tiff",
            "tile_y0": 0,
            "tile_x0": 0,
            "tile_h": int(tile.shape[0]),
            "tile_w": int(tile.shape[1]),
            "tile_score": 1.0,
            "tile_status": "processed",
            "mask_path": str(mask_path.resolve()),
            "raw_mask_path": str(raw_mask_path.resolve()),
        }
    ]
    _write_tile_manifest(tile_manifest_rows, artifacts["tile_manifest"])

    overlay_path = ""
    if debug:
        overlay_path = str(_save_overlay(tile, raw_cells, filtered_cells, artifacts["debug_dir"], name).resolve())

    if crops and measured:
        crop_dir = output_dir / "crops" / name
        for row in measured:
            fname = save_cell_crop(tile, row["centroid_y"], row["centroid_x"], int(row["label"]), crop_dir, name, crop_size)
            row["crop_filename"] = str((crop_dir / fname).resolve())
            row["crop_size_px"] = crop_size
            row["centroid_in_crop_y"] = crop_size // 2
            row["centroid_in_crop_x"] = crop_size // 2

    for row in measured:
        row["tile_name"] = f"{name}_full_image.tiff"
        row["tile_y0"] = 0
        row["tile_x0"] = 0
        row["tile_h"] = int(tile.shape[0])
        row["tile_w"] = int(tile.shape[1])
        row["tile_score"] = 1.0
        row["mask_path"] = str(mask_path.resolve())
        row["raw_mask_path"] = str(raw_mask_path.resolve())
        row["mask_label_id"] = int(row["label"])
        row["tile_manifest_path"] = str(artifacts["tile_manifest"].resolve())
        row["overlay_path"] = overlay_path
        row["source_image_path"] = str(image_path.resolve())

    _save_csv(measured, output_dir, name)
    return measured


def _process_tiled(
    image_path: Path,
    name: str,
    img_h: int,
    img_w: int,
    output_dir: Path,
    *,
    debug: bool,
    crops: bool,
    crop_size: int,
    tile_filter: str,
    settings: RuleSettings,
    manual_threshold: int | None,
    max_cells: int,
    max_tiles: int,
) -> list[dict]:
    t0 = time.time()
    store = tifffile.imread(str(image_path), aszarr=True)
    z = zarr.open(store, mode="r")
    if isinstance(z, zarr.Group):
        z = z["0"]
    if z.ndim == 3:
        z = z[0]

    artifacts = _artifact_paths(output_dir, name)
    crop_dir = output_dir / "crops" / name if crops else None
    tile_manifest_rows: list[dict] = []

    tile_coords = []
    for y0 in range(0, img_h, TILE_SIZE):
        for x0 in range(0, img_w, TILE_SIZE):
            y1 = min(y0 + TILE_SIZE, img_h)
            x1 = min(x0 + TILE_SIZE, img_w)
            if y1 - y0 >= EDGE_MARGIN * 2 and x1 - x0 >= EDGE_MARGIN * 2:
                tile_coords.append((y0, x0, y1, x1))

    scored = []
    skipped = 0
    skipped_yellow = 0
    for y0, x0, y1, x1 in tile_coords:
        cy, cx = (y0 + y1) // 2, (x0 + x1) // 2
        half = min(128, (y1 - y0) // 2, (x1 - x0) // 2)
        probe = _normalize_uint8(np.asarray(z[cy - half:cy + half, cx - half:cx + half]))
        cls, score = classify_tile(probe)
        tile_name = f"tile_y{y0:06d}_x{x0:06d}.tiff"
        if cls == "red":
            skipped += 1
            tile_manifest_rows.append(
                {
                    "tile_name": tile_name,
                    "tile_y0": y0,
                    "tile_x0": x0,
                    "tile_h": y1 - y0,
                    "tile_w": x1 - x0,
                    "tile_score": "",
                    "tile_status": "skipped_empty",
                    "mask_path": "",
                    "raw_mask_path": "",
                }
            )
            continue
        if tile_filter == "green" and cls == "yellow":
            skipped_yellow += 1
            tile_manifest_rows.append(
                {
                    "tile_name": tile_name,
                    "tile_y0": y0,
                    "tile_x0": x0,
                    "tile_h": y1 - y0,
                    "tile_w": x1 - x0,
                    "tile_score": "",
                    "tile_status": "skipped_borderline",
                    "mask_path": "",
                    "raw_mask_path": "",
                }
            )
            continue
        if tile_filter == "yellow" and cls != "yellow":
            tile_manifest_rows.append(
                {
                    "tile_name": tile_name,
                    "tile_y0": y0,
                    "tile_x0": x0,
                    "tile_h": y1 - y0,
                    "tile_w": x1 - x0,
                    "tile_score": "",
                    "tile_status": "skipped_filter",
                    "mask_path": "",
                    "raw_mask_path": "",
                }
            )
            continue
        scored.append((score, y0, x0, y1, x1, tile_name))

    scored.sort(reverse=True)
    if max_tiles > 0:
        scored = scored[:max_tiles]
    skip_detail = f"{skipped} empty"
    if skipped_yellow:
        skip_detail += f", {skipped_yellow} borderline"
    print(f"  {name}: {img_w}x{img_h} ({len(scored)} tiles, {skip_detail}, filter={tile_filter})")

    all_rows: list[dict] = []
    crop_counter = 0
    overlay_path = ""
    for idx, (score, y0, x0, y1, x1, tile_name) in enumerate(scored):
        tile = _normalize_uint8(np.asarray(z[y0:y1, x0:x1]))
        t_tile = time.time()
        measured, raw_cells, filtered_cells = process_tile(tile, y0, x0, settings, manual_threshold=manual_threshold)
        dt = time.time() - t_tile

        raw_mask_path = artifacts["mask_dir"] / f"{tile_name[:-5]}__raw_labels.tiff"
        mask_path = artifacts["mask_dir"] / f"{tile_name[:-5]}__filtered_labels.tiff"
        _write_label_mask(raw_cells, raw_mask_path)
        _write_label_mask(filtered_cells, mask_path)
        tile_manifest_rows.append(
            {
                "tile_name": tile_name,
                "tile_y0": y0,
                "tile_x0": x0,
                "tile_h": y1 - y0,
                "tile_w": x1 - x0,
                "tile_score": score,
                "tile_status": "processed",
                "mask_path": str(mask_path.resolve()),
                "raw_mask_path": str(raw_mask_path.resolve()),
            }
        )

        if debug and idx == 0:
            overlay_path = str(_save_overlay(tile, raw_cells, filtered_cells, artifacts["debug_dir"], name).resolve())

        for row in measured:
            row["tile_score"] = float(score)
            row["mask_path"] = str(mask_path.resolve())
            row["raw_mask_path"] = str(raw_mask_path.resolve())
            row["mask_label_id"] = int(row["label"])
            row["tile_manifest_path"] = str(artifacts["tile_manifest"].resolve())
            row["overlay_path"] = overlay_path
            row["source_image_path"] = str(image_path.resolve())

        if crop_dir and measured:
            for row in measured:
                crop_counter += 1
                local_y = row["centroid_y"] - y0
                local_x = row["centroid_x"] - x0
                fname = save_cell_crop(tile, local_y, local_x, crop_counter, crop_dir, name, crop_size)
                row["crop_filename"] = str((crop_dir / fname).resolve())
                row["crop_size_px"] = crop_size
                row["centroid_in_crop_y"] = crop_size // 2
                row["centroid_in_crop_x"] = crop_size // 2

        all_rows.extend(measured)
        if idx < 3 or (idx + 1) % 20 == 0:
            print(
                f"    tile {idx + 1}/{len(scored)}: ({y0},{x0}) score={score:.2f} "
                f"{len(measured)} cells [{dt:.1f}s] total={len(all_rows)}",
                flush=True,
            )
        if max_cells > 0 and len(all_rows) >= max_cells:
            print(f"    Reached {max_cells} cell cap after {idx + 1}/{len(scored)} tiles")
            break

    store.close()
    elapsed = time.time() - t0
    print(f"  {name}: {len(all_rows)} kept cells [{elapsed:.0f}s]")

    _write_tile_manifest(tile_manifest_rows, artifacts["tile_manifest"])
    _save_csv(all_rows, output_dir, name)
    return all_rows


def run_batch(
    output_dir: Path,
    *,
    debug: bool = False,
    crops: bool = True,
    crop_size: int = DEFAULT_CROP_SIZE,
    settings: RuleSettings | None = None,
    threshold_csv: Path | None = None,
) -> list[dict]:
    settings = settings or RuleSettings()
    threshold_map = load_threshold_map(threshold_csv)

    master_path = PROJECT / "master_image_metadata.csv"
    if not master_path.exists():
        print(f"ERROR: {master_path} not found")
        sys.exit(1)

    with open(master_path) as handle:
        rows = list(csv.DictReader(handle))

    analysis = [
        row
        for row in rows
        if row.get("analysis_set") == "True" and row.get("image_type") == "brightfield"
    ]
    print(f"Batch: {len(analysis)} brightfield images in analysis set")

    all_results = []
    for i, row in enumerate(analysis, start=1):
        image_path = PROJECT / "data" / "brightfield" / row["filename"]
        if not image_path.exists():
            print(f"  SKIP [{i}/{len(analysis)}] {row['filename']}: not found")
            continue
        print(f"\n[{i}/{len(analysis)}] {row['filename']} ({row.get('species', '?')})")
        try:
            measured = process_image(
                image_path,
                Path(output_dir) / "brightfield",
                debug=debug,
                crops=crops,
                crop_size=crop_size,
                image_type="brightfield",
                settings=settings,
                manual_threshold=threshold_map.get(row["filename"]),
            )
            for item in measured:
                item["filename"] = row["filename"]
                item["slide_id"] = row.get("slide_id", "")
                item["specimen_id"] = row.get("specimen_id", "")
                item["species"] = row.get("species", "")
                item["image_type"] = "brightfield"
            all_results.extend(measured)
        except Exception as exc:
            print(f"  ERROR: {exc}")
            import traceback

            traceback.print_exc()

    if all_results:
        combined_path = Path(output_dir) / "all_measurements.csv"
        with open(combined_path, "w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=all_results[0].keys())
            writer.writeheader()
            writer.writerows(all_results)
        print(f"\nTotal: {len(all_results)} cells -> {combined_path}")
    else:
        print("\nNo measurements generated")
    return all_results


def load_threshold_map(path: Path | None) -> dict[str, int]:
    if path is None or not Path(path).exists():
        return {}
    out: dict[str, int] = {}
    with open(path) as handle:
        for row in csv.DictReader(handle):
            fn = str(row.get("filename", "")).strip()
            if not fn or str(row.get("excluded", "")).strip().lower() == "true":
                continue
            try:
                out[fn] = int(row["threshold"])
            except (KeyError, TypeError, ValueError):
                continue
    return out


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("image", nargs="?", help="Single brightfield image to process")
    parser.add_argument("--output-dir", "-o", type=Path, default=None)
    parser.add_argument("--debug", action="store_true")
    parser.add_argument("--batch", action="store_true")
    parser.add_argument("--crops", action=argparse.BooleanOptionalAction, default=True)
    parser.add_argument("--crop-size", type=int, default=DEFAULT_CROP_SIZE)
    parser.add_argument("--tile-filter", choices=["green", "auto", "yellow", "all"], default="auto")
    parser.add_argument("--threshold-csv", type=Path, help="Per-image threshold CSV")
    parser.add_argument("--settings-json", type=Path, help="Rule-settings JSON override")
    parser.add_argument("--threshold", type=int, default=None, help="Manual nucleus threshold for single-image runs")
    parser.add_argument("--max-cells", type=int, default=0, help="Stop after this many kept cells (0 = no limit)")
    parser.add_argument("--max-tiles", type=int, default=0, help="Process only the top scored tiles (0 = no limit)")
    return parser


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()
    settings = load_settings(args.settings_json)
    if args.image:
        image_path = Path(args.image)
        threshold = args.threshold
        if threshold is None:
            threshold = load_threshold_map(args.threshold_csv).get(image_path.name)
        measured = process_image(
            image_path,
            args.output_dir,
            debug=args.debug,
            crops=args.crops,
            crop_size=args.crop_size,
            tile_filter=args.tile_filter,
            settings=settings,
            manual_threshold=threshold,
            max_cells=args.max_cells,
            max_tiles=args.max_tiles,
        )
        print(f"\n{len(measured)} cells measured")
        if measured:
            areas = [float(row["area_um2"]) for row in measured]
            print(f"  Area: {np.mean(areas):.1f} +/- {np.std(areas):.1f} um^2")
        return
    if not args.batch:
        parser.print_help()
        sys.exit(1)
    run_batch(
        args.output_dir or PROJECT / "output" / "segmentation_rule_based",
        debug=args.debug,
        crops=args.crops,
        crop_size=args.crop_size,
        settings=settings,
        threshold_csv=args.threshold_csv,
    )


if __name__ == "__main__":
    main()
