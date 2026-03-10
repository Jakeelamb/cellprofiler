#!/usr/bin/env python3
"""Segment isolated red blood cells and measure area/IOD.

Uses Cellpose (cpsam) for instance segmentation, then filters for
isolated cells and measures area + IOD.

For each image:
  1. Cellpose instance segmentation (grayscale)
  2. Filter: shape (solidity, circularity)
  3. Filter: keep only isolated cells (centroid distance to nearest neighbor)
  4. Measure area (px -> um^2) and IOD

Usage:
    python Cellsize_segmentation_cellpose_pipeline/scripts/segment_cells.py <image.tiff> [--output-dir DIR] [--debug]
    python Cellsize_segmentation_cellpose_pipeline/scripts/segment_cells.py --batch
"""

import argparse
import csv
import sys
import time
from pathlib import Path

import numpy as np
import tifffile
import zarr
from cellpose import models
from scipy.spatial import KDTree
from skimage import measure as sk_measure, segmentation

def _find_repo_root():
    p = Path(__file__).resolve().parent
    while p != p.parent:
        if (p / ".git").exists():
            return p
        p = p.parent
    raise RuntimeError("Not inside a git repo")

PROJECT = _find_repo_root()
PIXEL_SIZE_UM = 0.12
PIXEL_AREA_UM2 = PIXEL_SIZE_UM ** 2  # 0.0144 um^2

# Segmentation params
MAX_NUCLEUS_DIAMETER_PX = 200
MIN_AREA_PX = int(np.pi * (30 / 2) ** 2 * 0.5)  # ~353
PMOUNT_DIAMETER_PX = 13

# Shape filters
MIN_SOLIDITY = 0.7
MIN_CIRCULARITY = 0.4  # 4*pi*area/perimeter^2

# Isolation: reject if nearest neighbor centroid within this many pixels
NEIGHBOR_DISTANCE_PX = 50  # ~6um at 0.12um/px

# Tiling
TILE_SIZE = 4096
EDGE_MARGIN = MAX_NUCLEUS_DIAMETER_PX  # discard objects near tile edges
MAX_CELLS_PER_IMAGE = 0  # 0 = no limit
GPU_TILE_BATCH_IMAGES = 2
GPU_PATCH_BATCH_SIZE = 16
CPU_TILE_BATCH_IMAGES = 1
CPU_PATCH_BATCH_SIZE = 8

# Cellpose model (lazy-loaded singleton)
_model = None
_model_gpu = False
_model_diameter = None
_model_pretrained = None


def get_model(gpu=None, diameter=None, pretrained_model=None):
    """Get or create Cellpose model. Reuses across calls.

    First call with gpu=True/False creates the model. Subsequent calls
    reuse it regardless of gpu arg (singleton).
    """
    global _model, _model_gpu, _model_diameter, _model_pretrained
    requested_pretrained = pretrained_model or "cpsam"
    if _model is None or _model_pretrained != requested_pretrained:
        use_gpu = gpu if gpu is not None else False
        _model = models.CellposeModel(gpu=use_gpu, pretrained_model=requested_pretrained)
        _model_gpu = use_gpu
        _model_pretrained = requested_pretrained
    _model_diameter = diameter
    return _model


def _tile_batch_images():
    return GPU_TILE_BATCH_IMAGES if _model_gpu else CPU_TILE_BATCH_IMAGES


def _cellpose_patch_batch_size():
    return GPU_PATCH_BATCH_SIZE if _model_gpu else CPU_PATCH_BATCH_SIZE


def _normalize_mask_batch(masks):
    if isinstance(masks, list):
        return masks
    if getattr(masks, "ndim", 0) == 2:
        return [masks]
    if getattr(masks, "ndim", 0) == 3:
        return [masks[i] for i in range(masks.shape[0])]
    raise ValueError(f"Unexpected Cellpose mask output shape: {getattr(masks, 'shape', None)}")


def segment_cells(image, cellpose_batch_size=None):
    """Segment cells using Cellpose (auto-diameter).

    Returns (labeled_array, regionprops_list).
    """
    model = get_model()
    masks, _, _ = model.eval(
        image,
        diameter=_model_diameter,  # None = auto-detect
        min_size=MIN_AREA_PX,
        batch_size=cellpose_batch_size or _cellpose_patch_batch_size(),
    )
    props = sk_measure.regionprops(masks, intensity_image=image.astype(np.float32))
    return masks, props


def segment_cells_batch(images, cellpose_batch_size=None):
    """Segment a same-sized batch of grayscale tiles with one Cellpose call."""
    if not images:
        return []
    if len(images) == 1:
        masks, _ = segment_cells(images[0], cellpose_batch_size=cellpose_batch_size)
        return [masks]

    stacked = np.stack([img[..., np.newaxis] for img in images], axis=0)
    model = get_model()
    masks, _, _ = model.eval(
        stacked,
        diameter=_model_diameter,  # None = auto-detect
        min_size=MIN_AREA_PX,
        batch_size=cellpose_batch_size or _cellpose_patch_batch_size(),
    )
    return _normalize_mask_batch(masks)


def filter_shape(labeled, props):
    """Keep cells with good shape: round, solid. Cellpose handles size."""
    keep = set()
    for p in props:
        if p.solidity < MIN_SOLIDITY:
            continue
        circ = 4 * np.pi * p.area / (p.perimeter ** 2) if p.perimeter > 0 else 0
        if circ < MIN_CIRCULARITY:
            continue
        keep.add(p.label)

    filtered = np.where(np.isin(labeled, list(keep)), labeled, 0)
    return filtered, [p for p in props if p.label in keep]


def filter_isolated(labeled, props):
    """Keep cells whose nearest neighbor centroid is > NEIGHBOR_DISTANCE_PX away."""
    if len(props) < 2:
        return labeled, props

    centroids = np.array([p.centroid for p in props])
    tree = KDTree(centroids)
    dists, _ = tree.query(centroids, k=2)  # k=2: self + nearest other
    nearest = dists[:, 1]  # skip self (index 0)

    keep = {p.label for p, d in zip(props, nearest) if d > NEIGHBOR_DISTANCE_PX}
    filtered = np.where(np.isin(labeled, list(keep)), labeled, 0)
    return filtered, [p for p in props if p.label in keep]


def measure(image, labeled, props):
    """Measure area and IOD per Hardie et al. (2002).

    IOD = sum(log10(I_bg / I_pixel)) over object pixels.
    For brightfield: objects are cells (area measurement).
    For pmount: objects are nuclei (area + IOD for genome size).
    """
    bg_mask = labeled == 0
    i_bg = float(np.percentile(image[bg_mask], 95)) if bg_mask.any() else 255.0

    results = []
    for p in props:
        coords = p.coords
        pixel_vals = np.clip(image[coords[:, 0], coords[:, 1]].astype(np.float64), 1, None)
        od_per_pixel = np.log10(i_bg / pixel_vals)
        iod = float(np.sum(od_per_pixel))
        mean_od = float(np.mean(od_per_pixel))

        row = {
            "label": p.label,
            "area_px": p.area,
            "area_um2": p.area * PIXEL_AREA_UM2,
            "solidity": p.solidity,
            "circularity": 4 * np.pi * p.area / (p.perimeter ** 2) if p.perimeter > 0 else 0,
            "iod": iod,
            "mean_od": mean_od,
            "centroid_y": p.centroid[0],
            "centroid_x": p.centroid[1],
            "i_bg": i_bg,
        }
        results.append(row)
    return results


DEFAULT_CROP_SIZE = MAX_NUCLEUS_DIAMETER_PX * 2  # 400px


def _artifact_paths(output_dir, image_stem):
    root = Path(output_dir)
    return {
        "mask_dir": root / "masks" / image_stem,
        "tile_manifest": root / "tile_manifests" / f"{image_stem}_tile_manifest.csv",
        "debug_dir": root / "debug",
        "partial_measurements": root / f"{image_stem}_measurements.partial.csv",
    }


def _write_label_mask(label_image, mask_path):
    """Persist a label image with a stable integer dtype."""
    mask_path = Path(mask_path)
    mask_path.parent.mkdir(parents=True, exist_ok=True)
    tifffile.imwrite(str(mask_path), label_image.astype(np.uint32), compression="zlib")


def _write_tile_manifest(manifest_rows, manifest_path):
    """Write one manifest row per tile or full-image region."""
    manifest_path = Path(manifest_path)
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    if not manifest_rows:
        return
    fieldnames = [
        "tile_name", "tile_y0", "tile_x0", "tile_h", "tile_w",
        "tile_score", "tile_status", "mask_path", "raw_mask_path",
    ]
    with open(manifest_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(manifest_rows)


def _append_tile_manifest_rows(manifest_rows, manifest_path):
    """Append processed tile rows so interrupted runs can resume mid-image."""
    manifest_path = Path(manifest_path)
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    if not manifest_rows:
        return
    fieldnames = [
        "tile_name", "tile_y0", "tile_x0", "tile_h", "tile_w",
        "tile_score", "tile_status", "mask_path", "raw_mask_path",
    ]
    write_header = not manifest_path.exists() or manifest_path.stat().st_size == 0
    with open(manifest_path, "a", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        if write_header:
            writer.writeheader()
        writer.writerows(manifest_rows)


def _load_tile_manifest_rows(manifest_path):
    manifest_path = Path(manifest_path)
    if not manifest_path.exists() or manifest_path.stat().st_size == 0:
        return []
    with open(manifest_path, newline="") as f:
        return list(csv.DictReader(f))


def _append_rows(csv_path, rows):
    csv_path = Path(csv_path)
    if not rows:
        return
    csv_path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = list(rows[0].keys())
    write_header = not csv_path.exists() or csv_path.stat().st_size == 0
    with open(csv_path, "a", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        if write_header:
            writer.writeheader()
        writer.writerows(rows)


def _load_csv_rows(csv_path):
    csv_path = Path(csv_path)
    if not csv_path.exists() or csv_path.stat().st_size == 0:
        return []
    with open(csv_path, newline="") as f:
        return list(csv.DictReader(f))


def save_cell_crop(image, centroid_y, centroid_x, label, crop_dir, image_stem, crop_size=DEFAULT_CROP_SIZE):
    """Cut a square crop around a cell centroid and save as TIFF."""
    half = crop_size // 2
    cy, cx = int(round(centroid_y)), int(round(centroid_x))
    h, w = image.shape[:2]

    # Source region (may extend past image bounds)
    sy0, sy1 = cy - half, cy + half
    sx0, sx1 = cx - half, cx + half

    # Pad with background if needed
    pad_val = int(np.median(image))
    crop = np.full((crop_size, crop_size), pad_val, dtype=np.uint8)

    # Clamp to image bounds
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


def process_tile(tile, tile_y0, tile_x0):
    """Process one tile and return both measurements and label images."""
    raw_labeled, props = segment_cells(tile)
    return process_labeled_tile(tile, raw_labeled, tile_y0, tile_x0)


def process_labeled_tile(tile, raw_labeled, tile_y0, tile_x0):
    """Filter and measure one tile from a precomputed raw label image."""
    h, w = tile.shape[:2]
    props = sk_measure.regionprops(raw_labeled, intensity_image=tile.astype(np.float32))
    n_raw = len(props)
    labeled, props = filter_shape(raw_labeled, props)
    labeled, props = filter_isolated(labeled, props)

    # Discard objects near tile edges (may be partial)
    interior_props = [p for p in props
                      if EDGE_MARGIN < p.centroid[0] < h - EDGE_MARGIN
                      and EDGE_MARGIN < p.centroid[1] < w - EDGE_MARGIN]

    interior_labels = {p.label for p in interior_props}
    filtered_labeled = np.where(np.isin(labeled, list(interior_labels)), labeled, 0) if interior_labels else np.zeros_like(labeled)

    measurements = measure(tile, filtered_labeled, interior_props)
    for m in measurements:
        m["centroid_y"] += tile_y0
        m["centroid_x"] += tile_x0
        m["tile_name"] = f"tile_y{tile_y0:06d}_x{tile_x0:06d}.tiff"
        m["tile_y0"] = tile_y0
        m["tile_x0"] = tile_x0
        m["tile_h"] = h
        m["tile_w"] = w

    return measurements, n_raw, raw_labeled, filtered_labeled


def process_image(image_path, output_dir=None, debug=False, crops=True,
                  crop_size=DEFAULT_CROP_SIZE, selected_tiles=None,
                  tile_filter="auto", image_type="brightfield",
                  diameter_override=None, max_cells=None, pretrained_model=None):
    """Full pipeline for one image. Tiles automatically for large images.

    image_type: 'brightfield' = cell segmentation (auto diameter),
                'pmount' = nucleus segmentation (diameter=13, IOD for genome size).
    selected_tiles: optional list of [y0, x0, y1, x1] coords to process.
    tile_filter: 'green' = skip empty+borderline, 'auto' = skip empty only,
                 'yellow' = only borderline tiles, 'all' = no filtering.
    """
    global _model_diameter
    # Pmount defaults to fixed nucleus diameter; brightfield uses auto diameter.
    if diameter_override is not None:
        _model_diameter = diameter_override
    else:
        _model_diameter = PMOUNT_DIAMETER_PX if image_type == "pmount" else None
    get_model(pretrained_model=pretrained_model)

    image_path = Path(image_path)
    name = image_path.stem
    output_dir = Path(output_dir or PROJECT / "output" / "segmentation")
    output_dir.mkdir(parents=True, exist_ok=True)

    with tifffile.TiffFile(str(image_path)) as tif:
        if not tif.pages:
            raise ValueError(f"TIFF has no pages: {image_path}")
        page = tif.pages[0]
        img_h, img_w = page.shape[:2]

    if max_cells is None:
        max_cells = MAX_CELLS_PER_IMAGE

    needs_tiling = img_h > TILE_SIZE * 1.5 or img_w > TILE_SIZE * 1.5

    if not needs_tiling:
        return _process_small(image_path, name, output_dir, debug, crops, crop_size)
    else:
        return _process_tiled(image_path, name, img_h, img_w, output_dir, crops, crop_size,
                              selected_tiles=selected_tiles, tile_filter=tile_filter,
                              max_cells=max_cells)


def _process_small(image_path, name, output_dir, debug, crops=True, crop_size=DEFAULT_CROP_SIZE):
    """Direct processing for images that fit in RAM."""
    img = tifffile.imread(str(image_path))
    if img.ndim == 3:
        img = img[0]
    if img.dtype == np.uint16:
        img = (img >> 8).astype(np.uint8)

    print(f"  {name}: {img.shape[1]}x{img.shape[0]}")

    raw_labeled, props = segment_cells(img)
    n_raw = len(props)
    labeled, props = filter_shape(raw_labeled, props)
    n_shaped = len(props)
    labeled, props = filter_isolated(labeled, props)
    n_isolated = len(props)

    print(f"  {name}: {n_raw} raw -> {n_shaped} shaped -> {n_isolated} isolated")

    artifacts = _artifact_paths(output_dir, name)
    raw_mask_path = artifacts["mask_dir"] / f"{name}_raw_labels.tiff"
    mask_path = artifacts["mask_dir"] / f"{name}_filtered_labels.tiff"
    _write_label_mask(raw_labeled, raw_mask_path)
    _write_label_mask(labeled, mask_path)

    measurements = measure(img, labeled, props)
    tile_manifest_rows = [{
        "tile_name": f"{name}_full_image.tiff",
        "tile_y0": 0,
        "tile_x0": 0,
        "tile_h": int(img.shape[0]),
        "tile_w": int(img.shape[1]),
        "tile_score": 1.0,
        "tile_status": "processed",
        "mask_path": str(mask_path.resolve()),
        "raw_mask_path": str(raw_mask_path.resolve()),
    }]
    _write_tile_manifest(tile_manifest_rows, artifacts["tile_manifest"])

    if crops and measurements:
        crop_dir = output_dir / "crops" / name
        for m in measurements:
            fname = save_cell_crop(img, m["centroid_y"], m["centroid_x"],
                                   m["label"], crop_dir, name, crop_size)
            m["crop_filename"] = str(crop_dir / fname)
            m["crop_size_px"] = crop_size
            m["centroid_in_crop_y"] = crop_size // 2
            m["centroid_in_crop_x"] = crop_size // 2
        print(f"  {name}: saved {len(measurements)} crops to {crop_dir}")

    overlay_path = ""
    if debug:
        overlay_path = str(_save_overlay(img, labeled, measurements, artifacts["debug_dir"], name).resolve())

    for m in measurements:
        m["tile_name"] = f"{name}_full_image.tiff"
        m["tile_y0"] = 0
        m["tile_x0"] = 0
        m["tile_h"] = img.shape[0]
        m["tile_w"] = img.shape[1]
        m["tile_score"] = 1.0
        m["mask_path"] = str(mask_path.resolve())
        m["raw_mask_path"] = str(raw_mask_path.resolve())
        m["mask_label_id"] = m["label"]
        m["tile_manifest_path"] = str(artifacts["tile_manifest"].resolve())
        m["overlay_path"] = overlay_path

    _save_csv(measurements, output_dir, name)
    return measurements


def _process_tiled(image_path, name, img_h, img_w, output_dir, crops=True,
                   crop_size=DEFAULT_CROP_SIZE, selected_tiles=None,
                   tile_filter="auto", max_cells=0):
    """Zarr-based tiled processing for large images. Constant ~200MB RAM.

    Pre-scans all tiles, scores by sparseness (moderate std = isolated cells),
    then processes tiles in best-first order. Sparse tiles yield the most
    high-quality isolated cells per GPU-second.

    max_cells: stop after this many cells (0 = unlimited).
    selected_tiles: optional list of [y0, x0, y1, x1] to restrict processing.
    tile_filter: 'green' = high-content only, 'auto' = skip empties,
                 'yellow' = borderline only, 'all' = no filtering.
    """
    t0 = time.time()
    store = tifffile.imread(str(image_path), aszarr=True)
    z = zarr.open(store, mode='r')
    if z.ndim == 3:
        z = z[0]

    crop_dir = output_dir / "crops" / name if crops else None
    artifacts = _artifact_paths(output_dir, name)
    partial_measurement_path = artifacts["partial_measurements"]
    existing_measurements = _load_csv_rows(partial_measurement_path)
    existing_manifest_rows = _load_tile_manifest_rows(artifacts["tile_manifest"])
    processed_manifest_rows = {
        row["tile_name"]: row
        for row in existing_manifest_rows
        if row.get("tile_status") == "processed" and row.get("tile_name")
    }
    processed_tile_names = set(processed_manifest_rows)
    if processed_tile_names or existing_measurements:
        print(f"  {name}: resuming {len(processed_tile_names)} processed tiles "
              f"and {len(existing_measurements)} saved measurements")
    skipped_manifest_rows = []

    # Build tile coordinate list
    tile_coords = []
    for y0 in range(0, img_h, TILE_SIZE):
        for x0 in range(0, img_w, TILE_SIZE):
            y1 = min(y0 + TILE_SIZE, img_h)
            x1 = min(x0 + TILE_SIZE, img_w)
            if y1 - y0 >= EDGE_MARGIN * 2 and x1 - x0 >= EDGE_MARGIN * 2:
                tile_coords.append((y0, x0, y1, x1))

    # Filter to pre-selected tiles if provided
    total_tiles = len(tile_coords)
    if selected_tiles is not None:
        sel_set = {(t[0], t[1], t[2], t[3]) for t in selected_tiles}
        tile_coords = [t for t in tile_coords if t in sel_set]
        print(f"  {name}: {len(tile_coords)}/{total_tiles} tiles selected")

    # --- Phase 1: Pre-scan all tiles, classify and score by sparseness ---
    scored = []
    tiles_skipped = 0
    tiles_skipped_yellow = 0

    for y0, x0, y1, x1 in tile_coords:
        cy, cx = (y0 + y1) // 2, (x0 + x1) // 2
        half = min(128, (y1 - y0) // 2, (x1 - x0) // 2)
        probe = np.array(z[cy - half:cy + half, cx - half:cx + half])
        if probe.dtype == np.uint16:
            probe = (probe >> 8).astype(np.uint8)
        tmean, tstd = float(probe.mean()), float(probe.std())

        is_red = tmean > 230 or tmean < 10 or tstd < 5
        is_yellow = not is_red and (tstd < 15 or tmean > 200)

        tile_name = f"tile_y{y0:06d}_x{x0:06d}.tiff"
        if is_red:
            tiles_skipped += 1
            skipped_manifest_rows.append({
                "tile_name": tile_name,
                "tile_y0": y0,
                "tile_x0": x0,
                "tile_h": y1 - y0,
                "tile_w": x1 - x0,
                "tile_score": "",
                "tile_status": "skipped_empty",
                "mask_path": "",
                "raw_mask_path": "",
            })
            continue
        if tile_filter == "green" and is_yellow:
            tiles_skipped_yellow += 1
            skipped_manifest_rows.append({
                "tile_name": tile_name,
                "tile_y0": y0,
                "tile_x0": x0,
                "tile_h": y1 - y0,
                "tile_w": x1 - x0,
                "tile_score": "",
                "tile_status": "skipped_borderline",
                "mask_path": "",
                "raw_mask_path": "",
            })
            continue
        if tile_filter == "yellow" and not is_yellow:
            skipped_manifest_rows.append({
                "tile_name": tile_name,
                "tile_y0": y0,
                "tile_x0": x0,
                "tile_h": y1 - y0,
                "tile_w": x1 - x0,
                "tile_score": "",
                "tile_status": "skipped_filter",
                "mask_path": "",
                "raw_mask_path": "",
            })
            continue

        # Sparseness score: moderate std (~25) = isolated cells on clean bg.
        # High std = dense clusters (harder to isolate). Gaussian weight.
        score = float(np.exp(-((tstd - 18) ** 2) / (2 * 10 ** 2)))
        if is_yellow:
            score *= 0.5
        scored.append((score, y0, x0, y1, x1, tile_name))

    # Sort: sparsest tiles first (highest score = most isolated cells expected)
    scored.sort(reverse=True)

    skip_detail = f"{tiles_skipped} empty"
    if tiles_skipped_yellow:
        skip_detail += f", {tiles_skipped_yellow} borderline"
    print(f"  {name}: {img_w}x{img_h} ({len(scored)} tiles to process, "
          f"{skip_detail}, filter={tile_filter})")

    # --- Phase 2: Process tiles in sparseness-priority order ---
    measurements = list(existing_measurements)
    tiles_processed = len(processed_tile_names)
    crop_counter = 0
    pending_scored = [entry for entry in scored if entry[5] not in processed_tile_names]
    if crop_dir and measurements:
        crop_counter = len(measurements)
    batch_images = _tile_batch_images()
    budget_rows = {}

    for batch_start in range(0, len(pending_scored), batch_images):
        batch_entries = pending_scored[batch_start:batch_start + batch_images]
        tiles = []
        tile_infos = []
        for score, y0, x0, y1, x1, tile_name in batch_entries:
            tile = np.array(z[y0:y1, x0:x1])
            if tile.dtype == np.uint16:
                tile = (tile >> 8).astype(np.uint8)
            tiles.append(tile)
            tile_infos.append((score, y0, x0, y1, x1, tile_name))

        t_batch = time.time()
        raw_label_images = segment_cells_batch(tiles)
        batch_dt = time.time() - t_batch

        batch_measurements = []
        batch_manifest_rows = []
        for tile, raw_labeled, (score, y0, x0, y1, x1, tile_name) in zip(tiles, raw_label_images, tile_infos):
            tile_meas, _, raw_labeled, filtered_labeled = process_labeled_tile(tile, raw_labeled, y0, x0)

            raw_mask_path = artifacts["mask_dir"] / f"{tile_name[:-5]}__raw_labels.tiff"
            mask_path = artifacts["mask_dir"] / f"{tile_name[:-5]}__filtered_labels.tiff"
            _write_label_mask(raw_labeled, raw_mask_path)
            _write_label_mask(filtered_labeled, mask_path)
            manifest_row = {
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
            batch_manifest_rows.append(manifest_row)
            processed_manifest_rows[tile_name] = manifest_row

            for m in tile_meas:
                m["tile_score"] = score
                m["mask_path"] = str(mask_path.resolve())
                m["raw_mask_path"] = str(raw_mask_path.resolve())
                m["mask_label_id"] = m["label"]
                m["tile_manifest_path"] = str(artifacts["tile_manifest"].resolve())
                m["overlay_path"] = ""

            if crop_dir and tile_meas:
                for m in tile_meas:
                    crop_counter += 1
                    local_y = m["centroid_y"] - y0
                    local_x = m["centroid_x"] - x0
                    fname = save_cell_crop(tile, local_y, local_x,
                                           crop_counter, crop_dir, name, crop_size)
                    m["crop_filename"] = str(crop_dir / fname)
                    m["crop_size_px"] = crop_size
                    m["centroid_in_crop_y"] = crop_size // 2
                    m["centroid_in_crop_x"] = crop_size // 2

            batch_measurements.extend(tile_meas)
            measurements.extend(tile_meas)
            tiles_processed += 1

            if tiles_processed <= 3 or tiles_processed % 20 == 0:
                print(f"    tile {tiles_processed}/{len(scored)}: ({y0},{x0}) "
                      f"score={score:.2f} {len(tile_meas)} cells "
                      f"[batch {batch_dt:.1f}s/{len(tile_infos)} tiles] "
                      f"total={len(measurements)}", flush=True)

        _append_rows(partial_measurement_path, batch_measurements)
        _append_tile_manifest_rows(batch_manifest_rows, artifacts["tile_manifest"])

        if max_cells > 0 and len(measurements) >= max_cells:
            print(f"    Reached {max_cells} cell cap after "
                  f"{tiles_processed}/{len(scored)} tiles")
            for rem_score, rem_y0, rem_x0, rem_y1, rem_x1, rem_tile_name in pending_scored[batch_start + len(batch_entries):]:
                budget_rows[rem_tile_name] = {
                    "tile_name": rem_tile_name,
                    "tile_y0": rem_y0,
                    "tile_x0": rem_x0,
                    "tile_h": rem_y1 - rem_y0,
                    "tile_w": rem_x1 - rem_x0,
                    "tile_score": rem_score,
                    "tile_status": "not_processed_budget",
                    "mask_path": "",
                    "raw_mask_path": "",
                }
            break

    store.close()
    elapsed = time.time() - t0
    print(f"  {name}: {tiles_processed}/{len(scored)} tiles processed ({skip_detail}) "
          f"-> {len(measurements)} isolated [{elapsed:.0f}s]")
    if crop_dir:
        print(f"  {name}: saved {len(measurements)} crops to {crop_dir}")

    final_manifest_rows = list(skipped_manifest_rows)
    for _, _, _, _, _, tile_name in scored:
        if tile_name in processed_manifest_rows:
            final_manifest_rows.append(processed_manifest_rows[tile_name])
        elif tile_name in budget_rows:
            final_manifest_rows.append(budget_rows[tile_name])
    _write_tile_manifest(final_manifest_rows, artifacts["tile_manifest"])
    _save_csv(measurements, output_dir, name)
    if partial_measurement_path.exists():
        partial_measurement_path.unlink()
    return measurements


def _save_overlay(img, labeled, measurements, debug_dir, name):
    """Save debug overlay: red outlines, green centroids."""
    from PIL import Image as PILImage
    from scipy.ndimage import binary_dilation

    debug_dir = Path(debug_dir)
    debug_dir.mkdir(parents=True, exist_ok=True)

    rgb = np.stack([img, img, img], axis=-1)

    # Object outlines (red)
    outlines = segmentation.find_boundaries(labeled, mode='outer')
    outlines = binary_dilation(outlines, iterations=2)
    rgb[outlines] = [255, 0, 0]

    for m in measurements:
        y, x = int(m["centroid_y"]), int(m["centroid_x"])
        rgb[max(0, y-3):y+4, max(0, x-3):x+4] = [0, 255, 0]

    # Save as PNG (smaller, viewable)
    overlay_path = debug_dir / f"{name}_overlay.png"
    PILImage.fromarray(rgb).save(str(overlay_path))
    return overlay_path


def _save_csv(measurements, output_dir, name):
    """Save per-image measurements CSV."""
    if not measurements:
        return
    csv_path = output_dir / f"{name}_measurements.csv"
    with open(csv_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=measurements[0].keys())
        w.writeheader()
        w.writerows(measurements)
    print(f"  {name}: saved {len(measurements)} measurements to {csv_path.name}")


def _write_dataset_summary(all_results, output_dir):
    """Write per-species summary stats after batch run."""
    from collections import defaultdict
    species_data = defaultdict(list)
    for m in all_results:
        species_data[m.get("species", "unknown")].append(m)

    summary_path = output_dir / "dataset_summary.csv"
    rows = []
    for species, cells in sorted(species_data.items()):
        areas = [c["area_um2"] for c in cells]
        iods = [c["iod"] for c in cells]
        rows.append({
            "species": species,
            "count": len(cells),
            "mean_area_um2": f"{np.mean(areas):.2f}",
            "sd_area_um2": f"{np.std(areas):.2f}",
            "mean_iod": f"{np.mean(iods):.3f}",
            "iod_cv": f"{np.std(iods) / np.mean(iods):.3f}" if np.mean(iods) > 0 else "0",
        })
    rows.append({
        "species": "TOTAL",
        "count": len(all_results),
        "mean_area_um2": f"{np.mean([c['area_um2'] for c in all_results]):.2f}",
        "sd_area_um2": f"{np.std([c['area_um2'] for c in all_results]):.2f}",
        "mean_iod": f"{np.mean([c['iod'] for c in all_results]):.3f}",
        "iod_cv": f"{np.std([c['iod'] for c in all_results]) / np.mean([c['iod'] for c in all_results]):.3f}" if all_results else "0",
    })

    with open(summary_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=rows[0].keys())
        w.writeheader()
        w.writerows(rows)
    print(f"Dataset summary: {summary_path}")


def run_batch(output_dir, debug=False, crops=True, crop_size=DEFAULT_CROP_SIZE,
              diameter_override=None):
    """Process all images in analysis_set from master CSV."""
    master_path = PROJECT / "master_image_metadata.csv"
    if not master_path.exists():
        print(f"ERROR: {master_path} not found")
        sys.exit(1)

    with open(master_path) as f:
        rows = list(csv.DictReader(f))

    analysis = [r for r in rows if r.get("analysis_set") == "True"
                and r.get("image_type") in ("brightfield", "pmount")]
    print(f"Batch: {len(analysis)} images in analysis set")

    all_results = []
    for i, row in enumerate(analysis):
        img_type = row["image_type"]
        filename = row["filename"]
        img_path = PROJECT / "data" / img_type / filename

        if not img_path.exists():
            print(f"  SKIP [{i+1}/{len(analysis)}] {filename}: not found")
            continue

        print(f"\n[{i+1}/{len(analysis)}] {filename} ({row.get('species', '?')}, {img_type})")
        try:
            t0 = time.time()
            measurements = process_image(
                img_path,
                output_dir / img_type,
                debug,
                crops,
                crop_size,
                image_type=img_type,
                diameter_override=diameter_override,
            )
            for m in measurements:
                m["filename"] = filename
                m["slide_id"] = row.get("slide_id", "")
                m["specimen_id"] = row.get("specimen_id", "")
                m["species"] = row.get("species", "")
                m["image_type"] = img_type
            all_results.extend(measurements)
        except Exception as e:
            print(f"  ERROR: {e}")
            import traceback; traceback.print_exc()
            continue

    if all_results:
        combined_path = output_dir / "all_measurements.csv"
        with open(combined_path, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=all_results[0].keys())
            w.writeheader()
            w.writerows(all_results)
        print(f"\nTotal: {len(all_results)} nuclei from {len(analysis)} images")
        print(f"Saved: {combined_path}")
        _write_dataset_summary(all_results, output_dir)
    else:
        print("\nNo measurements generated")

    return all_results


def main():
    parser = argparse.ArgumentParser(description="Segment isolated RBCs")
    parser.add_argument("image", nargs="?", help="Single image to process")
    parser.add_argument("--output-dir", "-o", type=Path, default=None)
    parser.add_argument("--debug", "-d", action="store_true")
    parser.add_argument("--batch", action="store_true")
    parser.add_argument("--crops", action=argparse.BooleanOptionalAction, default=True,
                        help="Save individual cell crops (default: enabled)")
    parser.add_argument("--crop-size", type=int, default=DEFAULT_CROP_SIZE,
                        help=f"Crop side length in px (default: {DEFAULT_CROP_SIZE})")
    parser.add_argument("--gpu", action="store_true", help="Use GPU for Cellpose inference")
    parser.add_argument("--diameter", type=float, default=None,
                        help="Override Cellpose auto-diameter detection (pixels)")
    parser.add_argument("--image-type", choices=["brightfield", "pmount"],
                        default="brightfield",
                        help="Segmentation mode for single-image runs (default: brightfield)")
    args = parser.parse_args()

    # Initialize Cellpose model with CLI settings
    get_model(gpu=args.gpu, diameter=args.diameter)

    if args.image:
        results = process_image(args.image, args.output_dir, args.debug,
                                args.crops, args.crop_size,
                                image_type=args.image_type,
                                diameter_override=args.diameter)
        print(f"\n{len(results)} objects measured")
        if results:
            areas = [r["area_um2"] for r in results]
            iods = [r["iod"] for r in results]
            print(f"  Area: {np.mean(areas):.1f} +/- {np.std(areas):.1f} um^2")
            print(f"  IOD:  {np.mean(iods):.2f} +/- {np.std(iods):.2f}")
    elif args.batch:
        run_batch(args.output_dir or PROJECT / "output" / "segmentation",
                  args.debug, args.crops, args.crop_size,
                  diameter_override=args.diameter)
    else:
        parser.print_help()
        sys.exit(1)


if __name__ == "__main__":
    main()
