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
    python scripts/segment_cells.py <image.tiff> [--output-dir DIR] [--debug]
    python scripts/segment_cells.py --batch  # process all from master CSV
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
from skimage import measure, segmentation

PROJECT = Path(__file__).resolve().parent.parent
PIXEL_SIZE_UM = 0.12
PIXEL_AREA_UM2 = PIXEL_SIZE_UM ** 2  # 0.0144 um^2

# Segmentation params
MAX_NUCLEUS_DIAMETER_PX = 200
MIN_AREA_PX = int(np.pi * (30 / 2) ** 2 * 0.5)  # ~353

# Shape filters
MIN_SOLIDITY = 0.7
MIN_CIRCULARITY = 0.4  # 4*pi*area/perimeter^2

# Isolation: reject if nearest neighbor centroid within this many pixels
NEIGHBOR_DISTANCE_PX = 30  # ~3.6um at 0.12um/px

# Tiling
TILE_SIZE = 4096
EDGE_MARGIN = MAX_NUCLEUS_DIAMETER_PX  # discard objects near tile edges
MAX_CELLS_PER_IMAGE = 500  # stop after finding this many isolated cells

# Cellpose model (lazy-loaded singleton)
_model = None
_model_gpu = False
_model_diameter = None


def get_model(gpu=False, diameter=None):
    """Get or create Cellpose model. Reuses across calls."""
    global _model, _model_gpu, _model_diameter
    if _model is None or gpu != _model_gpu:
        _model = models.CellposeModel(gpu=gpu)
        _model_gpu = gpu
    _model_diameter = diameter
    return _model


def segment_nuclei(image):
    """Segment cells using Cellpose.

    Returns (labeled_array, regionprops_list).
    """
    model = get_model()
    masks, _, _ = model.eval(
        image,
        diameter=_model_diameter,  # None = auto-detect
        min_size=MIN_AREA_PX,
    )
    props = measure.regionprops(masks, intensity_image=image.astype(np.float32))
    return masks, props


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


_nuclei_model = None


def segment_nuclei_within_cells(image, cell_labeled):
    """Second-pass Cellpose on raw image with diameter=13 to find nuclei."""
    global _nuclei_model
    if _nuclei_model is None:
        _nuclei_model = models.CellposeModel(gpu=_model_gpu)
    nuc_masks, _, _ = _nuclei_model.eval(image, diameter=13, min_size=50)
    # Assign each nucleus to its parent cell by max overlap
    nuclei = np.zeros_like(cell_labeled)
    for nuc_label in np.unique(nuc_masks):
        if nuc_label == 0:
            continue
        nuc_mask = nuc_masks == nuc_label
        cell_vals = cell_labeled[nuc_mask]
        cell_vals = cell_vals[cell_vals > 0]
        if len(cell_vals) == 0:
            continue
        parent = np.bincount(cell_vals).argmax()
        nuclei[nuc_mask & (cell_labeled == parent)] = parent
    return nuclei


def measure_iod(image, labeled, props):
    """Calculate Integrated Optical Density per Hardie et al. (2002).

    IOD = sum(log10(I_bg / I_pixel)) over nucleus pixels.
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

        results.append({
            "label": p.label,
            "area_px": p.area,
            "area_um2": p.area * PIXEL_AREA_UM2,
            "centroid_y": p.centroid[0],
            "centroid_x": p.centroid[1],
            "solidity": p.solidity,
            "circularity": 4 * np.pi * p.area / (p.perimeter ** 2) if p.perimeter > 0 else 0,
            "mean_intensity": float(p.intensity_mean),
            "iod": iod,
            "mean_od": mean_od,
            "i_bg": i_bg,
        })
    return results


DEFAULT_CROP_SIZE = MAX_NUCLEUS_DIAMETER_PX * 2  # 400px


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
    """Process one tile, discard edge objects, return measurements in global coords."""
    h, w = tile.shape[:2]

    raw_labeled, props = segment_nuclei(tile)
    n_raw = len(props)
    labeled, props = filter_shape(raw_labeled, props)
    labeled, props = filter_isolated(labeled, props)

    # Discard objects near tile edges (may be partial)
    interior_props = [p for p in props
                      if EDGE_MARGIN < p.centroid[0] < h - EDGE_MARGIN
                      and EDGE_MARGIN < p.centroid[1] < w - EDGE_MARGIN]

    interior_labels = {p.label for p in interior_props}
    filtered_labeled = np.where(np.isin(labeled, list(interior_labels)), labeled, 0) if interior_labels else np.zeros_like(labeled)

    measurements = measure_iod(tile, filtered_labeled, interior_props)
    for m in measurements:
        m["centroid_y"] += tile_y0
        m["centroid_x"] += tile_x0

    return measurements, n_raw


def process_image(image_path, output_dir=None, debug=False, crops=True, crop_size=DEFAULT_CROP_SIZE):
    """Full pipeline for one image. Tiles automatically for large images."""
    image_path = Path(image_path)
    name = image_path.stem
    output_dir = Path(output_dir or PROJECT / "output" / "segmentation")
    output_dir.mkdir(parents=True, exist_ok=True)

    with tifffile.TiffFile(str(image_path)) as tif:
        page = tif.pages[0]
        img_h, img_w = page.shape[:2]

    needs_tiling = img_h > TILE_SIZE * 1.5 or img_w > TILE_SIZE * 1.5

    if not needs_tiling:
        return _process_small(image_path, name, output_dir, debug, crops, crop_size)
    else:
        return _process_tiled(image_path, name, img_h, img_w, output_dir, crops, crop_size)


def _process_small(image_path, name, output_dir, debug, crops=True, crop_size=DEFAULT_CROP_SIZE):
    """Direct processing for images that fit in RAM."""
    img = tifffile.imread(str(image_path))
    if img.ndim == 3:
        img = img[0]
    if img.dtype == np.uint16:
        img = (img >> 8).astype(np.uint8)

    print(f"  {name}: {img.shape[1]}x{img.shape[0]}")

    raw_labeled, props = segment_nuclei(img)
    n_raw = len(props)
    labeled, props = filter_shape(raw_labeled, props)
    n_shaped = len(props)
    labeled, props = filter_isolated(labeled, props)
    n_isolated = len(props)

    print(f"  {name}: {n_raw} raw -> {n_shaped} shaped -> {n_isolated} isolated")
    measurements = measure_iod(img, labeled, props)

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

    if debug:
        _save_overlay(img, labeled, measurements, output_dir / "debug", name)

    _save_csv(measurements, output_dir, name)
    return measurements


def _process_tiled(image_path, name, img_h, img_w, output_dir, crops=True, crop_size=DEFAULT_CROP_SIZE):
    """Zarr-based tiled processing for large images. Constant ~200MB RAM.

    Shuffles tile order for spatial sampling and stops after MAX_CELLS_PER_IMAGE.
    """
    import random

    t0 = time.time()
    store = tifffile.imread(str(image_path), aszarr=True)
    z = zarr.open(store, mode='r')
    if z.ndim == 3:
        z = z[0]

    crop_dir = output_dir / "crops" / name if crops else None

    # Build tile coordinate list, shuffle for spatial sampling
    tiles = []
    for y0 in range(0, img_h, TILE_SIZE):
        for x0 in range(0, img_w, TILE_SIZE):
            y1 = min(y0 + TILE_SIZE, img_h)
            x1 = min(x0 + TILE_SIZE, img_w)
            if y1 - y0 >= EDGE_MARGIN * 2 and x1 - x0 >= EDGE_MARGIN * 2:
                tiles.append((y0, x0, y1, x1))

    random.seed(42)  # reproducible
    random.shuffle(tiles)
    print(f"  {name}: {img_w}x{img_h} ({len(tiles)} tiles)")

    measurements = []
    tiles_processed = 0
    tiles_skipped = 0
    crop_counter = 0

    for y0, x0, y1, x1 in tiles:
        tile = np.array(z[y0:y1, x0:x1])
        if tile.dtype == np.uint16:
            tile = (tile >> 8).astype(np.uint8)

        tmean = tile.mean()
        tstd = tile.std()
        # Skip empty/uniform tiles
        if tmean > 230 or tmean < 10 or tstd < 5:
            tiles_skipped += 1
            continue

        t_tile = time.time()
        tile_meas, _ = process_tile(tile, y0, x0)
        dt = time.time() - t_tile

        # Save crops from tile data (using local tile coordinates)
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

        measurements.extend(tile_meas)
        tiles_processed += 1

        if tiles_processed <= 3 or tiles_processed % 20 == 0:
            elapsed = time.time() - t0
            print(f"    tile {tiles_processed}: ({y0},{x0}) "
                  f"{len(tile_meas)} cells [{dt:.1f}s] total={len(measurements)}", flush=True)

        if len(measurements) >= MAX_CELLS_PER_IMAGE:
            print(f"    Reached {MAX_CELLS_PER_IMAGE} cells, stopping early")
            break

    store.close()
    elapsed = time.time() - t0
    print(f"  {name}: {tiles_processed} tiles ({tiles_skipped} empty) "
          f"-> {len(measurements)} isolated [{elapsed:.0f}s]")
    if crop_dir:
        print(f"  {name}: saved {len(measurements)} crops to {crop_dir}")

    _save_csv(measurements, output_dir, name)
    return measurements


def _save_overlay(img, labeled, measurements, debug_dir, name):
    """Save debug overlay: red cell outlines, cyan nucleus outlines, green centroids."""
    from PIL import Image as PILImage
    from scipy.ndimage import binary_dilation

    debug_dir = Path(debug_dir)
    debug_dir.mkdir(parents=True, exist_ok=True)

    rgb = np.stack([img, img, img], axis=-1)

    # Cell outlines (red)
    outlines = segmentation.find_boundaries(labeled, mode='outer')
    outlines = binary_dilation(outlines, iterations=2)
    rgb[outlines] = [255, 0, 0]

    # Nucleus outlines (cyan)
    nuclei = segment_nuclei_within_cells(img, labeled)
    nuc_outlines = segmentation.find_boundaries(nuclei, mode='outer')
    nuc_outlines = binary_dilation(nuc_outlines, iterations=2)
    rgb[nuc_outlines] = [0, 255, 255]

    for m in measurements:
        y, x = int(m["centroid_y"]), int(m["centroid_x"])
        rgb[max(0, y-3):y+4, max(0, x-3):x+4] = [0, 255, 0]

    # Save as PNG (smaller, viewable)
    PILImage.fromarray(rgb).save(str(debug_dir / f"{name}_overlay.png"))


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
            "cell_count": len(cells),
            "mean_area_um2": f"{np.mean(areas):.2f}",
            "std_area_um2": f"{np.std(areas):.2f}",
            "mean_iod": f"{np.mean(iods):.3f}",
            "iod_cv": f"{np.std(iods) / np.mean(iods):.3f}" if np.mean(iods) > 0 else "0",
        })
    rows.append({
        "species": "TOTAL",
        "cell_count": len(all_results),
        "mean_area_um2": f"{np.mean([c['area_um2'] for c in all_results]):.2f}",
        "std_area_um2": f"{np.std([c['area_um2'] for c in all_results]):.2f}",
        "mean_iod": f"{np.mean([c['iod'] for c in all_results]):.3f}",
        "iod_cv": f"{np.std([c['iod'] for c in all_results]) / np.mean([c['iod'] for c in all_results]):.3f}" if all_results else "0",
    })

    with open(summary_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=rows[0].keys())
        w.writeheader()
        w.writerows(rows)
    print(f"Dataset summary: {summary_path}")


def run_batch(output_dir, debug=False, crops=True, crop_size=DEFAULT_CROP_SIZE):
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
            measurements = process_image(img_path, output_dir / img_type, debug, crops, crop_size)
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
    args = parser.parse_args()

    # Initialize Cellpose model with CLI settings
    get_model(gpu=args.gpu, diameter=args.diameter)

    if args.image:
        results = process_image(args.image, args.output_dir, args.debug,
                                args.crops, args.crop_size)
        print(f"\n{len(results)} nuclei measured")
        if results:
            areas = [r["area_um2"] for r in results]
            iods = [r["iod"] for r in results]
            print(f"  Area: {np.mean(areas):.1f} +/- {np.std(areas):.1f} um^2")
            print(f"  IOD:  {np.mean(iods):.2f} +/- {np.std(iods):.2f}")
    elif args.batch:
        run_batch(args.output_dir or PROJECT / "output" / "segmentation",
                  args.debug, args.crops, args.crop_size)
    else:
        parser.print_help()
        sys.exit(1)


if __name__ == "__main__":
    main()
