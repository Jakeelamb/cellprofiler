#!/usr/bin/env python3
"""Pure-Python nucleus IOD measurement — replaces ImageJ macro.

Hardie et al. (2002) protocol implemented in scipy/skimage:
  1. I_bg = 95th percentile of image (or cached value)
  2. Invert (dark nuclei -> bright)
  3. Gaussian blur sigma=4
  4. Otsu threshold -> binary
  5. Watershed (split touching nuclei)
  6. Connected components + size/circularity filter
  7. IOD = sum(log10(I_bg / I_pixel)) per nucleus on original image

Uses zarr random-access tiling (same as Cellpose pipeline) —
no disk I/O, no Java/ImageJ subprocess overhead. Tiles are scored
by sparseness and processed sparsest-first for fastest convergence.

Drop-in replacement for imagej_nucleus_iod.py: same process_image()
signature, same CSV columns, same save_csv().
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
import time
from collections import defaultdict
from pathlib import Path

import numpy as np
import tifffile
import zarr
from scipy import ndimage
from skimage.filters import gaussian, threshold_otsu
from skimage.measure import label, regionprops
from skimage.segmentation import watershed

PROJECT = Path(__file__).resolve().parent.parent

PIXEL_SIZE_UM = 0.12
PIXEL_AREA_UM2 = PIXEL_SIZE_UM ** 2  # 0.0144

# Hardie protocol defaults
BLUR_SIGMA = 4.0
MIN_AREA_PX = 500
MAX_AREA_PX = 5000
MIN_CIRCULARITY = 0.5

# Tiling
TILE_SIZE = 4096
EDGE_MARGIN = 200

CSV_COLUMNS = [
    "filename", "label", "area_px", "area_um2", "iod", "mean_od",
    "centroid_x", "centroid_y", "i_bg",
    "tile_name", "tile_y0", "tile_x0", "tile_height_px", "tile_width_px",
    "slide_id", "specimen_id", "species", "image_type",
]


def _normalize_uint8(arr: np.ndarray) -> np.ndarray:
    if arr.ndim == 3:
        arr = arr[0]
    if arr.dtype == np.uint16:
        arr = (arr >> 8).astype(np.uint8)
    elif arr.dtype != np.uint8:
        arr = np.clip(arr, 0, 255).astype(np.uint8)
    return arr


def classify_tile(probe: np.ndarray) -> tuple[str, float]:
    """Classify tile and return (class, sparseness_score)."""
    tmean, tstd = float(probe.mean()), float(probe.std())
    if tmean > 230 or tmean < 10 or tstd < 5:
        return "red", 0.0
    if tstd < 15 or tmean > 200:
        score = float(np.exp(-((tstd - 18) ** 2) / (2 * 10 ** 2))) * 0.5
        return "yellow", score
    score = float(np.exp(-((tstd - 18) ** 2) / (2 * 10 ** 2)))
    return "green", score


def segment_nuclei(tile: np.ndarray, manual_threshold: int | None = None) -> np.ndarray:
    """Hardie protocol: invert -> blur -> threshold -> watershed -> labeled.

    Returns labeled array (0=background, 1..N=nuclei).
    If manual_threshold is provided, uses it instead of Otsu.
    """
    # Invert (dark nuclei -> bright for thresholding)
    inverted = 255 - tile.astype(np.float32)

    # Gaussian blur
    blurred = gaussian(inverted, sigma=BLUR_SIGMA, preserve_range=True)

    # Threshold (manual override or Otsu)
    if manual_threshold is not None:
        thresh = manual_threshold
    else:
        try:
            thresh = threshold_otsu(blurred)
        except ValueError:
            return np.zeros(tile.shape[:2], dtype=np.int32)
    binary = blurred > thresh

    # Distance transform + watershed to split touching nuclei
    distance = ndimage.distance_transform_edt(binary)
    # Find local maxima as watershed seeds
    from skimage.feature import peak_local_max
    coords = peak_local_max(distance, min_distance=10, labels=binary)
    if len(coords) == 0:
        return label(binary)
    markers = np.zeros_like(binary, dtype=np.int32)
    for i, (y, x) in enumerate(coords, 1):
        markers[y, x] = i
    markers = ndimage.label(ndimage.binary_dilation(markers > 0, iterations=2))[0]
    labeled = watershed(-distance, markers, mask=binary)

    return labeled


def filter_nuclei(labeled: np.ndarray) -> list:
    """Filter by size and circularity, return regionprops list."""
    props = regionprops(labeled)
    kept = []
    for p in props:
        if p.area < MIN_AREA_PX or p.area > MAX_AREA_PX:
            continue
        circ = 4 * np.pi * p.area / (p.perimeter ** 2) if p.perimeter > 0 else 0
        if circ < MIN_CIRCULARITY:
            continue
        kept.append(p)
    return kept


def measure_iod(tile: np.ndarray, props: list, i_bg: float) -> list[dict]:
    """Measure IOD per nucleus on the original (uninverted) image."""
    results = []
    for p in props:
        coords = p.coords
        pixel_vals = np.clip(tile[coords[:, 0], coords[:, 1]].astype(np.float64), 1, None)
        od_per_pixel = np.log10(i_bg / pixel_vals)
        iod = float(np.sum(od_per_pixel))
        mean_od = float(np.mean(od_per_pixel))
        area_um2 = p.area * PIXEL_AREA_UM2

        results.append({
            "label": p.label,
            "area_px": p.area,
            "area_um2": round(area_um2, 4),
            "iod": round(iod, 6),
            "mean_od": round(mean_od, 6),
            "centroid_y": round(p.centroid[0], 2),
            "centroid_x": round(p.centroid[1], 2),
            "i_bg": i_bg,
        })
    return results


def process_tile(tile: np.ndarray, i_bg: float, tile_y0: int, tile_x0: int,
                 manual_threshold: int | None = None) -> list[dict]:
    """Segment + filter + measure one tile. Discard edge nuclei."""
    h, w = tile.shape[:2]
    labeled = segment_nuclei(tile, manual_threshold=manual_threshold)
    props = filter_nuclei(labeled)

    # Discard nuclei near tile edges
    interior = [p for p in props
                if EDGE_MARGIN < p.centroid[0] < h - EDGE_MARGIN
                and EDGE_MARGIN < p.centroid[1] < w - EDGE_MARGIN]

    measurements = measure_iod(tile, interior, i_bg)
    for m in measurements:
        m["centroid_y"] = round(m["centroid_y"] + tile_y0, 2)
        m["centroid_x"] = round(m["centroid_x"] + tile_x0, 2)
        m["tile_name"] = f"tile_y{tile_y0:06d}_x{tile_x0:06d}.tiff"
        m["tile_y0"] = tile_y0
        m["tile_x0"] = tile_x0
        m["tile_height_px"] = h
        m["tile_width_px"] = w
    return measurements


def compute_image_ibg(z, img_h: int, img_w: int,
                      sample_grid: int = 8, sample_box: int = 256) -> int:
    """95th percentile background from a grid of sample boxes (in-memory, no vips)."""
    box_h = min(sample_box, img_h)
    box_w = min(sample_box, img_w)
    y_positions = np.linspace(0, max(0, img_h - box_h), num=sample_grid, dtype=int)
    x_positions = np.linspace(0, max(0, img_w - box_w), num=sample_grid, dtype=int)

    hist = np.zeros(256, dtype=np.int64)
    total = 0
    for y0 in y_positions:
        for x0 in x_positions:
            chunk = np.array(z[int(y0):int(y0) + box_h, int(x0):int(x0) + box_w])
            chunk = _normalize_uint8(chunk)
            hist += np.bincount(chunk.ravel(), minlength=256)
            total += chunk.size

    if total == 0:
        return 255
    target = total * 0.95
    cumulative = np.cumsum(hist)
    return int(np.searchsorted(cumulative, target))


def process_image(image_path, tile_filter="auto", cached_i_bg=None, artifact_dir=None,
                  manual_threshold=None):
    """Process one image: zarr tiling, sparse-first scoring, Python segmentation.

    Drop-in replacement for imagej_nucleus_iod.process_image().
    If manual_threshold is provided, uses it instead of per-tile Otsu.
    """
    image_path = Path(image_path)
    name = image_path.stem

    store = tifffile.imread(str(image_path), aszarr=True)
    z = zarr.open(store, mode="r")
    if isinstance(z, zarr.Group):
        z = z["0"]
    if z.ndim == 3:
        z = z[0]
    img_h, img_w = z.shape[:2]

    # Background
    if cached_i_bg is not None:
        i_bg = cached_i_bg
    else:
        i_bg = compute_image_ibg(z, img_h, img_w)

    needs_tiling = img_h > TILE_SIZE * 1.5 or img_w > TILE_SIZE * 1.5

    if not needs_tiling:
        # Small image: process directly
        tile = _normalize_uint8(np.array(z[:, :]))
        store.close()
        th_label = f"thresh={manual_threshold}" if manual_threshold is not None else "Otsu"
        print(f"  {name}: {img_w}x{img_h} (direct, I_bg={i_bg}, {th_label})")
        labeled = segment_nuclei(tile, manual_threshold=manual_threshold)
        props = filter_nuclei(labeled)
        measurements = measure_iod(tile, props, i_bg)
        for m in measurements:
            m["tile_name"] = image_path.name
            m["tile_y0"] = 0
            m["tile_x0"] = 0
            m["tile_height_px"] = img_h
            m["tile_width_px"] = img_w
            m["filename"] = image_path.name
        print(f"  {name}: {len(measurements)} nuclei")
        _save_artifacts(artifact_dir, image_path, None, img_h, img_w)
        return measurements

    # Large image: zarr tiled processing with sparse-first scoring
    t0 = time.time()

    # Build and score tiles
    tile_coords = []
    for y0 in range(0, img_h, TILE_SIZE):
        for x0 in range(0, img_w, TILE_SIZE):
            y1 = min(y0 + TILE_SIZE, img_h)
            x1 = min(x0 + TILE_SIZE, img_w)
            if y1 - y0 >= EDGE_MARGIN * 2 and x1 - x0 >= EDGE_MARGIN * 2:
                tile_coords.append((y0, x0, y1, x1))

    scored = []
    skipped = 0
    for y0, x0, y1, x1 in tile_coords:
        cy, cx = (y0 + y1) // 2, (x0 + x1) // 2
        half = min(128, (y1 - y0) // 2, (x1 - x0) // 2)
        probe = _normalize_uint8(np.array(z[cy - half:cy + half, cx - half:cx + half]))
        cls, score = classify_tile(probe)

        if cls == "red":
            skipped += 1
            continue
        if tile_filter == "green" and cls == "yellow":
            skipped += 1
            continue
        scored.append((score, y0, x0, y1, x1))

    # Sort: sparsest first
    scored.sort(reverse=True)
    th_label = f"thresh={manual_threshold}" if manual_threshold is not None else "Otsu"
    print(f"  {name}: {img_w}x{img_h} (I_bg={i_bg}, {th_label}, {len(scored)} tiles, {skipped} skipped)")

    # Process tiles
    all_measurements = []
    tile_infos = []
    for idx, (score, y0, x0, y1, x1) in enumerate(scored):
        tile = _normalize_uint8(np.array(z[y0:y1, x0:x1]))
        t_tile = time.time()
        tile_meas = process_tile(tile, i_bg, y0, x0, manual_threshold=manual_threshold)
        dt = time.time() - t_tile
        all_measurements.extend(tile_meas)

        tile_infos.append({
            "tile_name": f"tile_y{y0:06d}_x{x0:06d}.tiff",
            "y0": y0, "x0": x0,
            "height_px": y1 - y0, "width_px": x1 - x0,
            "tile_class": "green",
        })

        if idx < 3 or (idx + 1) % 20 == 0:
            print(f"    tile {idx+1}/{len(scored)}: ({y0},{x0}) "
                  f"score={score:.2f} {len(tile_meas)} nuclei [{dt:.1f}s] "
                  f"total={len(all_measurements)}", flush=True)

    store.close()
    elapsed = time.time() - t0

    # Set filename on all rows
    for m in all_measurements:
        m["filename"] = image_path.name

    print(f"  {name}: {len(all_measurements)} nuclei after edge filtering [{elapsed:.0f}s]")

    _save_artifacts(artifact_dir, image_path, tile_infos, img_h, img_w)
    return all_measurements


def _save_artifacts(artifact_dir, image_path, tile_infos, img_h, img_w):
    """Write tile manifest to artifact dir if provided."""
    if artifact_dir is None:
        return
    artifact_dir = Path(artifact_dir)
    artifact_dir.mkdir(parents=True, exist_ok=True)
    if tile_infos:
        manifest_path = artifact_dir / "tile_manifest.csv"
        fieldnames = [
            "source_image", "image_height_px", "image_width_px",
            "tile_name", "y0", "x0", "height_px", "width_px", "tile_class",
        ]
        with open(manifest_path, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=fieldnames)
            w.writeheader()
            for t in tile_infos:
                w.writerow({
                    "source_image": image_path.name,
                    "image_height_px": img_h,
                    "image_width_px": img_w,
                    **t,
                })


def load_background_cache(csv_path):
    """Load per-image background estimates from an existing measurement CSV."""
    if csv_path is None:
        return {}
    csv_path = Path(csv_path)
    if not csv_path.exists():
        return {}
    vals = defaultdict(list)
    with open(csv_path) as f:
        for row in csv.DictReader(f):
            fname = row.get("filename")
            i_bg = row.get("i_bg")
            if not fname or i_bg in (None, ""):
                continue
            try:
                vals[fname].append(float(i_bg))
            except ValueError:
                continue
    return {fname: int(round(float(np.median(items)))) for fname, items in vals.items() if items}


def save_csv(measurements, csv_path, metadata=None):
    """Append measurements to output CSV."""
    if not measurements:
        return
    for m in measurements:
        if metadata:
            m.update(metadata)
    write_header = not csv_path.exists() or csv_path.stat().st_size == 0
    with open(csv_path, "a", newline="") as f:
        w = csv.DictWriter(f, fieldnames=CSV_COLUMNS, extrasaction="ignore")
        if write_header:
            w.writeheader()
        w.writerows(measurements)


def main():
    parser = argparse.ArgumentParser(description="Python nucleus IOD (Hardie et al. 2002)")
    parser.add_argument("image", nargs="?", help="Single image path")
    parser.add_argument("--batch", action="store_true")
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--image-type", choices=["brightfield", "pmount", "both"],
                        default="both")
    parser.add_argument("--limit", type=int, default=0)
    parser.add_argument("--tile-filter", choices=["green", "auto", "all"], default="auto")
    parser.add_argument("--output-dir", "-o", type=Path,
                        default=PROJECT / "output" / "nucleus_iod")
    parser.add_argument("--artifact-dir", type=Path)
    parser.add_argument("--background-cache", type=Path)
    args = parser.parse_args()

    if args.image:
        bg_cache = load_background_cache(args.background_cache)
        img_name = Path(args.image).name
        rows = process_image(
            args.image, tile_filter=args.tile_filter,
            cached_i_bg=bg_cache.get(img_name),
            artifact_dir=(args.artifact_dir / Path(args.image).stem) if args.artifact_dir else None,
        )
        args.output_dir.mkdir(parents=True, exist_ok=True)
        name = Path(args.image).stem
        csv_path = args.output_dir / f"{name}_nucleus_iod.csv"
        save_csv(rows, csv_path)
        print(f"\n{len(rows)} nuclei -> {csv_path}")
        if rows:
            iods = [float(r["iod"]) for r in rows]
            print(f"  IOD: {np.mean(iods):.2f} +/- {np.std(iods):.2f}")
        return

    if not args.batch:
        parser.print_help()
        sys.exit(1)

    # Batch mode
    master = PROJECT / "master_image_metadata.csv"
    if not master.exists():
        print(f"ERROR: {master} not found")
        sys.exit(1)

    with open(master) as f:
        all_rows = list(csv.DictReader(f))

    types = ["brightfield", "pmount"] if args.image_type == "both" else [args.image_type]
    jobs = []
    for r in all_rows:
        if r.get("analysis_set") != "True":
            continue
        if r.get("image_type") not in types:
            continue
        img_path = PROJECT / "data" / r["image_type"] / r["filename"]
        if not img_path.exists():
            print(f"  SKIP {r['filename']}: not found")
            continue
        jobs.append(r | {"path": str(img_path)})

    if args.limit > 0:
        jobs = jobs[:args.limit]

    bg_cache = load_background_cache(args.background_cache)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    csv_path = args.output_dir / "nucleus_iod_measurements.csv"

    completed = set()
    if args.resume and csv_path.exists():
        with open(csv_path) as f:
            for row in csv.DictReader(f):
                completed.add(row.get("filename", ""))
        print(f"Resuming: {len(completed)} images already done")

    remaining = [j for j in jobs if j["filename"] not in completed]
    print(f"Batch: {len(remaining)}/{len(jobs)} images ({', '.join(types)})")

    t_start = time.time()
    total_nuclei = 0

    for i, job in enumerate(remaining):
        print(f"\n[{len(completed) + i + 1}/{len(jobs)}] {job['filename']} "
              f"({job.get('species', '?')}, {job['image_type']})")
        try:
            rows = process_image(
                job["path"], tile_filter=args.tile_filter,
                cached_i_bg=bg_cache.get(job["filename"]),
                artifact_dir=(args.artifact_dir / Path(job["filename"]).stem) if args.artifact_dir else None,
            )
            save_csv(rows, csv_path, metadata={
                "slide_id": job.get("slide_id", ""),
                "specimen_id": job.get("specimen_id", ""),
                "species": job.get("species", ""),
                "image_type": job["image_type"],
            })
            total_nuclei += len(rows)
            completed.add(job["filename"])
        except Exception as e:
            print(f"  ERROR: {e}")
            import traceback
            traceback.print_exc()

    elapsed = time.time() - t_start
    print(f"\nDone: {total_nuclei} nuclei from {len(remaining)} images [{elapsed:.0f}s]")


if __name__ == "__main__":
    main()
