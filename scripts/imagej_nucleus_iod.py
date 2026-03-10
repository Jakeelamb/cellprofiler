#!/usr/bin/env python3
"""Run ImageJ headless for nucleus detection + IOD measurement.

Wraps nucleus_iod.ijm macro. Works on both brightfield and pmount images.
For whole-slide images, pre-tiles into a temp directory then runs ImageJ
on each tile, merging results with edge deduplication.

Hardie et al. (2002) protocol:
  Invert -> Gaussian blur -> Otsu -> Watershed -> Particle analysis -> IOD

Usage:
    python scripts/imagej_nucleus_iod.py data/brightfield/image.tiff
    python scripts/imagej_nucleus_iod.py --batch
    python scripts/imagej_nucleus_iod.py --batch --image-type pmount --limit 5
"""

import argparse
import csv
import json
import os
import shutil
import subprocess
import sys
import tempfile
import time
from pathlib import Path

import numpy as np
import tifffile
import zarr
from collections import defaultdict

PROJECT = Path(__file__).resolve().parent.parent
MACRO_PATH = PROJECT / "scripts" / "nucleus_iod.ijm"
DEFAULT_BACKGROUND_CACHE = PROJECT / "output" / "nucleus_iod" / "nucleus_iod_measurements.csv"

PIXEL_SIZE_UM = 0.12

# Hardie protocol defaults
BLUR_SIGMA = 4.0
MIN_AREA_PX = 500
MAX_AREA_PX = 5000
MIN_CIRCULARITY = 0.5

# Tiling
TILE_SIZE = 4096
EDGE_MARGIN = 200

# ImageJ binary
IMAGEJ_CMD = "imagej"
IMAGEJ_MEMORY_MB = 4096
VIPS_CMD = shutil.which("vips")

CSV_COLUMNS = [
    "filename", "label", "area_px", "area_um2", "iod", "mean_od",
    "centroid_x", "centroid_y", "i_bg",
    "tile_name", "tile_y0", "tile_x0", "tile_height_px", "tile_width_px",
    "mask_path", "roi_zip_path", "tile_manifest_path", "raw_imagej_results_path",
    "slide_id", "specimen_id", "species", "image_type",
]


def classify_tile(probe):
    """Fast tile classification."""
    tmean, tstd = float(probe.mean()), float(probe.std())
    if tmean > 230 or tmean < 10 or tstd < 5:
        return "red"
    if tstd < 15 or tmean > 200:
        return "yellow"
    return "green"


def _normalize_grayscale_uint8(arr):
    """Normalize extracted tile/region arrays to 2D uint8 grayscale."""
    if arr.ndim == 3:
        arr = arr[0]
    if arr.dtype == np.uint16:
        arr = (arr >> 8).astype(np.uint8)
    elif arr.dtype != np.uint8:
        arr = np.clip(arr, 0, 255).astype(np.uint8)
    return arr


def extract_region_vips(image_path, x0, y0, width, height):
    """Extract one region with libvips and return it as a uint8 numpy array."""
    with tempfile.NamedTemporaryFile(prefix="vips_crop_", suffix=".tif", delete=False) as tmp:
        tmp_path = Path(tmp.name)

    try:
        cmd = [
            VIPS_CMD, "crop",
            str(image_path), str(tmp_path),
            str(int(x0)), str(int(y0)), str(int(width)), str(int(height)),
        ]
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
        if result.returncode != 0:
            raise RuntimeError(f"vips crop failed for {image_path}: {result.stderr[:400]}")
        arr = tifffile.imread(str(tmp_path))
        return _normalize_grayscale_uint8(arr)
    finally:
        try:
            tmp_path.unlink(missing_ok=True)
        except Exception:
            pass


def extract_tiles_tifffile(image_path, tile_dir, tile_filter="auto"):
    """Extract tiles from a whole-slide image into tile_dir as individual TIFFs.

    Returns list of tile metadata dicts.
    """
    tile_dir.mkdir(parents=True, exist_ok=True)
    store = tifffile.imread(str(image_path), aszarr=True)
    z = zarr.open(store, mode='r')
    if isinstance(z, zarr.Group):
        z = z['0']
    if z.ndim == 3:
        z = z[0]

    img_h, img_w = z.shape[:2]
    tiles = []
    skipped = 0

    for y0 in range(0, img_h, TILE_SIZE):
        for x0 in range(0, img_w, TILE_SIZE):
            y1 = min(y0 + TILE_SIZE, img_h)
            x1 = min(x0 + TILE_SIZE, img_w)
            if y1 - y0 < EDGE_MARGIN * 2 or x1 - x0 < EDGE_MARGIN * 2:
                continue

            # Fast classification
            cy, cx = (y0 + y1) // 2, (x0 + x1) // 2
            half = min(128, (y1 - y0) // 2, (x1 - x0) // 2)
            probe = np.array(z[cy - half:cy + half, cx - half:cx + half])
            if probe.dtype == np.uint16:
                probe = (probe >> 8).astype(np.uint8)

            cls = classify_tile(probe)
            if cls == "red":
                skipped += 1
                continue
            if tile_filter == "green" and cls == "yellow":
                skipped += 1
                continue

            # Write tile
            tile = np.array(z[y0:y1, x0:x1])
            if tile.dtype == np.uint16:
                tile = (tile >> 8).astype(np.uint8)

            tile_name = f"tile_y{y0:06d}_x{x0:06d}.tiff"
            tifffile.imwrite(str(tile_dir / tile_name), tile)
            tiles.append({
                "tile_name": tile_name,
                "y0": y0,
                "x0": x0,
                "height_px": y1 - y0,
                "width_px": x1 - x0,
                "tile_class": cls,
            })

    store.close()
    print(f"  Extracted {len(tiles)} tiles ({skipped} skipped)")
    return tiles, img_h, img_w


def extract_tiles_vips(image_path, tile_dir, tile_filter="auto", img_h=None, img_w=None):
    """Extract tiles using libvips crop, which is reliable on large tiled TIFFs."""
    tile_dir.mkdir(parents=True, exist_ok=True)
    if img_h is None or img_w is None:
        with tifffile.TiffFile(str(image_path)) as tif:
            page = tif.pages[0]
            img_h, img_w = page.shape[:2]

    tiles = []
    skipped = 0

    for y0 in range(0, img_h, TILE_SIZE):
        for x0 in range(0, img_w, TILE_SIZE):
            y1 = min(y0 + TILE_SIZE, img_h)
            x1 = min(x0 + TILE_SIZE, img_w)
            if y1 - y0 < EDGE_MARGIN * 2 or x1 - x0 < EDGE_MARGIN * 2:
                continue

            tile = extract_region_vips(image_path, x0, y0, x1 - x0, y1 - y0)

            cy = (y1 - y0) // 2
            cx = (x1 - x0) // 2
            half = min(128, (y1 - y0) // 2, (x1 - x0) // 2)
            probe = tile[cy - half:cy + half, cx - half:cx + half]
            cls = classify_tile(probe)
            if cls == "red":
                skipped += 1
                continue
            if tile_filter == "green" and cls == "yellow":
                skipped += 1
                continue

            tile_name = f"tile_y{y0:06d}_x{x0:06d}.tiff"
            tifffile.imwrite(str(tile_dir / tile_name), tile)
            tiles.append({
                "tile_name": tile_name,
                "y0": y0,
                "x0": x0,
                "height_px": y1 - y0,
                "width_px": x1 - x0,
                "tile_class": cls,
            })

    print(f"  Extracted {len(tiles)} tiles ({skipped} skipped) [backend=vips]")
    return tiles, img_h, img_w


def extract_tiles(image_path, tile_dir, tile_filter="auto", img_h=None, img_w=None):
    """Extract tiles using the most reliable available backend."""
    if VIPS_CMD:
        return extract_tiles_vips(image_path, tile_dir, tile_filter=tile_filter, img_h=img_h, img_w=img_w)
    return extract_tiles_tifffile(image_path, tile_dir, tile_filter=tile_filter)


def compute_image_ibg(image_path, chunk_rows=2048, max_full_pixels=200_000_000,
                      sample_grid=8, sample_box=256):
    """Compute a single 95th-percentile background value for the source image.

    For very large images, use a regular grid of sample boxes rather than a full
    raster scan. This preserves an image-level background reference while keeping
    runtime practical for whole-slide images.
    """
    store = tifffile.imread(str(image_path), aszarr=True)
    z = zarr.open(store, mode="r")
    if isinstance(z, zarr.Group):
        z = z["0"]
    if z.ndim == 3:
        z = z[0]

    img_h, img_w = z.shape[:2]
    hist = np.zeros(256, dtype=np.int64)
    total = 0

    try:
        total_pixels = img_h * img_w
        if total_pixels <= max_full_pixels:
            for y0 in range(0, img_h, chunk_rows):
                y1 = min(y0 + chunk_rows, img_h)
                chunk = np.array(z[y0:y1, :])
                if chunk.dtype == np.uint16:
                    chunk = (chunk >> 8).astype(np.uint8)
                elif chunk.dtype != np.uint8:
                    chunk = np.clip(chunk, 0, 255).astype(np.uint8)

                hist += np.bincount(chunk.ravel(), minlength=256)
                total += chunk.size
        else:
            box_h = min(sample_box, img_h)
            box_w = min(sample_box, img_w)
            y_positions = np.linspace(0, max(0, img_h - box_h), num=sample_grid, dtype=int)
            x_positions = np.linspace(0, max(0, img_w - box_w), num=sample_grid, dtype=int)
            for y0 in y_positions:
                for x0 in x_positions:
                    chunk = np.array(z[y0:y0 + box_h, x0:x0 + box_w])
                    if chunk.dtype == np.uint16:
                        chunk = (chunk >> 8).astype(np.uint8)
                    elif chunk.dtype != np.uint8:
                        chunk = np.clip(chunk, 0, 255).astype(np.uint8)

                    hist += np.bincount(chunk.ravel(), minlength=256)
                    total += chunk.size

        if total == 0:
            raise RuntimeError(f"Could not compute I_bg for {image_path}")

        target = total * 0.95
        cumulative = np.cumsum(hist)
        i_bg = int(np.searchsorted(cumulative, target))
        return i_bg
    finally:
        store.close()


def compute_image_ibg_vips(image_path, sample_grid=8, sample_box=256):
    """Compute a 95th-percentile background from vips-cropped sample boxes."""
    with tifffile.TiffFile(str(image_path)) as tif:
        page = tif.pages[0]
        img_h, img_w = page.shape[:2]

    box_h = min(sample_box, img_h)
    box_w = min(sample_box, img_w)
    y_positions = np.linspace(0, max(0, img_h - box_h), num=sample_grid, dtype=int)
    x_positions = np.linspace(0, max(0, img_w - box_w), num=sample_grid, dtype=int)

    hist = np.zeros(256, dtype=np.int64)
    total = 0
    for y0 in y_positions:
        for x0 in x_positions:
            chunk = extract_region_vips(image_path, x0, y0, box_w, box_h)
            hist += np.bincount(chunk.ravel(), minlength=256)
            total += chunk.size

    if total == 0:
        raise RuntimeError(f"Could not compute I_bg for {image_path}")

    target = total * 0.95
    cumulative = np.cumsum(hist)
    return int(np.searchsorted(cumulative, target))


def load_background_cache(csv_path):
    """Load per-image background estimates from an existing measurement CSV."""
    if csv_path is None:
        return {}
    csv_path = Path(csv_path)
    if not csv_path.exists():
        return {}

    vals = defaultdict(list)
    with open(csv_path) as f:
        reader = csv.DictReader(f)
        for row in reader:
            fname = row.get("filename")
            i_bg = row.get("i_bg")
            if not fname or i_bg in (None, ""):
                continue
            try:
                vals[fname].append(float(i_bg))
            except ValueError:
                continue

    out = {}
    for fname, items in vals.items():
        if items:
            out[fname] = int(round(float(np.median(items))))
    return out


def run_imagej_macro(input_dir, output_file, i_bg_override=None, artifact_dir=None,
                     manual_threshold=None):
    """Run the ImageJ macro via direct Java invocation.

    Bypasses the Arch `imagej` wrapper which has issues with argument passing
    and xvfb interaction under subprocess.
    """
    bg_arg = "" if i_bg_override is None else str(i_bg_override)
    artifact_arg = "" if artifact_dir is None else str(artifact_dir)
    thresh_arg = "" if manual_threshold is None else str(manual_threshold)
    args_str = (
        f"{input_dir}|{output_file}|{MIN_AREA_PX}|{MAX_AREA_PX}|"
        f"{MIN_CIRCULARITY}|{BLUR_SIGMA}|{bg_arg}|{artifact_arg}|{thresh_arg}"
    )

    ij_jar = Path("/usr/share/imagej/ij.jar")
    if not ij_jar.exists():
        raise FileNotFoundError(f"ImageJ jar not found at {ij_jar}")

    cmd = [
        "java", f"-Xmx{IMAGEJ_MEMORY_MB}m",
        "-cp", f"{ij_jar}:/usr/share/imagej/lib/*",
        "ij.ImageJ",
        "-ijpath", str(Path.home() / ".imagej"),
        "-batch", str(MACRO_PATH),
        args_str,
    ]

    print(f"  Running: java -batch {MACRO_PATH.name} ...")
    result = subprocess.run(cmd, capture_output=True, text=True, timeout=3600)

    for line in result.stdout.strip().split("\n"):
        if line.strip():
            print(f"  [IJ] {line}")

    if result.returncode != 0:
        print(f"  ImageJ stderr: {result.stderr[:500]}")
        raise RuntimeError(f"ImageJ exited with code {result.returncode}")

    return True


def write_tile_manifest(tiles, manifest_path, image_path, image_height_px, image_width_px):
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    with open(manifest_path, "w", newline="") as f:
        fieldnames = [
            "source_image", "image_height_px", "image_width_px",
            "tile_name", "y0", "x0", "height_px", "width_px", "tile_class",
        ]
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        for tile in tiles:
            w.writerow({
                "source_image": image_path.name,
                "image_height_px": image_height_px,
                "image_width_px": image_width_px,
                **tile,
            })


def read_imagej_csv(csv_path):
    """Read the CSV produced by the ImageJ macro."""
    if not csv_path.exists() or csv_path.stat().st_size == 0:
        return []
    with open(csv_path) as f:
        return list(csv.DictReader(f))


def _tile_artifact_paths(image_artifact_dir: Path | None, tile_name: str) -> dict[str, str]:
    if image_artifact_dir is None:
        return {
            "mask_path": "",
            "roi_zip_path": "",
            "tile_manifest_path": "",
            "raw_imagej_results_path": "",
        }

    tile_base = Path(tile_name).stem
    macro_dir = image_artifact_dir / "macro_artifacts"
    mask_path = macro_dir / f"{tile_base}__mask.tif"
    roi_zip_path = macro_dir / f"{tile_base}__rois.zip"
    tile_manifest_path = image_artifact_dir / "tile_manifest.csv"
    raw_results_path = image_artifact_dir / "raw_imagej_results.csv"
    return {
        "mask_path": str(mask_path) if mask_path.exists() else "",
        "roi_zip_path": str(roi_zip_path) if roi_zip_path.exists() else "",
        "tile_manifest_path": str(tile_manifest_path) if tile_manifest_path.exists() else "",
        "raw_imagej_results_path": str(raw_results_path) if raw_results_path.exists() else "",
    }


def remap_tile_coords(rows, tiles_info):
    """Convert tile-local coordinates to global image coordinates.

    tiles_info: dict mapping tile_filename -> (y0, x0)
    """
    remapped = []
    for row in rows:
        tile_fname = row["filename"]
        if tile_fname not in tiles_info:
            continue

        tile_meta = tiles_info[tile_fname]
        y0 = tile_meta["y0"]
        x0 = tile_meta["x0"]
        local_cx = float(row["centroid_x"])
        local_cy = float(row["centroid_y"])

        # Skip objects near tile edges
        tile_h = tile_meta["height_px"]
        tile_w = tile_meta["width_px"]
        if (local_cy < EDGE_MARGIN or local_cy > tile_h - EDGE_MARGIN or
                local_cx < EDGE_MARGIN or local_cx > tile_w - EDGE_MARGIN):
            continue

        row["centroid_x"] = str(round(local_cx + x0, 2))
        row["centroid_y"] = str(round(local_cy + y0, 2))
        row["tile_name"] = tile_fname
        row["tile_y0"] = str(y0)
        row["tile_x0"] = str(x0)
        row["tile_height_px"] = str(tile_h)
        row["tile_width_px"] = str(tile_w)
        remapped.append(row)

    return remapped


def process_image(image_path, tile_filter="auto", cached_i_bg=None, artifact_dir=None,
                  manual_threshold=None):
    """Process one image: tile if needed, run ImageJ, collect results."""
    image_path = Path(image_path)
    name = image_path.stem
    image_artifact_dir = Path(artifact_dir) if artifact_dir else None
    macro_artifact_dir = None
    if image_artifact_dir is not None:
        image_artifact_dir.mkdir(parents=True, exist_ok=True)
        macro_artifact_dir = image_artifact_dir / "macro_artifacts"
        macro_artifact_dir.mkdir(parents=True, exist_ok=True)

    with tifffile.TiffFile(str(image_path)) as tif:
        page = tif.pages[0]
        img_h, img_w = page.shape[:2]

    needs_tiling = img_h > TILE_SIZE * 1.5 or img_w > TILE_SIZE * 1.5
    if cached_i_bg is not None:
        global_i_bg = cached_i_bg
    else:
        global_i_bg = compute_image_ibg_vips(image_path) if VIPS_CMD else compute_image_ibg(image_path)

    with tempfile.TemporaryDirectory(prefix="ij_nuclei_") as tmpdir:
        tmpdir = Path(tmpdir)
        output_csv = tmpdir / "results.csv"

        if not needs_tiling:
            # Small image: copy directly, run ImageJ on it
            img = tifffile.imread(str(image_path))
            if img.ndim == 3:
                img = img[0]
            if img.dtype == np.uint16:
                img = (img >> 8).astype(np.uint8)

            if image_artifact_dir is not None:
                input_dir = image_artifact_dir / "input_copy"
                if input_dir.exists():
                    shutil.rmtree(input_dir)
            else:
                input_dir = tmpdir / "input"
            input_dir.mkdir(parents=True, exist_ok=True)
            tifffile.imwrite(str(input_dir / f"{name}.tiff"), img)
            print(f"  {name}: {img_w}x{img_h} (direct, I_bg={global_i_bg})")
            run_imagej_macro(
                str(input_dir),
                str(output_csv),
                i_bg_override=global_i_bg,
                artifact_dir=macro_artifact_dir,
                manual_threshold=manual_threshold,
            )

            rows = read_imagej_csv(output_csv)
            if image_artifact_dir is not None and output_csv.exists():
                shutil.copy2(output_csv, image_artifact_dir / "raw_imagej_results.csv")
            # Fix filename to original
            for r in rows:
                r["filename"] = image_path.name
                r["tile_name"] = image_path.name
                r["tile_y0"] = "0"
                r["tile_x0"] = "0"
                r["tile_height_px"] = str(img_h)
                r["tile_width_px"] = str(img_w)
                r.update(_tile_artifact_paths(image_artifact_dir, image_path.name))
            return rows

        # Large image: tile, run ImageJ, remap coordinates
        print(f"  {name}: {img_w}x{img_h} (tiled, shared I_bg={global_i_bg})")
        if image_artifact_dir is not None:
            tile_dir = image_artifact_dir / "tiles_input"
            if tile_dir.exists():
                shutil.rmtree(tile_dir)
        else:
            tile_dir = tmpdir / "tiles"
        tile_dir.mkdir(parents=True, exist_ok=True)

        tiles, _, _ = extract_tiles(image_path, tile_dir, tile_filter, img_h=img_h, img_w=img_w)
        if not tiles:
            print(f"  {name}: no content tiles found")
            return []
        if image_artifact_dir is not None:
            write_tile_manifest(tiles, image_artifact_dir / "tile_manifest.csv", image_path, img_h, img_w)

        run_imagej_macro(
            str(tile_dir),
            str(output_csv),
            i_bg_override=global_i_bg,
            artifact_dir=macro_artifact_dir,
            manual_threshold=manual_threshold,
        )

        rows = read_imagej_csv(output_csv)
        if image_artifact_dir is not None and output_csv.exists():
            shutil.copy2(output_csv, image_artifact_dir / "raw_imagej_results.csv")
        tiles_info = {t["tile_name"]: t for t in tiles}
        rows = remap_tile_coords(rows, tiles_info)

        # Set filename to original image
        for r in rows:
            r["filename"] = image_path.name
            r.update(_tile_artifact_paths(image_artifact_dir, r["tile_name"]))

        print(f"  {name}: {len(rows)} nuclei after edge filtering")
        return rows


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
    global MIN_AREA_PX, MAX_AREA_PX, MIN_CIRCULARITY, BLUR_SIGMA

    parser = argparse.ArgumentParser(description="ImageJ nucleus IOD measurement (Hardie et al. 2002)")
    parser.add_argument("image", nargs="?", help="Single image path")
    parser.add_argument("--batch", action="store_true", help="Process from master CSV")
    parser.add_argument("--resume", action="store_true", help="Skip images already in output CSV")
    parser.add_argument("--image-type", choices=["brightfield", "pmount", "both"],
                        default="both", help="Which image types to process")
    parser.add_argument("--limit", type=int, default=0)
    parser.add_argument("--tile-filter", choices=["green", "auto", "all"], default="auto")
    parser.add_argument("--output-dir", "-o", type=Path,
                        default=PROJECT / "output" / "nucleus_iod")
    parser.add_argument("--artifact-dir", type=Path,
                        help="Persistent per-image artifact root for copied inputs, masks, ROIs, and tile manifests")
    parser.add_argument("--background-cache", type=Path,
                        default=DEFAULT_BACKGROUND_CACHE,
                        help="Existing measurement CSV used to seed per-image I_bg values")
    parser.add_argument("--min-size", type=int, default=MIN_AREA_PX)
    parser.add_argument("--max-size", type=int, default=MAX_AREA_PX)
    parser.add_argument("--min-circ", type=float, default=MIN_CIRCULARITY)
    parser.add_argument("--blur-sigma", type=float, default=BLUR_SIGMA)
    args = parser.parse_args()

    MIN_AREA_PX = args.min_size
    MAX_AREA_PX = args.max_size
    MIN_CIRCULARITY = args.min_circ
    BLUR_SIGMA = args.blur_sigma

    # Check ImageJ is available
    if not shutil.which(IMAGEJ_CMD):
        print(f"ERROR: '{IMAGEJ_CMD}' not found in PATH")
        sys.exit(1)

    args.output_dir.mkdir(parents=True, exist_ok=True)
    bg_cache = load_background_cache(args.background_cache)
    if args.background_cache:
        print(f"Loaded cached I_bg for {len(bg_cache)} images from {args.background_cache}")

    if args.image:
        img_name = Path(args.image).name
        image_artifact_dir = None
        if args.artifact_dir:
            image_artifact_dir = args.artifact_dir / Path(args.image).stem
        rows = process_image(args.image, tile_filter=args.tile_filter,
                             cached_i_bg=bg_cache.get(img_name),
                             artifact_dir=image_artifact_dir)
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

    csv_path = args.output_dir / "nucleus_iod_measurements.csv"
    progress_path = args.output_dir / "imagej_progress.json"

    # Resume: find already-completed filenames from output CSV
    completed = set()
    if args.resume and csv_path.exists():
        with open(csv_path) as f:
            for row in csv.DictReader(f):
                completed.add(row.get("filename", ""))
        print(f"Resuming: {len(completed)} images already done")

    remaining = [j for j in jobs if j["filename"] not in completed]
    print(f"Batch: {len(remaining)}/{len(jobs)} images to process ({', '.join(types)})")

    t_start = time.time()
    total_nuclei = 0
    processed = 0

    for i, job in enumerate(remaining):
        print(f"\n[{len(completed) + processed + 1}/{len(jobs)}] {job['filename']} ({job.get('species', '?')}, {job['image_type']})")
        try:
            rows = process_image(job["path"], tile_filter=args.tile_filter,
                                 cached_i_bg=bg_cache.get(job["filename"]),
                                 artifact_dir=(args.artifact_dir / Path(job["filename"]).stem) if args.artifact_dir else None)
            meta = {
                "slide_id": job.get("slide_id", ""),
                "specimen_id": job.get("specimen_id", ""),
                "species": job.get("species", ""),
                "image_type": job["image_type"],
            }
            save_csv(rows, csv_path, metadata=meta)
            total_nuclei += len(rows)
            processed += 1
            completed.add(job["filename"])
            print(f"  {len(rows)} nuclei (total={total_nuclei})")

            # Checkpoint progress
            if processed % 3 == 0:
                _save_progress(progress_path, completed, total_nuclei, t_start, len(jobs))
        except Exception as e:
            print(f"  ERROR: {e}")
            import traceback; traceback.print_exc()

    _save_progress(progress_path, completed, total_nuclei, t_start, len(jobs))
    elapsed = time.time() - t_start
    print(f"\nDone: {total_nuclei} nuclei from {processed} images [{elapsed:.0f}s]")
    if csv_path.exists():
        print(f"Results: {csv_path}")


def _save_progress(path, completed, total_nuclei, t_start, total_jobs):
    elapsed = time.time() - t_start
    data = {
        "completed": sorted(completed),
        "count": len(completed),
        "total_jobs": total_jobs,
        "total_nuclei": total_nuclei,
        "elapsed_s": round(elapsed),
        "updated": time.strftime("%Y-%m-%d %H:%M:%S"),
    }
    tmp = path.with_suffix(".tmp")
    with open(tmp, "w") as f:
        json.dump(data, f, indent=2)
    tmp.rename(path)


if __name__ == "__main__":
    main()
