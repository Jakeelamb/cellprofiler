#!/usr/bin/env python3
"""Prepare tile images and backfilled cell measurements for mixed Cellpose+YOLO linkage runs."""

from __future__ import annotations

import argparse
import csv
import json
import subprocess
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
import tifffile
import zarr
from PIL import Image

PROJECT = Path(__file__).resolve().parent.parent

import sys

sys.path.insert(0, str(PROJECT / "src"))

from cellprofiler_tools.pipeline_runs import resolve_image_path  # noqa: E402
from imagej_nucleus_iod import VIPS_CMD  # noqa: E402


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--cell-csv",
        type=Path,
        default=PROJECT / "output" / "runs" / "full_dataset_v1" / "cell_size_segmentation" / "all_measurements.csv",
    )
    parser.add_argument(
        "--cell-run-manifest",
        type=Path,
        default=PROJECT / "output" / "runs" / "full_dataset_v1" / "cell_size_segmentation" / "run_manifest.json",
    )
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--image-type", default="brightfield", choices=["brightfield"])
    parser.add_argument("--tile-format", choices=["png", "tiff"], default="png")
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


def normalize_uint8(tile: np.ndarray) -> np.ndarray:
    arr = np.asarray(tile)
    if arr.ndim == 3 and arr.shape[0] in (1, 3, 4):
        arr = arr[0]
    if arr.ndim == 3 and arr.shape[-1] in (3, 4):
        arr = arr[..., 0]
    arr = np.squeeze(arr)
    if arr.dtype == np.uint16:
        arr = (arr >> 8).astype(np.uint8)
    elif arr.dtype != np.uint8:
        arr = np.clip(arr, 0, 255).astype(np.uint8)
    return arr


def extract_region_vips_fast(source_path: Path, x0: int, y0: int, w: int, h: int) -> np.ndarray:
    if not VIPS_CMD:
        raise RuntimeError("libvips is not available")
    with tempfile.NamedTemporaryFile(prefix="vips_crop_", suffix=".tif", delete=False) as tmp:
        tmp_path = Path(tmp.name)
    try:
        cmd = [
            VIPS_CMD,
            "crop",
            str(source_path),
            str(tmp_path),
            str(int(x0)),
            str(int(y0)),
            str(int(w)),
            str(int(h)),
        ]
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
        if result.returncode != 0:
            raise RuntimeError(f"vips crop failed for {source_path}: {result.stderr[:400]}")
        return normalize_uint8(np.asarray(Image.open(tmp_path)))
    finally:
        tmp_path.unlink(missing_ok=True)


def load_cells(path: Path, image_type: str) -> pd.DataFrame:
    df = pd.read_csv(path, keep_default_na=False)
    if "image_type" in df.columns:
        df = df[df["image_type"].astype(str) == image_type].copy()
    if df.empty:
        raise SystemExit(f"no cell rows found for image_type={image_type}")
    return df.reset_index(drop=True)


def backfill_cell_traceability(cells: pd.DataFrame, run_manifest_path: Path) -> pd.DataFrame:
    out = cells.copy()
    if "source_image_path" not in out.columns:
        out["source_image_path"] = ""
    if "run_manifest_path" not in out.columns:
        out["run_manifest_path"] = ""

    empty_source = out["source_image_path"].astype(str).str.strip().eq("")
    if empty_source.any():
        unique_pairs = (
            out.loc[empty_source, ["filename", "image_type"]]
            .drop_duplicates()
            .itertuples(index=False)
        )
        lookup = {
            (str(filename), str(image_type)): str(
                resolve_image_path({"filename": str(filename), "image_type": str(image_type)}).resolve()
            )
            for filename, image_type in unique_pairs
        }
        out.loc[empty_source, "source_image_path"] = [
            lookup[(str(filename), str(image_type))]
            for filename, image_type in out.loc[empty_source, ["filename", "image_type"]].itertuples(index=False)
        ]

    empty_manifest = out["run_manifest_path"].astype(str).str.strip().eq("")
    if empty_manifest.any():
        out.loc[empty_manifest, "run_manifest_path"] = str(run_manifest_path.resolve())

    return out


def unique_tiles(cells: pd.DataFrame) -> pd.DataFrame:
    cols = [
        "filename",
        "slide_id",
        "specimen_id",
        "species",
        "image_type",
        "source_image_path",
        "tile_name",
        "tile_y0",
        "tile_x0",
        "tile_h",
        "tile_w",
        "tile_score",
        "tile_manifest_path",
    ]
    keep_cols = [col for col in cols if col in cells.columns]
    tiles = cells[keep_cols].drop_duplicates().copy()
    numeric_cols = ["tile_y0", "tile_x0", "tile_h", "tile_w", "tile_score"]
    for col in numeric_cols:
        if col in tiles.columns:
            tiles[col] = pd.to_numeric(tiles[col], errors="coerce")
    tiles = tiles.sort_values(["filename", "tile_y0", "tile_x0"]).reset_index(drop=True)
    return tiles


def _load_image_array(source_path: Path):
    try:
        store = tifffile.imread(str(source_path), aszarr=True)
        arr = zarr.open(store, mode="r")
        if getattr(arr, "ndim", 0) == 3:
            arr = arr[0]
        return arr, store
    except Exception:
        with tifffile.TiffFile(str(source_path)) as tif:
            page = tif.pages[0]
            arr = page.asarray(out="memmap")
        return arr, None


def extract_tiles(tiles: pd.DataFrame, tiles_dir: Path, tile_format: str) -> list[dict[str, object]]:
    manifest_rows: list[dict[str, object]] = []
    grouped = list(tiles.groupby("source_image_path", sort=True))
    for idx, (source_path_str, group) in enumerate(grouped, start=1):
        source_path = Path(str(source_path_str))
        require_exists(source_path)
        print(
            f"[{idx}/{len(grouped)}] extracting {len(group)} tiles from {source_path.name}",
            flush=True,
        )
        image_array = None
        store = None
        if not VIPS_CMD:
            image_array, store = _load_image_array(source_path)
        try:
            for row in group.itertuples(index=False):
                y0 = int(row.tile_y0)
                x0 = int(row.tile_x0)
                h = int(row.tile_h)
                w = int(row.tile_w)
                if VIPS_CMD:
                    tile = extract_region_vips_fast(source_path, x0, y0, w, h)
                else:
                    tile = np.array(image_array[y0 : y0 + h, x0 : x0 + w])
                    tile = normalize_uint8(tile)

                image_stem = Path(str(row.filename)).stem
                tile_stem = Path(str(row.tile_name)).stem
                tile_path = tiles_dir / image_stem / f"{tile_stem}.{tile_format}"
                tile_path.parent.mkdir(parents=True, exist_ok=True)
                if tile_format == "png":
                    Image.fromarray(tile).save(tile_path)
                else:
                    tifffile.imwrite(str(tile_path), tile, compression="zlib")

                manifest_rows.append(
                    {
                        "filename": str(row.filename),
                        "slide_id": str(getattr(row, "slide_id", "")),
                        "specimen_id": str(getattr(row, "specimen_id", "")),
                        "species": str(getattr(row, "species", "")),
                        "image_type": str(getattr(row, "image_type", "")),
                        "image_path": str(tile_path.resolve()),
                        "source_image_path": str(source_path.resolve()),
                        "tile_name": str(row.tile_name),
                        "tile_y0": y0,
                        "tile_x0": x0,
                        "tile_h": h,
                        "tile_w": w,
                        "tile_score": float(getattr(row, "tile_score", 1.0) or 1.0),
                        "cell_tile_manifest_path": str(getattr(row, "tile_manifest_path", "")),
                    }
                )
        finally:
            if store is not None and hasattr(store, "close"):
                store.close()
    return manifest_rows


def write_tile_manifest(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    args = parse_args()
    require_exists(args.cell_csv)
    require_exists(args.cell_run_manifest)

    args.output_dir.mkdir(parents=True, exist_ok=True)
    clear_dir(args.output_dir)
    tiles_dir = args.output_dir / "tiles"
    tiles_dir.mkdir(parents=True, exist_ok=True)

    cells = load_cells(args.cell_csv, args.image_type)
    cells = backfill_cell_traceability(cells, args.cell_run_manifest)
    tiles = unique_tiles(cells)
    tile_rows = extract_tiles(tiles, tiles_dir, args.tile_format)
    if not tile_rows:
        raise SystemExit("no tile rows extracted")

    backfilled_csv = args.output_dir / "cell_measurements_backfilled.csv"
    cells.to_csv(backfilled_csv, index=False)

    tile_manifest = args.output_dir / "tile_manifest.csv"
    write_tile_manifest(tile_manifest, tile_rows)

    summary = {
        "cell_csv": str(args.cell_csv.resolve()),
        "cell_run_manifest": str(args.cell_run_manifest.resolve()),
        "output_dir": str(args.output_dir.resolve()),
        "image_type": args.image_type,
        "tile_format": args.tile_format,
        "n_cell_rows": int(len(cells)),
        "n_unique_tiles": int(len(tiles)),
        "n_unique_images": int(tiles["filename"].nunique()),
        "cell_measurements_backfilled_csv": str(backfilled_csv.resolve()),
        "tile_manifest_csv": str(tile_manifest.resolve()),
    }
    (args.output_dir / "summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
