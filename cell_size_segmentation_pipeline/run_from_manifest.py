#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
import sys
import time
from pathlib import Path

PROJECT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT / "src"))
sys.path.insert(0, str(PROJECT / "Cellsize_segmentation_cellpose_pipeline" / "scripts"))
sys.path.insert(0, str(PROJECT / "scripts"))

from cellprofiler_tools.pipeline_runs import (  # noqa: E402
    current_git_hash,
    load_manifest_rows,
    resolve_image_path,
)
from cellprofiler_tools.convergence import interleave_by_species  # noqa: E402


def _load_backend(name: str):
    if name == "cellpose":
        from segment_cells import get_model, process_image  # noqa: E402

        return {
            "initialize": lambda *, gpu, diameter, settings_json, threshold_map, cellpose_model: get_model(
                gpu=gpu,
                diameter=diameter,
                pretrained_model=cellpose_model,
            ),
            "process_image": process_image,
            "supports_threshold_map": False,
            "supports_settings_json": False,
        }

    if name == "rules":
        from rule_based_cell_size import load_settings, load_threshold_map, process_image  # noqa: E402

        def initialize(*, gpu, diameter, settings_json, threshold_map, cellpose_model):
            _ = gpu
            _ = diameter
            _ = cellpose_model
            settings = load_settings(settings_json)
            thresholds = load_threshold_map(threshold_map)
            return settings, thresholds

        return {
            "initialize": initialize,
            "process_image": process_image,
            "supports_threshold_map": True,
            "supports_settings_json": True,
        }

    raise ValueError(f"Unknown cell-size backend: {name}")

CSV_COLUMNS = [
    "label", "seed_label", "area_px", "area_um2", "solidity", "circularity",
    "iod", "mean_od", "centroid_y", "centroid_x", "i_bg",
    "tile_name", "tile_y0", "tile_x0", "tile_h", "tile_w", "tile_score",
    "mask_path", "raw_mask_path", "mask_label_id", "tile_manifest_path", "overlay_path",
    "nucleus_area_px_seed", "nucleus_area_um2_seed", "nc_ratio_seed",
    "distance_over_cell_radius_seed", "cell_extent", "cell_edge_touch",
    "cell_threshold", "nucleus_threshold", "source_image_path",
    "run_manifest_path",
    "filename", "slide_id", "specimen_id", "species", "image_type",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--gpu", action="store_true")
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--debug", action="store_true")
    parser.add_argument("--no-crops", action="store_true")
    parser.add_argument("--tile-filter", choices=["green", "auto", "yellow", "all"], default="auto")
    parser.add_argument("--max-cells-image", type=int, default=500)
    parser.add_argument(
        "--max-cells-species",
        type=int,
        default=0,
        help="Stop scheduling more images for a species once this many cells have been collected; 0 disables the cap.",
    )
    parser.add_argument("--max-tiles-image", type=int, default=0)
    parser.add_argument("--diameter", type=float, default=None)
    parser.add_argument("--backend", choices=["rules", "cellpose"], default="rules")
    parser.add_argument(
        "--cellpose-model",
        type=Path,
        help="Custom Cellpose model path for the cellpose backend; defaults to cpsam if omitted.",
    )
    parser.add_argument("--settings-json", type=Path, help="Backend settings JSON (rules backend)")
    parser.add_argument("--threshold-csv", type=Path, help="Per-image threshold CSV (rules backend)")
    parser.add_argument("--image-type", choices=["brightfield", "pmount", "all"],
                        default="brightfield",
                        help="Process only this image type (default: brightfield)")
    return parser.parse_args()


def append_rows(csv_path: Path, rows: list[dict]) -> None:
    if not rows:
        return
    write_header = not csv_path.exists() or csv_path.stat().st_size == 0
    with open(csv_path, "a", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=CSV_COLUMNS, extrasaction="ignore")
        if write_header:
            writer.writeheader()
        writer.writerows(rows)


def write_index(index_path: Path, index_rows: list[dict]) -> None:
    with open(index_path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=index_rows[0].keys() if index_rows else [
            "filename", "image_type", "species", "slide_id", "specimen_id",
            "status", "source_image_path", "measurement_csv_path", "tile_manifest_path",
            "mask_root", "overlay_path",
        ])
        writer.writeheader()
        writer.writerows(index_rows)


def build_run_manifest(args: argparse.Namespace, rows: list[dict]) -> dict:
    return {
        "timestamp": time.strftime("%Y-%m-%dT%H:%M:%S"),
        "git_hash": current_git_hash(),
        "manifest_path": str(args.manifest.resolve()),
        "output_dir": str(args.output_dir),
        "image_count": len(rows),
        "backend": args.backend,
        "gpu": args.gpu,
        "tile_filter": args.tile_filter,
        "max_cells_image": args.max_cells_image,
        "max_cells_species": args.max_cells_species,
        "max_tiles_image": args.max_tiles_image,
        "diameter_override": args.diameter,
        "cellpose_model": str(args.cellpose_model.resolve()) if args.cellpose_model else "",
        "settings_json": str(args.settings_json.resolve()) if args.settings_json else "",
        "threshold_csv": str(args.threshold_csv.resolve()) if args.threshold_csv else "",
        "crops_enabled": not args.no_crops,
        "traceability_note": "Current core writes measurement CSVs, per-image/tile label masks, and tile manifests. Overlay generation remains optional.",
        "dry_run": args.dry_run,
    }


def main() -> None:
    args = parse_args()
    backend = _load_backend(args.backend)
    rows = load_manifest_rows(args.manifest)
    if args.backend == "rules":
        before = len(rows)
        rows = [r for r in rows if r["image_type"] == "brightfield"]
        if before != len(rows):
            print(f"Cell-size rules backend: restricted to {len(rows)} brightfield images")
    if args.image_type != "all":
        rows = [r for r in rows if r["image_type"] == args.image_type]
        print(f"Cell-size: filtered to {len(rows)} {args.image_type} images")
    if args.max_cells_species > 0:
        rows = interleave_by_species(rows)
        print(f"Cell-size: interleaved images by species with cap={args.max_cells_species}")
    args.output_dir.mkdir(parents=True, exist_ok=True)
    combined_csv = args.output_dir / "all_measurements.csv"
    run_manifest_path = args.output_dir / "run_manifest.json"
    index_path = args.output_dir / "image_index.csv"
    index_rows = []
    species_counts: dict[str, int] = {}
    run_manifest_path.write_text(json.dumps(build_run_manifest(args, rows), indent=2))

    if args.resume and args.max_cells_species > 0 and combined_csv.exists():
        with combined_csv.open(newline="") as handle:
            for row in csv.DictReader(handle):
                species = str(row.get("species", "")).strip()
                if species:
                    species_counts[species] = species_counts.get(species, 0) + 1

    state = None
    if not args.dry_run:
        state = backend["initialize"](
            gpu=args.gpu,
            diameter=args.diameter,
            settings_json=args.settings_json,
            threshold_map=args.threshold_csv,
            cellpose_model=str(args.cellpose_model.resolve()) if args.cellpose_model else None,
        )
    settings = thresholds = None
    if isinstance(state, tuple) and len(state) == 2:
        settings, thresholds = state

    for row in rows:
        image_path = resolve_image_path(row)
        image_output_dir = args.output_dir / row["image_type"]
        stem = Path(row["filename"]).stem
        measurement_csv = image_output_dir / f"{stem}_measurements.csv"
        overlay_path = image_output_dir / "debug" / f"{stem}_overlay.png"
        tile_manifest_path = image_output_dir / "tile_manifests" / f"{stem}_tile_manifest.csv"
        mask_root = image_output_dir / "masks" / stem
        status = "planned"
        species = str(row.get("species", "")).strip()

        if args.max_cells_species > 0 and species and species_counts.get(species, 0) >= args.max_cells_species:
            status = "skipped_species_cap"
            index_rows.append(
                {
                    "filename": row["filename"],
                    "image_type": row["image_type"],
                    "species": row["species"],
                    "slide_id": row["slide_id"],
                    "specimen_id": row["specimen_id"],
                    "status": status,
                    "source_image_path": str(image_path),
                    "measurement_csv_path": str(measurement_csv),
                    "tile_manifest_path": str(tile_manifest_path),
                    "mask_root": str(mask_root),
                    "overlay_path": str(overlay_path),
                }
            )
            write_index(index_path, index_rows)
            continue

        if args.resume and measurement_csv.exists():
            status = "skipped_existing"
            index_rows.append(
                {
                    "filename": row["filename"],
                    "image_type": row["image_type"],
                    "species": row["species"],
                    "slide_id": row["slide_id"],
                    "specimen_id": row["specimen_id"],
                    "status": status,
                    "source_image_path": str(image_path),
                    "measurement_csv_path": str(measurement_csv),
                    "tile_manifest_path": str(tile_manifest_path),
                    "mask_root": str(mask_root),
                    "overlay_path": str(overlay_path),
                }
            )
            write_index(index_path, index_rows)
            continue

        if not args.dry_run:
            try:
                kwargs = {
                    "debug": args.debug,
                    "crops": not args.no_crops,
                    "tile_filter": args.tile_filter,
                    "image_type": row["image_type"],
                    "max_cells": args.max_cells_image,
                }
                if args.backend == "cellpose":
                    kwargs["diameter_override"] = args.diameter
                    kwargs["pretrained_model"] = str(args.cellpose_model.resolve()) if args.cellpose_model else None
                else:
                    kwargs["max_tiles"] = args.max_tiles_image
                    kwargs["settings"] = settings
                    kwargs["manual_threshold"] = thresholds.get(row["filename"]) if thresholds else None

                measurements = backend["process_image"](
                    image_path,
                    image_output_dir,
                    **kwargs,
                )
                for measurement in measurements:
                    measurement["filename"] = row["filename"]
                    measurement["slide_id"] = row["slide_id"]
                    measurement["specimen_id"] = row["specimen_id"]
                    measurement["species"] = row["species"]
                    measurement["image_type"] = row["image_type"]
                    measurement["source_image_path"] = str(image_path.resolve())
                    measurement["run_manifest_path"] = str(run_manifest_path.resolve())
                append_rows(combined_csv, measurements)
                if args.max_cells_species > 0 and species:
                    species_counts[species] = species_counts.get(species, 0) + len(measurements)
                status = "completed"
            except Exception as e:
                print(f"  ERROR {row['filename']}: {e}")
                status = "error"

        index_rows.append(
            {
                "filename": row["filename"],
                "image_type": row["image_type"],
                "species": row["species"],
                "slide_id": row["slide_id"],
                "specimen_id": row["specimen_id"],
                "status": status,
                "source_image_path": str(image_path),
                "measurement_csv_path": str(measurement_csv),
                "tile_manifest_path": str(tile_manifest_path),
                "mask_root": str(mask_root),
                "overlay_path": str(overlay_path),
            }
        )
        write_index(index_path, index_rows)
    print(f"Cell-size image index: {index_path}")
    print(f"Cell-size manifest: {args.output_dir / 'run_manifest.json'}")


if __name__ == "__main__":
    main()
