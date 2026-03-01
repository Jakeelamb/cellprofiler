#!/usr/bin/env python3
"""Run Cellpose segmentation on brightfield crops and keep isolated cells only.

This workflow is intentionally restricted to small brightfield crop images,
not whole-slide images.

Outputs (timestamped run directory):
  - per-image folders with raw image, masks, and overlay visualization
  - object-level metrics CSV with keep/reject reason
  - image-level QC summary CSV
  - concise markdown run report
"""

from __future__ import annotations

import argparse
import json
from collections import Counter
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
import tifffile
from cellpose import models
from skimage import measure, segmentation
from skimage.color import gray2rgb
from skimage.io import imsave

PROJECT = Path(__file__).resolve().parent.parent


@dataclass
class Thresholds:
    min_cell_area: int = 700
    max_cell_area: int = 20000
    min_solidity: float = 0.85
    max_eccentricity: float = 0.92
    min_extent: float = 0.40
    min_nucleus_fraction: float = 0.03
    max_nucleus_fraction: float = 0.55


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Cellpose segmentation for brightfield crops with clump filtering."
    )
    parser.add_argument(
        "--input-glob",
        action="append",
        default=[
            "data/tuning/brightfield/*.tif*",
            "data/tuning/subset/*brightfield*.tif*",
        ],
        help="Input glob(s). Repeat for multiple sets.",
    )
    parser.add_argument(
        "--output-root",
        default="results",
        help="Root directory for timestamped run output.",
    )
    parser.add_argument(
        "--run-name",
        default=None,
        help="Optional run name suffix.",
    )
    parser.add_argument(
        "--max-images",
        type=int,
        default=6,
        help="Maximum number of crop images to process.",
    )
    parser.add_argument(
        "--cell-diameter",
        type=float,
        default=62.0,
        help="Cellpose diameter for cell-scale pass.",
    )
    parser.add_argument(
        "--nucleus-diameter",
        type=float,
        default=12.0,
        help="Cellpose diameter for nucleus-scale pass.",
    )
    parser.add_argument(
        "--cell-min-size",
        type=int,
        default=250,
        help="Cellpose min_size for cell masks.",
    )
    parser.add_argument(
        "--nucleus-min-size",
        type=int,
        default=20,
        help="Cellpose min_size for nuclei masks.",
    )
    parser.add_argument(
        "--cellprob-threshold",
        type=float,
        default=0.0,
        help="Cellpose cellprob_threshold for both passes.",
    )
    parser.add_argument(
        "--invert",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Invert brightfield image before Cellpose inference.",
    )
    return parser.parse_args()


def resolve_inputs(globs: list[str], max_images: int) -> list[Path]:
    paths: set[Path] = set()
    for pattern in globs:
        paths.update(PROJECT.glob(pattern))
    images = sorted(p for p in paths if p.is_file())
    if not images:
        raise FileNotFoundError(f"No input images found for globs: {globs}")
    return images[:max_images]


def load_grayscale(path: Path) -> np.ndarray:
    arr = tifffile.imread(path)
    if arr.ndim == 3:
        arr = arr[0]
    if arr.dtype == np.uint16:
        arr = (arr >> 8).astype(np.uint8)
    elif arr.dtype != np.uint8:
        arr = np.clip(arr, 0, 255).astype(np.uint8)
    return arr


def run_cellpose(
    model: models.CellposeModel,
    image: np.ndarray,
    diameter: float,
    min_size: int,
    cellprob_threshold: float,
    invert: bool,
) -> np.ndarray:
    masks, _, _ = model.eval(
        image,
        diameter=diameter,
        min_size=min_size,
        cellprob_threshold=cellprob_threshold,
        invert=invert,
    )
    return masks.astype(np.int32, copy=False)


def _paint_boundary(rgb: np.ndarray, mask: np.ndarray, color: tuple[int, int, int]) -> None:
    boundary = segmentation.find_boundaries(mask > 0, mode="outer")
    rgb[boundary] = color


def build_keep_masks(
    cell_mask: np.ndarray,
    nucleus_mask: np.ndarray,
    thresholds: Thresholds,
) -> tuple[np.ndarray, np.ndarray, pd.DataFrame, dict[str, float | int]]:
    cell_props = measure.regionprops(cell_mask)
    nucleus_props = measure.regionprops(nucleus_mask)

    nucleus_area_by_id = {p.label: int(p.area) for p in nucleus_props}
    object_rows: list[dict[str, float | int | str]] = []
    kept_cells = np.zeros_like(cell_mask, dtype=np.int32)
    kept_nuclei = np.zeros_like(nucleus_mask, dtype=np.int32)
    reject_counter: Counter[str] = Counter()

    next_id = 1
    for c in cell_props:
        cell_id = int(c.label)
        coords = c.coords
        nucleus_ids = nucleus_mask[coords[:, 0], coords[:, 1]]
        nucleus_ids = nucleus_ids[nucleus_ids > 0]
        unique_nucleus_ids = np.unique(nucleus_ids)
        nucleus_count = int(unique_nucleus_ids.size)
        total_nucleus_overlap_area = int(nucleus_ids.size)
        nucleus_fraction = total_nucleus_overlap_area / float(c.area)

        reason = "kept"
        if (
            c.bbox[0] == 0
            or c.bbox[1] == 0
            or c.bbox[2] >= cell_mask.shape[0]
            or c.bbox[3] >= cell_mask.shape[1]
        ):
            reason = "border_touching"
        elif c.area < thresholds.min_cell_area or c.area > thresholds.max_cell_area:
            reason = "cell_area_out_of_range"
        elif c.solidity < thresholds.min_solidity:
            reason = "cell_low_solidity"
        elif c.eccentricity > thresholds.max_eccentricity:
            reason = "cell_high_eccentricity"
        elif c.extent < thresholds.min_extent:
            reason = "cell_low_extent"
        elif nucleus_count != 1:
            reason = "nucleus_overlap_or_missing"
        elif (
            nucleus_fraction < thresholds.min_nucleus_fraction
            or nucleus_fraction > thresholds.max_nucleus_fraction
        ):
            reason = "nucleus_fraction_out_of_range"

        if reason == "kept":
            kept_cells[cell_mask == cell_id] = next_id
            kept_nuclei[np.isin(nucleus_mask, unique_nucleus_ids)] = next_id
            kept_label = next_id
            next_id += 1
        else:
            kept_label = 0
            reject_counter[reason] += 1

        object_rows.append(
            {
                "cell_label_raw": cell_id,
                "kept_label": kept_label,
                "status": "kept" if reason == "kept" else "rejected",
                "reject_reason": reason if reason != "kept" else "",
                "cell_area_px": int(c.area),
                "cell_solidity": float(c.solidity),
                "cell_eccentricity": float(c.eccentricity),
                "cell_extent": float(c.extent),
                "cell_centroid_y": float(c.centroid[0]),
                "cell_centroid_x": float(c.centroid[1]),
                "nucleus_count_in_cell": nucleus_count,
                "nucleus_overlap_area_px": total_nucleus_overlap_area,
                "nucleus_fraction_of_cell": float(nucleus_fraction),
                "single_nucleus_area_px": int(nucleus_area_by_id.get(int(unique_nucleus_ids[0]), 0))
                if nucleus_count == 1
                else 0,
            }
        )

    object_df = pd.DataFrame(object_rows)
    summary = {
        "raw_cell_count": int(cell_mask.max()),
        "raw_nucleus_count": int(nucleus_mask.max()),
        "kept_count": int(kept_cells.max()),
        "rejected_count": int(max(cell_mask.max(), 0) - max(kept_cells.max(), 0)),
        "rejected_overlap_count": int(reject_counter.get("nucleus_overlap_or_missing", 0)),
        "rejected_border_count": int(reject_counter.get("border_touching", 0)),
        "raw_cell_mask_coverage": float((cell_mask > 0).mean()),
        "kept_cell_mask_coverage": float((kept_cells > 0).mean()),
    }
    if len(object_df) > 0 and summary["kept_count"] > 0:
        kept = object_df[object_df["status"] == "kept"]
        summary["kept_area_min_px"] = int(kept["cell_area_px"].min())
        summary["kept_area_max_px"] = int(kept["cell_area_px"].max())
        summary["kept_area_median_px"] = float(kept["cell_area_px"].median())
    else:
        summary["kept_area_min_px"] = 0
        summary["kept_area_max_px"] = 0
        summary["kept_area_median_px"] = 0.0
    return kept_cells, kept_nuclei, object_df, summary


def save_outputs(
    image: np.ndarray,
    cell_mask: np.ndarray,
    nucleus_mask: np.ndarray,
    kept_cells: np.ndarray,
    kept_nuclei: np.ndarray,
    image_dir: Path,
) -> None:
    image_dir.mkdir(parents=True, exist_ok=True)
    raw_png = image_dir / "raw_crop.png"
    overlay_png = image_dir / "overlay_kept_vs_rejected.png"
    full_overlay_png = image_dir / "overlay_raw_masks.png"

    tifffile.imwrite(image_dir / "cell_mask_raw.tiff", cell_mask.astype(np.uint16))
    tifffile.imwrite(image_dir / "nucleus_mask_raw.tiff", nucleus_mask.astype(np.uint16))
    tifffile.imwrite(image_dir / "cell_mask_kept.tiff", kept_cells.astype(np.uint16))
    tifffile.imwrite(image_dir / "nucleus_mask_kept.tiff", kept_nuclei.astype(np.uint16))

    imsave(raw_png, image, check_contrast=False)

    base = gray2rgb(image)
    overlay = base.copy()
    _paint_boundary(overlay, cell_mask, (255, 70, 70))
    _paint_boundary(overlay, nucleus_mask, (70, 170, 255))
    imsave(full_overlay_png, overlay, check_contrast=False)

    keep_overlay = base.copy()
    rejected_mask = np.where((cell_mask > 0) & (kept_cells == 0), cell_mask, 0)
    _paint_boundary(keep_overlay, rejected_mask, (255, 50, 50))
    _paint_boundary(keep_overlay, kept_cells, (40, 220, 80))
    _paint_boundary(keep_overlay, kept_nuclei, (80, 220, 255))
    imsave(overlay_png, keep_overlay, check_contrast=False)


def build_run_report(
    run_dir: Path,
    image_summary: pd.DataFrame,
    object_df: pd.DataFrame,
    args: argparse.Namespace,
    thresholds: Thresholds,
) -> None:
    lines = [
        "# Cellpose Brightfield Crop Run Report",
        "",
        f"- Timestamp: {datetime.now().isoformat(timespec='seconds')}",
        f"- Images processed: {len(image_summary)}",
        "- Pretrained model: `cpsam` (Cellpose v4 bundled model)",
        f"- Invert input: `{args.invert}`",
        f"- Diameters: nuclei={args.nucleus_diameter}, cells={args.cell_diameter}",
        "",
        "## QC Totals",
    ]
    total_raw_cells = int(image_summary["raw_cell_count"].sum())
    total_kept = int(image_summary["kept_count"].sum())
    total_rejected_overlap = int(image_summary["rejected_overlap_count"].sum())
    lines.extend(
        [
            f"- Raw cell detections: {total_raw_cells}",
            f"- Kept isolated cells: {total_kept}",
            f"- Rejected for overlap/missing nucleus: {total_rejected_overlap}",
            f"- Kept fraction: {total_kept / total_raw_cells:.3f}" if total_raw_cells > 0 else "- Kept fraction: 0.000",
            "",
            "## Keep/Reject Heuristics",
            f"- cell area in [{thresholds.min_cell_area}, {thresholds.max_cell_area}] px",
            f"- cell solidity >= {thresholds.min_solidity}",
            f"- cell eccentricity <= {thresholds.max_eccentricity}",
            f"- cell extent >= {thresholds.min_extent}",
            f"- exactly one nucleus overlap per kept cell",
            (
                f"- nucleus overlap fraction in "
                f"[{thresholds.min_nucleus_fraction}, {thresholds.max_nucleus_fraction}]"
            ),
            "",
            "## Top Reject Reasons",
        ]
    )
    if len(object_df) > 0:
        rejected = object_df[object_df["status"] == "rejected"]
        if len(rejected) > 0:
            for reason, count in rejected["reject_reason"].value_counts().head(6).items():
                lines.append(f"- {reason}: {int(count)}")
        else:
            lines.append("- none")
    else:
        lines.append("- no objects detected")

    (run_dir / "run_report.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    args = parse_args()
    inputs = resolve_inputs(args.input_glob, args.max_images)
    thresholds = Thresholds()

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    run_name = f"cellpose_brightfield_crops_{timestamp}"
    if args.run_name:
        run_name = f"{run_name}_{args.run_name}"
    run_dir = PROJECT / args.output_root / run_name
    images_dir = run_dir / "images"
    run_dir.mkdir(parents=True, exist_ok=True)

    config = {
        "timestamp": timestamp,
        "inputs": [str(p) for p in inputs],
        "cellpose_model": "cpsam",
        "invert": args.invert,
        "nucleus_diameter": args.nucleus_diameter,
        "cell_diameter": args.cell_diameter,
        "nucleus_min_size": args.nucleus_min_size,
        "cell_min_size": args.cell_min_size,
        "cellprob_threshold": args.cellprob_threshold,
        "thresholds": thresholds.__dict__,
    }
    (run_dir / "config.json").write_text(json.dumps(config, indent=2), encoding="utf-8")

    print(f"Running Cellpose brightfield crop workflow on {len(inputs)} images")
    print(f"Output: {run_dir}")
    print("Model: cpsam (published Cellpose pretrained model)")

    model = models.CellposeModel(gpu=False, pretrained_model="cpsam")

    summary_rows: list[dict[str, float | int | str]] = []
    object_rows: list[pd.DataFrame] = []
    for idx, path in enumerate(inputs, start=1):
        print(f"[{idx}/{len(inputs)}] {path}")
        image = load_grayscale(path)

        nucleus_raw = run_cellpose(
            model=model,
            image=image,
            diameter=args.nucleus_diameter,
            min_size=args.nucleus_min_size,
            cellprob_threshold=args.cellprob_threshold,
            invert=args.invert,
        )
        cell_raw = run_cellpose(
            model=model,
            image=image,
            diameter=args.cell_diameter,
            min_size=args.cell_min_size,
            cellprob_threshold=args.cellprob_threshold,
            invert=args.invert,
        )
        kept_cells, kept_nuclei, object_df, summary = build_keep_masks(
            cell_mask=cell_raw,
            nucleus_mask=nucleus_raw,
            thresholds=thresholds,
        )
        object_df.insert(0, "image_name", path.stem)
        image_out_dir = images_dir / path.stem
        save_outputs(
            image=image,
            cell_mask=cell_raw,
            nucleus_mask=nucleus_raw,
            kept_cells=kept_cells,
            kept_nuclei=kept_nuclei,
            image_dir=image_out_dir,
        )
        summary_row = {"image_name": path.stem, "input_path": str(path), **summary}
        summary_rows.append(summary_row)
        object_rows.append(object_df)
        print(
            "  raw cells={raw_cell_count}, kept={kept_count}, overlap rejects={rejected_overlap_count}".format(
                **summary
            )
        )

    summary_df = pd.DataFrame(summary_rows)
    objects_df = pd.concat(object_rows, ignore_index=True) if object_rows else pd.DataFrame()
    summary_df.to_csv(run_dir / "summary_metrics.csv", index=False)
    objects_df.to_csv(run_dir / "object_metrics.csv", index=False)
    build_run_report(
        run_dir=run_dir,
        image_summary=summary_df,
        object_df=objects_df,
        args=args,
        thresholds=thresholds,
    )
    print(f"Wrote summary_metrics.csv ({len(summary_df)} rows)")
    print(f"Wrote object_metrics.csv ({len(objects_df)} rows)")
    print(f"Run report: {run_dir / 'run_report.md'}")


if __name__ == "__main__":
    main()
