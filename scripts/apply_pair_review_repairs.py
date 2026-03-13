#!/usr/bin/env python3
"""Apply approved pair-review repairs and emit patched linkage tables."""

from __future__ import annotations

import argparse
import csv
import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import tifffile
from skimage.measure import regionprops

PROJECT = Path(__file__).resolve().parent.parent
import sys

sys.path.insert(0, str(PROJECT / "scripts"))

from build_cell_nucleus_linkage_report import annotate_pair_qc  # noqa: E402
from nucleus_iod_python import PIXEL_AREA_UM2  # noqa: E402
from run_linked_pair_review import (  # noqa: E402
    DEFAULT_PAIRS_CSV,
    DEFAULT_OUTPUT_DIR as DEFAULT_REVIEW_DIR,
    as_bool,
    as_float,
    as_int,
    build_repair_candidate,
    clean_text,
    load_label_array,
    load_label_slices,
    load_raw_array,
    normalize_repair_mode,
    review_key_for_row,
)

DEFAULT_LINKAGE_DIR = DEFAULT_PAIRS_CSV.parent
DEFAULT_CELL_LINKAGE_CSV = DEFAULT_LINKAGE_DIR / "cell_linkage_summary.csv.gz"
DEFAULT_DECISIONS_CSV = DEFAULT_REVIEW_DIR / "decisions.csv"
DEFAULT_OUTPUT_DIR = DEFAULT_REVIEW_DIR / "applied_repairs"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--linked-pairs-csv", type=Path, default=DEFAULT_PAIRS_CSV)
    parser.add_argument("--cell-linkage-csv", type=Path, default=DEFAULT_CELL_LINKAGE_CSV)
    parser.add_argument("--decisions-csv", type=Path, default=DEFAULT_DECISIONS_CSV)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--max-repairs", type=int, default=0, help="Optional cap for smoke tests.")
    return parser.parse_args()


def require_exists(path: Path) -> None:
    if not path.exists():
        raise FileNotFoundError(path)


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


def expand_slice(obj_slice: tuple[slice, ...], shape: tuple[int, int], pad: int) -> tuple[slice, slice]:
    y0 = max(0, obj_slice[0].start - pad)
    y1 = min(shape[0], obj_slice[0].stop + pad)
    x0 = max(0, obj_slice[1].start - pad)
    x1 = min(shape[1], obj_slice[1].stop + pad)
    return slice(y0, y1), slice(x0, x1)


def path_output_for_mask(output_dir: Path, kind: str, original_path: str) -> Path:
    src = Path(original_path)
    return output_dir / f"{kind}_masks" / src.parent.name / src.name


def load_pairs(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, low_memory=False)
    for col in df.columns:
        if col.startswith("keep_") or col in {"has_cell_match", "physical_pair_ok", "one_to_one_cell"}:
            df[col] = df[col].map(as_bool)
    df["review_key"] = [review_key_for_row(row) for row in df.to_dict(orient="records")]
    return df


def load_cell_linkage(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, low_memory=False)


def load_decisions(path: Path) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame(columns=["review_key", "decision", "repair_mode", "note", "updated_at"])
    rows: list[dict[str, Any]] = []
    with path.open(newline="") as handle:
        for row in csv.DictReader(handle):
            key = clean_text(row.get("review_key"))
            if not key:
                continue
            rows.append(
                {
                    "review_key": key,
                    "review_decision": clean_text(row.get("decision")) or "unlabeled",
                    "review_repair_mode": normalize_repair_mode(row.get("repair_mode")) or "",
                    "review_note": clean_text(row.get("note")),
                    "review_updated_at": clean_text(row.get("updated_at")),
                }
            )
    if not rows:
        return pd.DataFrame(columns=["review_key", "review_decision", "review_repair_mode", "review_note", "review_updated_at"])
    return pd.DataFrame(rows)


@dataclass
class RepairResult:
    row_index: int
    review_key: str
    kind: str
    applied: bool
    original_mask_path: str
    repaired_mask_path: str
    original_area_px: float
    repaired_area_px: float
    original_area_um2: float
    repaired_area_um2: float
    original_centroid_x: float
    original_centroid_y: float
    repaired_centroid_x: float
    repaired_centroid_y: float
    original_solidity: float
    repaired_solidity: float
    original_circularity: float
    repaired_circularity: float
    original_iod: float
    repaired_iod: float
    original_mean_od: float
    repaired_mean_od: float


def measure_binary_metrics(
    binary: np.ndarray,
    *,
    offset_y: int,
    offset_x: int,
) -> dict[str, float]:
    labels = binary.astype(np.uint8)
    props = regionprops(labels)
    if not props:
        raise ValueError("Repair candidate is empty.")
    prop = max(props, key=lambda item: item.area)
    area_px = float(prop.area)
    perimeter = float(prop.perimeter)
    circularity = 0.0 if perimeter <= 0 else float(4.0 * math.pi * area_px / (perimeter * perimeter))
    return {
        "area_px": area_px,
        "area_um2": float(area_px * PIXEL_AREA_UM2),
        "solidity": float(prop.solidity),
        "circularity": circularity,
        "centroid_y": float(prop.centroid[0] + offset_y),
        "centroid_x": float(prop.centroid[1] + offset_x),
        "equiv_radius_px": float(math.sqrt(area_px / math.pi)) if area_px > 0 else np.nan,
    }


def measure_repaired_iod(
    binary: np.ndarray,
    *,
    raw_path: str,
    local_slice: tuple[slice, slice],
    i_bg: float,
) -> dict[str, float]:
    gray = grayscale_uint8(load_raw_array(raw_path))
    coords = np.argwhere(binary)
    y = coords[:, 0] + local_slice[0].start
    x = coords[:, 1] + local_slice[1].start
    pixel_vals = np.clip(gray[y, x].astype(np.float64), 1, None)
    od_per_pixel = np.log10(max(i_bg, 1.0) / pixel_vals)
    iod = float(np.sum(od_per_pixel))
    mean_od = float(np.mean(od_per_pixel)) if len(od_per_pixel) else np.nan
    return {"iod": iod, "mean_od": mean_od}


def apply_repair_to_tile(
    *,
    kind: str,
    original_path: str,
    tile_rows: list[tuple[int, pd.Series]],
    output_dir: Path,
) -> tuple[str, list[RepairResult]]:
    source = load_label_array(original_path)
    updated = np.array(source, copy=True)
    slices = load_label_slices(original_path)
    pad = 18 if kind == "cell" else 10
    results: list[RepairResult] = []

    for row_index, row in tile_rows:
        label = as_int(row["mask_label_id"] if kind == "cell" else row["nucleus_label"], default=0)
        if label <= 0 or label - 1 >= len(slices) or slices[label - 1] is None:
            continue
        local_slice = expand_slice(slices[label - 1], updated.shape, pad)
        original_binary = source[local_slice] == label
        if not original_binary.any():
            continue
        repaired_binary = build_repair_candidate(original_binary, kind, row.to_dict(), force=True)
        occupancy = updated[local_slice]
        repaired_binary = repaired_binary & ((occupancy == 0) | (occupancy == label))
        if not repaired_binary.any():
            continue

        local_updated = np.array(updated[local_slice], copy=True)
        local_updated[local_updated == label] = 0
        local_updated[repaired_binary] = label
        updated[local_slice] = local_updated

        original_metrics = {
            "area_px": as_float(row["cell_area_px"] if kind == "cell" else row["nuc_area_px"]),
            "area_um2": as_float(row["cell_area_um2"] if kind == "cell" else row["nuc_area_um2"]),
            "centroid_x": as_float(row["cell_centroid_x"] if kind == "cell" else row["nuc_centroid_x"]),
            "centroid_y": as_float(row["cell_centroid_y"] if kind == "cell" else row["nuc_centroid_y"]),
            "solidity": as_float(row["cell_solidity"] if kind == "cell" else np.nan),
            "circularity": as_float(row["cell_circularity"] if kind == "cell" else np.nan),
            "iod": as_float(row["cell_iod"] if kind == "cell" else row["nuc_iod"]),
            "mean_od": as_float(row["cell_mean_od"] if kind == "cell" else row["nuc_mean_od"]),
        }
        repaired_metrics = measure_binary_metrics(
            repaired_binary,
            offset_y=local_slice[0].start,
            offset_x=local_slice[1].start,
        )
        if kind == "nucleus":
            repaired_metrics.update(
                measure_repaired_iod(
                    repaired_binary,
                    raw_path=clean_text(row["nucleus_source_image_path"]),
                    local_slice=local_slice,
                    i_bg=as_float(row["nuc_i_bg"]),
                )
            )
        else:
            repaired_metrics["iod"] = original_metrics["iod"]
            repaired_metrics["mean_od"] = original_metrics["mean_od"]

        repaired_path = path_output_for_mask(output_dir, kind, original_path)
        results.append(
            RepairResult(
                row_index=row_index,
                review_key=clean_text(row["review_key"]),
                kind=kind,
                applied=True,
                original_mask_path=original_path,
                repaired_mask_path=str(repaired_path.resolve()),
                original_area_px=original_metrics["area_px"],
                repaired_area_px=repaired_metrics["area_px"],
                original_area_um2=original_metrics["area_um2"],
                repaired_area_um2=repaired_metrics["area_um2"],
                original_centroid_x=original_metrics["centroid_x"],
                original_centroid_y=original_metrics["centroid_y"],
                repaired_centroid_x=repaired_metrics["centroid_x"],
                repaired_centroid_y=repaired_metrics["centroid_y"],
                original_solidity=original_metrics["solidity"],
                repaired_solidity=repaired_metrics["solidity"],
                original_circularity=original_metrics["circularity"],
                repaired_circularity=repaired_metrics["circularity"],
                original_iod=original_metrics["iod"],
                repaired_iod=repaired_metrics["iod"],
                original_mean_od=original_metrics["mean_od"],
                repaired_mean_od=repaired_metrics["mean_od"],
            )
        )

    repaired_path = path_output_for_mask(output_dir, kind, original_path)
    repaired_path.parent.mkdir(parents=True, exist_ok=True)
    tifffile.imwrite(repaired_path, updated, compression="zlib")
    return str(repaired_path.resolve()), results


def apply_manual_overrides(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    include = out["review_decision"].eq("keep") | (
        out["review_decision"].eq("repair") & out["review_repair_succeeded"].map(as_bool)
    )
    exclude = out["review_decision"].isin({"discard", "maybe"})
    include_ok = include & out["physical_pair_ok"] & out["one_to_one_cell"]
    out.loc[include_ok, ["keep_mask_pair", "keep_trim_5_95", "keep_strict_core"]] = True
    out.loc[exclude, ["keep_mask_pair", "keep_trim_5_95", "keep_strict_core", "keep_ultra_core"]] = False
    out["manual_report_include"] = np.where(
        include_ok,
        True,
        np.where(exclude, False, out["auto_keep_strict_core"]),
    )
    out["manual_report_decision_applied"] = np.where(
        out["review_decision"].isin({"keep", "repair", "discard", "maybe"}),
        True,
        False,
    )
    out["keep_strict_core"] = out["manual_report_include"]
    return out


def main() -> None:
    args = parse_args()
    require_exists(args.linked_pairs_csv)
    require_exists(args.cell_linkage_csv)
    args.output_dir.mkdir(parents=True, exist_ok=True)

    pairs = load_pairs(args.linked_pairs_csv)
    cell_linkage = load_cell_linkage(args.cell_linkage_csv)
    decisions = load_decisions(args.decisions_csv)

    pairs = pairs.merge(decisions, on="review_key", how="left")
    for col, default in {
        "review_decision": "unlabeled",
        "review_repair_mode": "",
        "review_note": "",
        "review_updated_at": "",
    }.items():
        if col not in pairs.columns:
            pairs[col] = default
        pairs[col] = pairs[col].fillna(default)

    if args.max_repairs > 0:
        repair_keys = list(
            pairs.loc[pairs["review_decision"].eq("repair"), "review_key"].drop_duplicates().head(args.max_repairs)
        )
        pairs = pairs[
            ~pairs["review_decision"].eq("repair") | pairs["review_key"].isin(repair_keys)
        ].copy()

    for col in [
        "keep_mask_pair",
        "keep_trim_5_95",
        "keep_strict_core",
        "keep_ultra_core",
        "physical_pair_ok",
        "flag_summary",
        "pair_core_distance",
    ]:
        if col in pairs.columns:
            pairs[f"original_{col}"] = pairs[col]

    pairs["cell_repair_applied"] = False
    pairs["nucleus_repair_applied"] = False
    pairs["cell_mask_path_original"] = pairs["cell_mask_path"]
    pairs["nucleus_mask_path_original"] = pairs["nucleus_mask_path"]

    cell_jobs: dict[str, list[tuple[int, pd.Series]]] = {}
    nucleus_jobs: dict[str, list[tuple[int, pd.Series]]] = {}
    for row_index, row in pairs.loc[pairs["review_decision"].eq("repair")].iterrows():
        mode = normalize_repair_mode(row["review_repair_mode"]) or "both"
        if mode in {"cell", "both"} and clean_text(row["cell_mask_path"]):
            cell_jobs.setdefault(clean_text(row["cell_mask_path"]), []).append((row_index, row))
        if mode in {"nucleus", "both"} and clean_text(row["nucleus_mask_path"]):
            nucleus_jobs.setdefault(clean_text(row["nucleus_mask_path"]), []).append((row_index, row))

    cell_path_map: dict[str, str] = {}
    nucleus_path_map: dict[str, str] = {}
    repair_results: list[RepairResult] = []

    for original_path, rows in cell_jobs.items():
        repaired_path, results = apply_repair_to_tile(
            kind="cell",
            original_path=original_path,
            tile_rows=rows,
            output_dir=args.output_dir,
        )
        cell_path_map[original_path] = repaired_path
        repair_results.extend(results)

    for original_path, rows in nucleus_jobs.items():
        repaired_path, results = apply_repair_to_tile(
            kind="nucleus",
            original_path=original_path,
            tile_rows=rows,
            output_dir=args.output_dir,
        )
        nucleus_path_map[original_path] = repaired_path
        repair_results.extend(results)

    if cell_path_map:
        pairs["cell_mask_path"] = pairs["cell_mask_path"].map(lambda p: cell_path_map.get(clean_text(p), clean_text(p)))
        cell_linkage["cell_mask_path"] = cell_linkage["cell_mask_path"].map(lambda p: cell_path_map.get(clean_text(p), clean_text(p)))
    if nucleus_path_map:
        pairs["nucleus_mask_path"] = pairs["nucleus_mask_path"].map(lambda p: nucleus_path_map.get(clean_text(p), clean_text(p)))

    for result in repair_results:
        idx = result.row_index
        if result.kind == "cell":
            pairs.at[idx, "cell_repair_applied"] = True
            pairs.at[idx, "cell_area_px"] = result.repaired_area_px
            pairs.at[idx, "cell_area_um2"] = round(result.repaired_area_um2, 4)
            pairs.at[idx, "cell_centroid_x"] = round(result.repaired_centroid_x, 2)
            pairs.at[idx, "cell_centroid_y"] = round(result.repaired_centroid_y, 2)
            pairs.at[idx, "cell_solidity"] = round(result.repaired_solidity, 6)
            pairs.at[idx, "cell_circularity"] = round(result.repaired_circularity, 6)
            pairs.at[idx, "cell_equiv_radius_px"] = round(math.sqrt(result.repaired_area_px / math.pi), 6)
            if "cell_object_id" in pairs.columns:
                cell_mask = cell_linkage["cell_object_id"].eq(pairs.at[idx, "cell_object_id"])
                cell_linkage.loc[cell_mask, "cell_area_px"] = result.repaired_area_px
                cell_linkage.loc[cell_mask, "cell_area_um2"] = round(result.repaired_area_um2, 4)
                cell_linkage.loc[cell_mask, "cell_centroid_x"] = round(result.repaired_centroid_x, 2)
                cell_linkage.loc[cell_mask, "cell_centroid_y"] = round(result.repaired_centroid_y, 2)
                cell_linkage.loc[cell_mask, "cell_solidity"] = round(result.repaired_solidity, 6)
                cell_linkage.loc[cell_mask, "cell_circularity"] = round(result.repaired_circularity, 6)
                cell_linkage.loc[cell_mask, "cell_equiv_radius_px"] = round(math.sqrt(result.repaired_area_px / math.pi), 6)
        else:
            pairs.at[idx, "nucleus_repair_applied"] = True
            pairs.at[idx, "nuc_area_px"] = result.repaired_area_px
            pairs.at[idx, "nuc_area_um2"] = round(result.repaired_area_um2, 4)
            pairs.at[idx, "nuc_centroid_x"] = round(result.repaired_centroid_x, 2)
            pairs.at[idx, "nuc_centroid_y"] = round(result.repaired_centroid_y, 2)
            pairs.at[idx, "nuc_iod"] = round(result.repaired_iod, 6)
            pairs.at[idx, "nuc_mean_od"] = round(result.repaired_mean_od, 6)

    pairs["review_repair_mode"] = pairs["review_repair_mode"].map(normalize_repair_mode)
    pairs["review_repair_mode"] = np.where(
        pairs["review_decision"].eq("repair") & pairs["review_repair_mode"].eq(""),
        "both",
        pairs["review_repair_mode"],
    )
    pairs["review_repair_succeeded"] = np.where(
        pairs["review_decision"].eq("repair"),
        np.select(
            [
                pairs["review_repair_mode"].eq("cell"),
                pairs["review_repair_mode"].eq("nucleus"),
                pairs["review_repair_mode"].eq("both"),
            ],
            [
                pairs["cell_repair_applied"],
                pairs["nucleus_repair_applied"],
                pairs["cell_repair_applied"] & pairs["nucleus_repair_applied"],
            ],
            default=False,
        ),
        False,
    )

    pairs["centroid_distance_px"] = np.sqrt(
        (pd.to_numeric(pairs["cell_centroid_x"], errors="coerce") - pd.to_numeric(pairs["nuc_centroid_x"], errors="coerce")) ** 2
        + (pd.to_numeric(pairs["cell_centroid_y"], errors="coerce") - pd.to_numeric(pairs["nuc_centroid_y"], errors="coerce")) ** 2
    )
    pairs["distance_over_cell_radius"] = pairs["centroid_distance_px"] / pd.to_numeric(pairs["cell_equiv_radius_px"], errors="coerce")
    pairs["distance_over_cell_radius"] = pairs["distance_over_cell_radius"].replace([np.inf, -np.inf], np.nan)
    pairs["nc_area_ratio"] = pd.to_numeric(pairs["nuc_area_um2"], errors="coerce") / pd.to_numeric(pairs["cell_area_um2"], errors="coerce")
    pairs["cytoplasm_area_um2"] = pd.to_numeric(pairs["cell_area_um2"], errors="coerce") - pd.to_numeric(pairs["nuc_area_um2"], errors="coerce")
    pairs["physical_pair_ok"] = (
        pairs["has_cell_match"].map(as_bool)
        & pd.to_numeric(pairs["cell_area_um2"], errors="coerce").gt(0)
        & pd.to_numeric(pairs["nuc_area_um2"], errors="coerce").gt(0)
        & pd.to_numeric(pairs["cell_area_um2"], errors="coerce").gt(pd.to_numeric(pairs["nuc_area_um2"], errors="coerce"))
        & pd.to_numeric(pairs["cytoplasm_area_um2"], errors="coerce").gt(0)
    )

    auto_qc = annotate_pair_qc(pairs, cell_linkage)
    for col in ["keep_mask_pair", "keep_trim_5_95", "keep_strict_core", "keep_ultra_core", "flag_summary", "pair_core_distance"]:
        auto_qc[f"auto_{col}"] = auto_qc[col]
    patched = apply_manual_overrides(auto_qc)

    repair_manifest = pd.DataFrame([result.__dict__ for result in repair_results])
    if repair_manifest.empty:
        repair_manifest = pd.DataFrame(
            columns=[
                "row_index",
                "review_key",
                "kind",
                "applied",
                "original_mask_path",
                "repaired_mask_path",
                "original_area_px",
                "repaired_area_px",
                "original_area_um2",
                "repaired_area_um2",
                "original_centroid_x",
                "original_centroid_y",
                "repaired_centroid_x",
                "repaired_centroid_y",
                "original_solidity",
                "repaired_solidity",
                "original_circularity",
                "repaired_circularity",
                "original_iod",
                "repaired_iod",
                "original_mean_od",
                "repaired_mean_od",
            ]
        )

    reviewed_pairs_csv = args.output_dir / "linked_nucleus_pairs_reviewed.csv.gz"
    patched.to_csv(reviewed_pairs_csv, index=False)
    reviewed_cell_linkage_csv = args.output_dir / "cell_linkage_summary_reviewed.csv.gz"
    cell_linkage.to_csv(reviewed_cell_linkage_csv, index=False)
    repair_manifest_csv = args.output_dir / "repair_manifest.csv"
    repair_manifest.to_csv(repair_manifest_csv, index=False)

    summary = {
        "linked_pairs_csv": str(args.linked_pairs_csv.resolve()),
        "decisions_csv": str(args.decisions_csv.resolve()),
        "n_pairs_total": int(len(patched)),
        "n_reviewed_pairs": int(patched["review_decision"].ne("unlabeled").sum()),
        "n_repair_decisions": int(patched["review_decision"].eq("repair").sum()),
        "n_discard_decisions": int(patched["review_decision"].eq("discard").sum()),
        "n_keep_decisions": int(patched["review_decision"].eq("keep").sum()),
        "n_maybe_decisions": int(patched["review_decision"].eq("maybe").sum()),
        "n_cell_repairs_applied": int(patched["cell_repair_applied"].sum()),
        "n_nucleus_repairs_applied": int(patched["nucleus_repair_applied"].sum()),
        "n_successful_repair_decisions": int(patched["review_repair_succeeded"].sum()),
        "n_effective_keep_strict_core": int(patched["keep_strict_core"].sum()),
        "reviewed_pairs_output": str(reviewed_pairs_csv.resolve()),
        "reviewed_cell_linkage_output": str(reviewed_cell_linkage_csv.resolve()),
        "repair_manifest_output": str(repair_manifest_csv.resolve()),
    }
    (args.output_dir / "summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
