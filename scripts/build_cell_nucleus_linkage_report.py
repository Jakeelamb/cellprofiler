#!/usr/bin/env python3
"""Link brightfield cell masks to nucleus/IOD measurements and build a QC report."""

from __future__ import annotations

import argparse
import html
import json
import zipfile
from collections import OrderedDict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import tifffile
from roifile import ImagejRoi
from scipy import stats as sp_stats
from skimage.draw import polygon

PROJECT = Path(__file__).resolve().parent.parent
DEFAULT_CELL_CSV = PROJECT / "output" / "runs" / "full_dataset_v1" / "cell_size_segmentation" / "all_measurements.csv"
ENRICHED_NUCLEUS_SNAPSHOT = PROJECT / "output" / "cell_nucleus_linkage" / "threshold_tuned_v1_nucleus_measurements_enriched.csv.gz"
DEFAULT_NUCLEUS_CSV = (
    ENRICHED_NUCLEUS_SNAPSHOT
    if ENRICHED_NUCLEUS_SNAPSHOT.exists()
    else PROJECT / "output" / "nucleus_iod" / "nucleus_iod_measurements.csv"
)
DEFAULT_OUTPUT_DIR = PROJECT / "output" / "cell_nucleus_linkage"
TILE_SIZE = 4096
TRIM_LO = 5
TRIM_HI = 95
STRICT_THRESHOLDS = {
    "strict_core": {"z": 2.5, "combined": 4.5},
    "ultra_core": {"z": 2.0, "combined": 3.5},
}
WINDOW_RADIUS = 1
MASK_CACHE_SIZE = 8


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cell-csv", type=Path, default=DEFAULT_CELL_CSV)
    parser.add_argument("--nucleus-csv", type=Path, default=DEFAULT_NUCLEUS_CSV)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--image-type", default="brightfield", choices=["brightfield"])
    parser.add_argument(
        "--window-radius",
        type=int,
        default=WINDOW_RADIUS,
        help="Fallback half-width for unique-label mask sampling if the rounded centroid lands outside a cell label.",
    )
    return parser.parse_args()


def require_exists(path: Path) -> None:
    if not path.exists():
        raise FileNotFoundError(path)


def load_measurements(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, keep_default_na=False)


def normalize_cells(df: pd.DataFrame, image_type: str) -> pd.DataFrame:
    out = df.copy()
    out = out[out["image_type"] == image_type].copy()
    numeric_cols = [
        "label", "area_px", "area_um2", "solidity", "circularity", "iod", "mean_od",
        "centroid_y", "centroid_x", "i_bg", "tile_y0", "tile_x0", "tile_h", "tile_w",
        "tile_score", "mask_label_id",
    ]
    for col in numeric_cols:
        if col in out.columns:
            out[col] = pd.to_numeric(out[col], errors="coerce")
    out["tile_name"] = out["tile_name"].astype(str)
    out["mask_path"] = out["mask_path"].astype(str)
    out["raw_mask_path"] = out["raw_mask_path"].astype(str)
    if "tile_manifest_path" not in out.columns:
        out["tile_manifest_path"] = ""
    if "source_image_path" not in out.columns:
        out["source_image_path"] = ""
    if "run_manifest_path" not in out.columns:
        out["run_manifest_path"] = ""
    out["tile_manifest_path"] = out["tile_manifest_path"].astype(str)
    out["source_image_path"] = out["source_image_path"].astype(str)
    out["run_manifest_path"] = out["run_manifest_path"].astype(str)
    out = out.rename(columns={"mask_path": "cell_mask_path", "raw_mask_path": "cell_raw_mask_path"})
    out = out.rename(
        columns={
            "tile_manifest_path": "cell_tile_manifest_path",
            "source_image_path": "cell_source_image_path",
            "run_manifest_path": "cell_run_manifest_path",
        }
    )
    out["cell_object_id"] = (
        out["filename"].astype(str)
        + "::"
        + out["tile_name"].astype(str)
        + "::"
        + out["mask_label_id"].fillna(-1).astype(int).astype(str)
    )
    out["cell_equiv_radius_px"] = np.sqrt(out["area_px"].clip(lower=0) / np.pi)
    return out.reset_index(drop=True)


def normalize_nuclei(df: pd.DataFrame, image_type: str) -> pd.DataFrame:
    out = df.copy()
    out = out[out["image_type"] == image_type].copy()
    numeric_cols = [
        "label", "area_px", "area_um2", "iod", "mean_od", "centroid_x", "centroid_y", "i_bg",
        "tile_y0", "tile_x0", "tile_height_px", "tile_width_px",
    ]
    for col in numeric_cols:
        if col in out.columns:
            out[col] = pd.to_numeric(out[col], errors="coerce")

    if "tile_y0" not in out.columns or out["tile_y0"].isna().all():
        out["tile_y0"] = (np.floor(out["centroid_y"] / TILE_SIZE) * TILE_SIZE).astype(int)
    else:
        out["tile_y0"] = out["tile_y0"].fillna(np.floor(out["centroid_y"] / TILE_SIZE) * TILE_SIZE).astype(int)

    if "tile_x0" not in out.columns or out["tile_x0"].isna().all():
        out["tile_x0"] = (np.floor(out["centroid_x"] / TILE_SIZE) * TILE_SIZE).astype(int)
    else:
        out["tile_x0"] = out["tile_x0"].fillna(np.floor(out["centroid_x"] / TILE_SIZE) * TILE_SIZE).astype(int)

    if "tile_height_px" not in out.columns:
        out["tile_height_px"] = TILE_SIZE
    out["tile_height_px"] = pd.to_numeric(out["tile_height_px"], errors="coerce").fillna(TILE_SIZE).astype(int)

    if "tile_width_px" not in out.columns:
        out["tile_width_px"] = TILE_SIZE
    out["tile_width_px"] = pd.to_numeric(out["tile_width_px"], errors="coerce").fillna(TILE_SIZE).astype(int)

    if "tile_name" not in out.columns or out["tile_name"].astype(str).eq("").all():
        out["tile_name"] = (
            "tile_y"
            + out["tile_y0"].astype(int).astype(str).str.zfill(6)
            + "_x"
            + out["tile_x0"].astype(int).astype(str).str.zfill(6)
            + ".tiff"
        )
    else:
        empty = out["tile_name"].astype(str).eq("")
        out.loc[empty, "tile_name"] = (
            "tile_y"
            + out.loc[empty, "tile_y0"].astype(int).astype(str).str.zfill(6)
            + "_x"
            + out.loc[empty, "tile_x0"].astype(int).astype(str).str.zfill(6)
            + ".tiff"
        )

    out["nucleus_object_id"] = (
        out["filename"].astype(str)
        + "::"
        + out["tile_name"].astype(str)
        + "::"
        + out["label"].fillna(-1).astype(int).astype(str)
    )
    if "mask_path" in out.columns:
        out = out.rename(columns={"mask_path": "nucleus_mask_path"})
    else:
        out["nucleus_mask_path"] = ""
    if "roi_zip_path" not in out.columns:
        out["roi_zip_path"] = ""
    if "tile_manifest_path" not in out.columns:
        out["tile_manifest_path"] = ""
    if "source_image_path" not in out.columns:
        out["source_image_path"] = ""
    if "run_manifest_path" not in out.columns:
        out["run_manifest_path"] = ""
    if "raw_imagej_results_path" not in out.columns:
        out["raw_imagej_results_path"] = ""
    out["nucleus_mask_path"] = out["nucleus_mask_path"].astype(str)
    out["roi_zip_path"] = out["roi_zip_path"].astype(str)
    out["tile_manifest_path"] = out["tile_manifest_path"].astype(str)
    out["source_image_path"] = out["source_image_path"].astype(str)
    out["run_manifest_path"] = out["run_manifest_path"].astype(str)
    out["raw_imagej_results_path"] = out["raw_imagej_results_path"].astype(str)
    out = out.rename(
        columns={
            "tile_manifest_path": "nucleus_tile_manifest_path",
            "source_image_path": "nucleus_source_image_path",
            "run_manifest_path": "nucleus_run_manifest_path",
        }
    )
    return out.reset_index(drop=True)


class MaskCache:
    def __init__(self, max_size: int = MASK_CACHE_SIZE) -> None:
        self.max_size = max_size
        self._store: OrderedDict[str, np.ndarray] = OrderedDict()

    def get(self, path: str) -> np.ndarray:
        key = str(path)
        arr = self._store.pop(key, None)
        if arr is None:
            arr = tifffile.imread(key)
        self._store[key] = arr
        while len(self._store) > self.max_size:
            self._store.popitem(last=False)
        return arr


class RoiZipCache:
    def __init__(self, max_size: int = MASK_CACHE_SIZE) -> None:
        self.max_size = max_size
        self._store: OrderedDict[str, dict[int, np.ndarray]] = OrderedDict()

    def get(self, path: str) -> dict[int, np.ndarray]:
        key = str(path)
        mapping = self._store.pop(key, None)
        if mapping is None:
            mapping = {}
            with zipfile.ZipFile(key) as handle:
                for name in sorted(handle.namelist()):
                    try:
                        roi_idx = int(name.split("-", 1)[0])
                    except ValueError:
                        continue
                    roi = ImagejRoi.frombytes(handle.read(name))
                    mapping[roi_idx] = roi.coordinates()
        self._store[key] = mapping
        while len(self._store) > self.max_size:
            self._store.popitem(last=False)
        return mapping


def trim_mask(series: pd.Series) -> pd.Series:
    if len(series) < 10:
        return pd.Series(True, index=series.index)
    lo, hi = np.percentile(series.to_numpy(dtype=float), [TRIM_LO, TRIM_HI])
    return (series >= lo) & (series <= hi)


def robust_z_scores(values: np.ndarray) -> np.ndarray:
    vals = np.asarray(values, dtype=float)
    if len(vals) < 3:
        return np.zeros(len(vals), dtype=float)
    med = np.median(vals)
    mad = np.median(np.abs(vals - med))
    if mad == 0:
        return np.zeros(len(vals), dtype=float)
    return 0.6745 * (vals - med) / mad


def lookup_single_label_window(mask: np.ndarray, y: int, x: int, radius: int) -> tuple[int, int]:
    y0 = max(0, y - radius)
    y1 = min(mask.shape[0], y + radius + 1)
    x0 = max(0, x - radius)
    x1 = min(mask.shape[1], x + radius + 1)
    window = mask[y0:y1, x0:x1]
    positive = window[window > 0]
    if positive.size == 0:
        return 0, 0
    labels, counts = np.unique(positive, return_counts=True)
    if len(labels) == 1:
        return int(labels[0]), int(counts[0])
    return 0, int(counts.max())


def roi_overlap_label(cell_mask: np.ndarray, coords: np.ndarray) -> tuple[int, int, int, float, int]:
    rr, cc = polygon(coords[:, 1], coords[:, 0], shape=cell_mask.shape)
    roi_px = int(len(rr))
    if roi_px == 0:
        return 0, 0, 0, 0.0, 0
    vals = cell_mask[rr, cc]
    positive = vals[vals > 0]
    if positive.size == 0:
        return 0, roi_px, 0, 0.0, 0
    labels, counts = np.unique(positive, return_counts=True)
    best_idx = int(np.argmax(counts))
    best_label = int(labels[best_idx])
    best_px = int(counts[best_idx])
    return (
        best_label if len(labels) == 1 else 0,
        roi_px,
        best_px,
        float(best_px / roi_px) if roi_px else 0.0,
        int(len(labels)),
    )


def mask_overlap_label(
    cell_mask: np.ndarray,
    nucleus_mask: np.ndarray,
    nucleus_label: int,
) -> tuple[int, int, int, float, int]:
    nucleus_binary = nucleus_mask == int(nucleus_label)
    nucleus_px = int(nucleus_binary.sum())
    if nucleus_px == 0:
        return 0, 0, 0, 0.0, 0
    vals = cell_mask[nucleus_binary]
    positive = vals[vals > 0]
    if positive.size == 0:
        return 0, nucleus_px, 0, 0.0, 0
    labels, counts = np.unique(positive, return_counts=True)
    best_idx = int(np.argmax(counts))
    best_label = int(labels[best_idx])
    best_px = int(counts[best_idx])
    return (
        best_label if len(labels) == 1 else 0,
        nucleus_px,
        best_px,
        float(best_px / nucleus_px) if nucleus_px else 0.0,
        int(len(labels)),
    )


def build_tile_metadata(cells: pd.DataFrame) -> pd.DataFrame:
    tile_cols = ["filename", "tile_name", "cell_mask_path", "tile_y0", "tile_x0", "tile_h", "tile_w"]
    tiles = cells[tile_cols].drop_duplicates().copy()
    tiles["cell_mask_path"] = tiles["cell_mask_path"].astype(str)
    tiles = tiles.rename(
        columns={
            "tile_y0": "cell_tile_y0",
            "tile_x0": "cell_tile_x0",
            "tile_h": "cell_tile_h",
            "tile_w": "cell_tile_w",
        }
    )
    return tiles


def link_nuclei_to_cells(
    cells: pd.DataFrame,
    nuclei: pd.DataFrame,
    window_radius: int,
) -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    tile_meta = build_tile_metadata(cells)
    merged = nuclei.merge(tile_meta, on=["filename", "tile_name"], how="left", validate="many_to_one")
    merged["has_cell_tile"] = merged["cell_mask_path"].notna() & merged["cell_mask_path"].astype(str).ne("")
    merged["has_nucleus_mask"] = merged["nucleus_mask_path"].notna() & merged["nucleus_mask_path"].astype(str).ne("")
    merged["has_nucleus_roi"] = merged["roi_zip_path"].notna() & merged["roi_zip_path"].astype(str).ne("")
    merged["mask_sample_label"] = 0
    merged["window_label_count"] = 0
    merged["link_method"] = "unmatched"
    merged["nucleus_mask_pixel_count"] = 0
    merged["nucleus_mask_best_overlap_px"] = 0
    merged["nucleus_mask_best_overlap_fraction"] = 0.0
    merged["nucleus_mask_overlap_label_count"] = 0
    merged["roi_pixel_count"] = 0
    merged["roi_best_overlap_px"] = 0
    merged["roi_best_overlap_fraction"] = 0.0
    merged["roi_overlap_label_count"] = 0

    mask_cache = MaskCache()
    roi_cache = RoiZipCache()
    mask_unique_hits = 0
    roi_unique_hits = 0
    direct_hits = 0
    window_hits = 0

    mask_groups = merged[merged["has_cell_tile"] & merged["has_nucleus_mask"]].groupby(
        ["cell_mask_path", "nucleus_mask_path"]
    )
    for (cell_mask_path, nucleus_mask_path), idx in mask_groups.groups.items():
        cell_mask = mask_cache.get(cell_mask_path)
        nucleus_mask = mask_cache.get(nucleus_mask_path)
        for row_idx in idx:
            try:
                nucleus_label = int(merged.at[row_idx, "label"])
            except (TypeError, ValueError):
                continue
            label, nucleus_px, best_px, best_frac, label_count = mask_overlap_label(
                cell_mask,
                nucleus_mask,
                nucleus_label,
            )
            merged.at[row_idx, "nucleus_mask_pixel_count"] = nucleus_px
            merged.at[row_idx, "nucleus_mask_best_overlap_px"] = best_px
            merged.at[row_idx, "nucleus_mask_best_overlap_fraction"] = best_frac
            merged.at[row_idx, "nucleus_mask_overlap_label_count"] = label_count
            if label > 0:
                mask_unique_hits += 1
                merged.at[row_idx, "mask_sample_label"] = label
                merged.at[row_idx, "link_method"] = "nucleus_mask_overlap_unique"

    roi_groups = merged[
        merged["has_cell_tile"] & merged["has_nucleus_roi"] & merged["mask_sample_label"].eq(0)
    ].groupby(
        ["cell_mask_path", "roi_zip_path"]
    )
    for (cell_mask_path, roi_zip_path), idx in roi_groups.groups.items():
        cell_mask = mask_cache.get(cell_mask_path)
        roi_map = roi_cache.get(roi_zip_path)
        for row_idx in idx:
            try:
                nucleus_label = int(merged.at[row_idx, "label"])
            except (TypeError, ValueError):
                continue
            coords = roi_map.get(nucleus_label)
            if coords is None or len(coords) == 0:
                continue
            label, roi_px, best_px, best_frac, label_count = roi_overlap_label(cell_mask, coords)
            merged.at[row_idx, "roi_pixel_count"] = roi_px
            merged.at[row_idx, "roi_best_overlap_px"] = best_px
            merged.at[row_idx, "roi_best_overlap_fraction"] = best_frac
            merged.at[row_idx, "roi_overlap_label_count"] = label_count
            if label > 0:
                roi_unique_hits += 1
                merged.at[row_idx, "mask_sample_label"] = label
                merged.at[row_idx, "link_method"] = "roi_overlap_unique"

    grouped = merged[(merged["has_cell_tile"]) & (merged["mask_sample_label"] == 0)].groupby("cell_mask_path")
    for cell_mask_path, idx in grouped.groups.items():
        mask = mask_cache.get(cell_mask_path)
        sub = merged.loc[idx, ["centroid_y", "centroid_x", "tile_y0", "tile_x0"]].copy()
        local_y = np.rint(sub["centroid_y"].to_numpy(dtype=float) - sub["tile_y0"].to_numpy(dtype=float)).astype(int)
        local_x = np.rint(sub["centroid_x"].to_numpy(dtype=float) - sub["tile_x0"].to_numpy(dtype=float)).astype(int)
        in_bounds = (
            (local_y >= 0)
            & (local_y < mask.shape[0])
            & (local_x >= 0)
            & (local_x < mask.shape[1])
        )

        labels = np.zeros(len(idx), dtype=np.int64)
        if in_bounds.any():
            labels[in_bounds] = mask[local_y[in_bounds], local_x[in_bounds]].astype(np.int64)
        hit_mask = labels > 0
        if hit_mask.any():
            direct_hits += int(hit_mask.sum())
            merged.loc[idx[hit_mask], "mask_sample_label"] = labels[hit_mask]
            merged.loc[idx[hit_mask], "link_method"] = "mask_centroid"

        miss_positions = np.where(~hit_mask & in_bounds)[0]
        if len(miss_positions) and window_radius > 0:
            for pos in miss_positions:
                label, count = lookup_single_label_window(mask, local_y[pos], local_x[pos], window_radius)
                if label > 0:
                    window_hits += 1
                    merged.loc[idx[pos], "mask_sample_label"] = label
                    merged.loc[idx[pos], "window_label_count"] = count
                    merged.loc[idx[pos], "link_method"] = "mask_window_unique"

    cell_lookup = cells[
        [
            "cell_object_id", "filename", "tile_name", "cell_mask_path", "mask_label_id",
            "area_px", "area_um2", "solidity", "circularity", "iod", "mean_od",
            "centroid_y", "centroid_x", "i_bg", "cell_equiv_radius_px",
            "cell_tile_manifest_path", "cell_source_image_path", "cell_run_manifest_path",
        ]
    ].copy()
    cell_lookup = cell_lookup.rename(
        columns={
            "area_px": "cell_area_px",
            "area_um2": "cell_area_um2",
            "solidity": "cell_solidity",
            "circularity": "cell_circularity",
            "iod": "cell_iod",
            "mean_od": "cell_mean_od",
            "centroid_y": "cell_centroid_y",
            "centroid_x": "cell_centroid_x",
            "i_bg": "cell_i_bg",
        }
    )

    linked = merged.merge(
        cell_lookup,
        left_on=["filename", "tile_name", "cell_mask_path", "mask_sample_label"],
        right_on=["filename", "tile_name", "cell_mask_path", "mask_label_id"],
        how="left",
    )

    linked = linked.rename(
        columns={
            "area_px": "nuc_area_px",
            "area_um2": "nuc_area_um2",
            "iod": "nuc_iod",
            "mean_od": "nuc_mean_od",
            "centroid_x": "nuc_centroid_x",
            "centroid_y": "nuc_centroid_y",
            "i_bg": "nuc_i_bg",
            "label": "nucleus_label",
        }
    )
    linked["has_cell_match"] = linked["cell_object_id"].notna() & linked["cell_object_id"].astype(str).ne("")

    linked["centroid_distance_px"] = np.sqrt(
        (linked["cell_centroid_x"] - linked["nuc_centroid_x"]) ** 2
        + (linked["cell_centroid_y"] - linked["nuc_centroid_y"]) ** 2
    )
    linked.loc[~linked["has_cell_match"], "centroid_distance_px"] = np.nan
    linked["distance_over_cell_radius"] = linked["centroid_distance_px"] / linked["cell_equiv_radius_px"]
    linked["distance_over_cell_radius"] = linked["distance_over_cell_radius"].replace([np.inf, -np.inf], np.nan)
    linked["nc_area_ratio"] = linked["nuc_area_um2"] / linked["cell_area_um2"]
    linked["cytoplasm_area_um2"] = linked["cell_area_um2"] - linked["nuc_area_um2"]
    linked["physical_pair_ok"] = (
        linked["has_cell_match"]
        & linked["cell_area_um2"].gt(0)
        & linked["nuc_area_um2"].gt(0)
        & linked["cell_area_um2"].gt(linked["nuc_area_um2"])
        & linked["cytoplasm_area_um2"].gt(0)
    )

    cell_linkage = cell_lookup.merge(
        linked.groupby("cell_object_id")
        .agg(
            nuclei_hitting_cell=("nucleus_object_id", "nunique"),
            primary_nucleus_object_id=("nucleus_object_id", lambda s: next(iter(pd.unique(s)), "")),
            primary_link_method=("link_method", lambda s: next(iter(pd.unique(s)), "")),
        )
        .reset_index(),
        on="cell_object_id",
        how="left",
    )
    cell_linkage["nuclei_hitting_cell"] = cell_linkage["nuclei_hitting_cell"].fillna(0).astype(int)
    cell_linkage["primary_nucleus_object_id"] = cell_linkage["primary_nucleus_object_id"].fillna("")
    cell_linkage["primary_link_method"] = cell_linkage["primary_link_method"].fillna("")
    cell_linkage["cell_linkage_status"] = np.select(
        [
            cell_linkage["nuclei_hitting_cell"].eq(0),
            cell_linkage["nuclei_hitting_cell"].eq(1),
            cell_linkage["nuclei_hitting_cell"].gt(1),
        ],
        ["no_nucleus", "one_to_one", "multi_nucleus"],
        default="no_nucleus",
    )

    summary = {
        "n_cells_total": int(len(cells)),
        "n_nuclei_total": int(len(nuclei)),
        "n_nuclei_with_candidate_tile": int(linked["has_cell_tile"].sum()),
        "n_nuclei_with_mask_artifact": int(linked["has_nucleus_mask"].sum()),
        "n_nuclei_with_roi_artifact": int(linked["has_nucleus_roi"].sum()),
        "n_nucleus_mask_unique_hits": int(mask_unique_hits),
        "n_roi_overlap_unique_hits": int(roi_unique_hits),
        "n_direct_mask_hits": int(direct_hits),
        "n_window_mask_hits": int(window_hits),
        "n_matched_nuclei": int(linked["has_cell_match"].sum()),
        "n_unmatched_nuclei": int((~linked["has_cell_match"]).sum()),
        "n_cells_one_to_one": int((cell_linkage["cell_linkage_status"] == "one_to_one").sum()),
        "n_cells_multi_nucleus": int((cell_linkage["cell_linkage_status"] == "multi_nucleus").sum()),
        "n_cells_without_nucleus": int((cell_linkage["cell_linkage_status"] == "no_nucleus").sum()),
    }
    return linked, cell_linkage, summary


def annotate_pair_qc(linked: pd.DataFrame, cell_linkage: pd.DataFrame) -> pd.DataFrame:
    one_to_one_cells = set(cell_linkage.loc[cell_linkage["cell_linkage_status"] == "one_to_one", "cell_object_id"])
    out = linked.copy()
    out["one_to_one_cell"] = out["cell_object_id"].isin(one_to_one_cells)
    out["keep_mask_pair"] = out["physical_pair_ok"] & out["one_to_one_cell"]
    out["flag_summary"] = ""

    valid = out["keep_mask_pair"].copy()
    out["keep_trim_5_95"] = False
    out["keep_strict_core"] = False
    out["keep_ultra_core"] = False
    out["z_log_cell_area"] = np.nan
    out["z_log_nuc_area"] = np.nan
    out["z_log_nuc_iod"] = np.nan
    out["z_log_nc_ratio"] = np.nan
    out["pair_core_distance"] = np.nan

    frames = []
    for _, grp in out.groupby(["filename", "species"], sort=False):
        g = grp.copy()
        eligible = g["keep_mask_pair"]
        if not eligible.any():
            frames.append(g)
            continue

        sub = g.loc[eligible].copy()
        trim_cell = trim_mask(sub["cell_area_um2"])
        trim_nuc = trim_mask(sub["nuc_area_um2"])
        trim_iod = trim_mask(sub["nuc_iod"])
        trim_ratio = trim_mask(sub["nc_area_ratio"])
        keep_trim = trim_cell & trim_nuc & trim_iod & trim_ratio

        log_cell = np.log(sub["cell_area_um2"].to_numpy(dtype=float))
        log_nuc = np.log(sub["nuc_area_um2"].to_numpy(dtype=float))
        log_iod = np.log(sub["nuc_iod"].to_numpy(dtype=float))
        log_ratio = np.log(sub["nc_area_ratio"].to_numpy(dtype=float))

        z_cell = robust_z_scores(log_cell)
        z_nuc = robust_z_scores(log_nuc)
        z_iod = robust_z_scores(log_iod)
        z_ratio = robust_z_scores(log_ratio)
        combined = np.sqrt(z_cell ** 2 + z_nuc ** 2 + z_iod ** 2 + z_ratio ** 2)

        sub["keep_trim_5_95"] = keep_trim.to_numpy()
        sub["z_log_cell_area"] = z_cell
        sub["z_log_nuc_area"] = z_nuc
        sub["z_log_nuc_iod"] = z_iod
        sub["z_log_nc_ratio"] = z_ratio
        sub["pair_core_distance"] = combined

        for tier_name, spec in STRICT_THRESHOLDS.items():
            keep = (
                sub["keep_trim_5_95"].to_numpy()
                & (np.abs(z_cell) <= spec["z"])
                & (np.abs(z_nuc) <= spec["z"])
                & (np.abs(z_iod) <= spec["z"])
                & (np.abs(z_ratio) <= spec["z"])
                & (combined <= spec["combined"])
            )
            sub[f"keep_{tier_name}"] = keep

        g.loc[sub.index, [
            "keep_trim_5_95", "keep_strict_core", "keep_ultra_core",
            "z_log_cell_area", "z_log_nuc_area", "z_log_nuc_iod",
            "z_log_nc_ratio", "pair_core_distance",
        ]] = sub[
            [
                "keep_trim_5_95", "keep_strict_core", "keep_ultra_core",
                "z_log_cell_area", "z_log_nuc_area", "z_log_nuc_iod",
                "z_log_nc_ratio", "pair_core_distance",
            ]
        ]
        frames.append(g)

    out = pd.concat(frames, ignore_index=True)

    issues = []
    for _, row in out.iterrows():
        flags = []
        if not row["has_cell_match"]:
            flags.append("no_cell_match")
        if row["has_cell_match"] and not row["one_to_one_cell"]:
            flags.append("multi_nucleus_cell")
        if row["has_cell_match"] and not row["physical_pair_ok"]:
            flags.append("physical_pair_fail")
        if row["keep_mask_pair"] and not row["keep_trim_5_95"]:
            flags.append("trim_outlier")
        if row["keep_trim_5_95"] and not row["keep_strict_core"]:
            flags.append("strict_core_outlier")
        if row["keep_strict_core"] and not row["keep_ultra_core"]:
            flags.append("ultra_core_outlier")
        issues.append("; ".join(flags))
    out["flag_summary"] = issues
    return out


def _corr_pair(df: pd.DataFrame, x: str, y: str, method: str) -> float:
    sub = df[[x, y]].dropna()
    if len(sub) < 3:
        return np.nan
    if method == "pearson":
        return float(sp_stats.pearsonr(sub[x], sub[y]).statistic)
    return float(sp_stats.spearmanr(sub[x], sub[y]).statistic)


def fisher_mean_correlation(values: pd.Series, weights: pd.Series) -> float:
    vals = pd.to_numeric(values, errors="coerce").to_numpy(dtype=float)
    wts = pd.to_numeric(weights, errors="coerce").to_numpy(dtype=float)
    keep = np.isfinite(vals) & np.isfinite(wts) & (wts > 0)
    if not keep.any():
        return np.nan
    clipped = np.clip(vals[keep], -0.999999, 0.999999)
    return float(np.tanh(np.average(np.arctanh(clipped), weights=wts[keep])))


def summarize_images(pairs: pd.DataFrame, cells: pd.DataFrame, nuclei: pd.DataFrame) -> pd.DataFrame:
    cell_counts = cells.groupby("filename").size().rename("n_cells_total")
    nuc_counts = nuclei.groupby("filename").size().rename("n_nuclei_total")
    rows = []

    for (filename, species), grp in pairs.groupby(["filename", "species"], sort=False):
        strict = grp[grp["keep_strict_core"]].copy()
        row = {
            "filename": filename,
            "species": species,
            "slide_id": grp["slide_id"].iloc[0],
            "specimen_id": grp["specimen_id"].iloc[0],
            "image_type": grp["image_type"].iloc[0],
            "n_cells_total": int(cell_counts.get(filename, 0)),
            "n_nuclei_total": int(nuc_counts.get(filename, 0)),
            "n_nuclei_with_cell_match": int(grp["has_cell_match"].sum()),
            "n_pairs_mask": int(grp["keep_mask_pair"].sum()),
            "n_pairs_trim_5_95": int(grp["keep_trim_5_95"].sum()),
            "n_pairs_strict_core": int(grp["keep_strict_core"].sum()),
            "n_pairs_ultra_core": int(grp["keep_ultra_core"].sum()),
        }
        row["nucleus_to_cell_match_rate"] = (
            row["n_nuclei_with_cell_match"] / row["n_nuclei_total"] if row["n_nuclei_total"] else np.nan
        )
        row["strict_pair_rate_from_cells"] = (
            row["n_pairs_strict_core"] / row["n_cells_total"] if row["n_cells_total"] else np.nan
        )
        row["strict_pair_rate_from_nuclei"] = (
            row["n_pairs_strict_core"] / row["n_nuclei_total"] if row["n_nuclei_total"] else np.nan
        )
        row["median_cell_area_um2_strict"] = strict["cell_area_um2"].median() if len(strict) else np.nan
        row["median_nuc_area_um2_strict"] = strict["nuc_area_um2"].median() if len(strict) else np.nan
        row["median_nuc_iod_strict"] = strict["nuc_iod"].median() if len(strict) else np.nan
        row["median_nc_ratio_strict"] = strict["nc_area_ratio"].median() if len(strict) else np.nan
        row["pearson_cell_vs_nuc_area_strict"] = _corr_pair(strict, "cell_area_um2", "nuc_area_um2", "pearson")
        row["spearman_cell_vs_nuc_area_strict"] = _corr_pair(strict, "cell_area_um2", "nuc_area_um2", "spearman")
        row["pearson_cell_vs_nuc_iod_strict"] = _corr_pair(strict, "cell_area_um2", "nuc_iod", "pearson")
        row["spearman_cell_vs_nuc_iod_strict"] = _corr_pair(strict, "cell_area_um2", "nuc_iod", "spearman")
        row["pearson_nuc_area_vs_iod_strict"] = _corr_pair(strict, "nuc_area_um2", "nuc_iod", "pearson")
        row["spearman_nuc_area_vs_iod_strict"] = _corr_pair(strict, "nuc_area_um2", "nuc_iod", "spearman")
        rows.append(row)

    return pd.DataFrame(rows).sort_values(["species", "filename"]).reset_index(drop=True)


def annotate_image_qc(image_summary: pd.DataFrame) -> tuple[pd.DataFrame, dict]:
    out = image_summary.copy()
    valid = out[out["n_pairs_strict_core"] > 0].copy()
    pair_rate_threshold = float(valid["strict_pair_rate_from_cells"].quantile(0.75)) if len(valid) else np.nan
    spearman_threshold = 0.10

    out["flag_low_pair_rate"] = False
    if np.isfinite(pair_rate_threshold):
        out["flag_low_pair_rate"] = out["strict_pair_rate_from_cells"] < pair_rate_threshold

    out["flag_weak_cell_nucleus_correlation"] = out["spearman_cell_vs_nuc_area_strict"] < spearman_threshold
    out["analysis_ready_image"] = (
        out["n_pairs_strict_core"].gt(0)
        & (~out["flag_low_pair_rate"])
        & (~out["flag_weak_cell_nucleus_correlation"])
    )
    return out, {
        "analysis_ready_pair_rate_threshold": pair_rate_threshold,
        "analysis_ready_spearman_threshold": spearman_threshold,
    }


def summarize_species(image_summary: pd.DataFrame, pairs: pd.DataFrame) -> pd.DataFrame:
    strict_pairs = pairs[pairs["keep_strict_core"]].copy()
    rows = []
    for species, grp in image_summary.groupby("species", sort=True):
        strict_species = strict_pairs[strict_pairs["species"] == species]
        row = {
            "species": species,
            "n_images": int(grp["filename"].nunique()),
            "n_specimens": int(grp["specimen_id"].nunique()),
            "n_analysis_ready_images": int(grp["analysis_ready_image"].sum()),
            "n_cells_total": int(grp["n_cells_total"].sum()),
            "n_nuclei_total": int(grp["n_nuclei_total"].sum()),
            "n_pairs_strict_core": int(grp["n_pairs_strict_core"].sum()),
            "n_pairs_strict_core_analysis_ready": int(grp.loc[grp["analysis_ready_image"], "n_pairs_strict_core"].sum()),
            "mean_image_median_cell_area_um2": grp["median_cell_area_um2_strict"].mean(),
            "mean_image_median_nuc_area_um2": grp["median_nuc_area_um2_strict"].mean(),
            "mean_image_median_nuc_iod": grp["median_nuc_iod_strict"].mean(),
            "mean_image_median_nc_ratio": grp["median_nc_ratio_strict"].mean(),
            "median_spearman_cell_vs_nuc_area_strict": grp["spearman_cell_vs_nuc_area_strict"].median(),
            "median_spearman_cell_vs_nuc_iod_strict": grp["spearman_cell_vs_nuc_iod_strict"].median(),
            "median_spearman_nuc_area_vs_iod_strict": grp["spearman_nuc_area_vs_iod_strict"].median(),
            "pooled_pearson_cell_vs_nuc_area_strict": _corr_pair(strict_species, "cell_area_um2", "nuc_area_um2", "pearson"),
            "pooled_spearman_cell_vs_nuc_area_strict": _corr_pair(strict_species, "cell_area_um2", "nuc_area_um2", "spearman"),
            "pooled_pearson_cell_vs_nuc_iod_strict": _corr_pair(strict_species, "cell_area_um2", "nuc_iod", "pearson"),
            "pooled_spearman_cell_vs_nuc_iod_strict": _corr_pair(strict_species, "cell_area_um2", "nuc_iod", "spearman"),
        }
        rows.append(row)
    return pd.DataFrame(rows).sort_values("species").reset_index(drop=True)


def global_summary(
    pairs: pd.DataFrame,
    image_summary: pd.DataFrame,
    cell_linkage: pd.DataFrame,
    link_summary: dict,
    image_qc_config: dict,
) -> dict:
    strict = pairs[pairs["keep_strict_core"]].copy()
    trim = pairs[pairs["keep_trim_5_95"]].copy()
    mask_pairs = pairs[pairs["keep_mask_pair"]].copy()
    ready_images = image_summary[image_summary["analysis_ready_image"]].copy()
    ready_files = set(ready_images["filename"])
    ready_strict = strict[strict["filename"].isin(ready_files)].copy()
    return {
        **link_summary,
        **image_qc_config,
        "n_pairs_mask": int(len(mask_pairs)),
        "n_pairs_trim_5_95": int(len(trim)),
        "n_pairs_strict_core": int(len(strict)),
        "n_pairs_ultra_core": int(pairs["keep_ultra_core"].sum()),
        "n_images_with_strict_pairs": int(image_summary["n_pairs_strict_core"].gt(0).sum()),
        "n_species_with_strict_pairs": int(image_summary.loc[image_summary["n_pairs_strict_core"].gt(0), "species"].nunique()),
        "n_analysis_ready_images": int(image_summary["analysis_ready_image"].sum()),
        "n_analysis_ready_pairs": int(len(ready_strict)),
        "median_strict_nc_ratio": float(strict["nc_area_ratio"].median()) if len(strict) else np.nan,
        "pearson_cell_vs_nuc_area_strict": _corr_pair(strict, "cell_area_um2", "nuc_area_um2", "pearson"),
        "spearman_cell_vs_nuc_area_strict": _corr_pair(strict, "cell_area_um2", "nuc_area_um2", "spearman"),
        "pearson_cell_vs_nuc_iod_strict": _corr_pair(strict, "cell_area_um2", "nuc_iod", "pearson"),
        "spearman_cell_vs_nuc_iod_strict": _corr_pair(strict, "cell_area_um2", "nuc_iod", "spearman"),
        "pearson_nuc_area_vs_iod_strict": _corr_pair(strict, "nuc_area_um2", "nuc_iod", "pearson"),
        "spearman_nuc_area_vs_iod_strict": _corr_pair(strict, "nuc_area_um2", "nuc_iod", "spearman"),
        "pearson_cell_vs_nuc_area_analysis_ready": _corr_pair(ready_strict, "cell_area_um2", "nuc_area_um2", "pearson"),
        "spearman_cell_vs_nuc_area_analysis_ready": _corr_pair(ready_strict, "cell_area_um2", "nuc_area_um2", "spearman"),
        "pearson_cell_vs_nuc_iod_analysis_ready": _corr_pair(ready_strict, "cell_area_um2", "nuc_iod", "pearson"),
        "spearman_cell_vs_nuc_iod_analysis_ready": _corr_pair(ready_strict, "cell_area_um2", "nuc_iod", "spearman"),
        "median_image_pearson_cell_vs_nuc_area_analysis_ready": float(ready_images["pearson_cell_vs_nuc_area_strict"].median()) if len(ready_images) else np.nan,
        "median_image_spearman_cell_vs_nuc_area_analysis_ready": float(ready_images["spearman_cell_vs_nuc_area_strict"].median()) if len(ready_images) else np.nan,
        "weighted_image_pearson_cell_vs_nuc_area_analysis_ready": fisher_mean_correlation(
            ready_images["pearson_cell_vs_nuc_area_strict"],
            ready_images["n_pairs_strict_core"],
        ),
        "weighted_image_spearman_cell_vs_nuc_area_analysis_ready": fisher_mean_correlation(
            ready_images["spearman_cell_vs_nuc_area_strict"],
            ready_images["n_pairs_strict_core"],
        ),
        "median_image_pearson_cell_vs_nuc_iod_analysis_ready": float(ready_images["pearson_cell_vs_nuc_iod_strict"].median()) if len(ready_images) else np.nan,
        "median_image_spearman_cell_vs_nuc_iod_analysis_ready": float(ready_images["spearman_cell_vs_nuc_iod_strict"].median()) if len(ready_images) else np.nan,
        "weighted_image_pearson_cell_vs_nuc_iod_analysis_ready": fisher_mean_correlation(
            ready_images["pearson_cell_vs_nuc_iod_strict"],
            ready_images["n_pairs_strict_core"],
        ),
        "weighted_image_spearman_cell_vs_nuc_iod_analysis_ready": fisher_mean_correlation(
            ready_images["spearman_cell_vs_nuc_iod_strict"],
            ready_images["n_pairs_strict_core"],
        ),
        "n_cells_one_to_one": int((cell_linkage["cell_linkage_status"] == "one_to_one").sum()),
    }


def apply_plot_style() -> None:
    plt.rcParams.update(
        {
            "figure.facecolor": "#f6efe8",
            "axes.facecolor": "#fffaf4",
            "axes.edgecolor": "#9a897a",
            "axes.labelcolor": "#2d261f",
            "xtick.color": "#2d261f",
            "ytick.color": "#2d261f",
            "font.family": "DejaVu Serif",
            "savefig.facecolor": "#f6efe8",
        }
    )


def plot_retention(summary: dict, output_path: Path) -> None:
    apply_plot_style()
    labels = ["Cells", "Nuclei", "Mask pairs", "Trimmed", "Strict", "Ultra"]
    values = [
        summary["n_cells_total"],
        summary["n_nuclei_total"],
        summary["n_pairs_mask"],
        summary["n_pairs_trim_5_95"],
        summary["n_pairs_strict_core"],
        summary["n_pairs_ultra_core"],
    ]
    colors = ["#6386a8", "#a66d5b", "#657a4f", "#ab8841", "#8a5d96", "#5d7d74"]
    fig, ax = plt.subplots(figsize=(9, 5))
    ax.bar(labels, values, color=colors, alpha=0.9)
    ax.set_ylabel("Count")
    ax.set_title("Retention across mask-link and trimming tiers")
    ax.grid(axis="y", linestyle=":", linewidth=0.9)
    for x, val in enumerate(values):
        ax.text(x, val + max(values) * 0.015, f"{val:,}", ha="center", va="bottom", fontsize=10)
    fig.tight_layout()
    fig.savefig(output_path, dpi=180, bbox_inches="tight")
    plt.close(fig)


def _sample_for_plot(df: pd.DataFrame, max_points: int = 12000) -> pd.DataFrame:
    if len(df) <= max_points:
        return df
    return df.sample(max_points, random_state=0)


def plot_scatter(df: pd.DataFrame, x: str, y: str, xlabel: str, ylabel: str, title: str, output_path: Path) -> None:
    apply_plot_style()
    plot_df = _sample_for_plot(df)
    fig, ax = plt.subplots(figsize=(7.5, 6))
    ax.scatter(plot_df[x], plot_df[y], s=10, alpha=0.25, color="#4c72b0", edgecolors="none")
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.grid(linestyle=":", linewidth=0.9)
    if len(plot_df) >= 3:
        r = _corr_pair(plot_df, x, y, "pearson")
        rho = _corr_pair(plot_df, x, y, "spearman")
        ax.text(0.04, 0.96, f"Pearson r={r:.3f}\nSpearman rho={rho:.3f}", transform=ax.transAxes, va="top")
    fig.tight_layout()
    fig.savefig(output_path, dpi=180, bbox_inches="tight")
    plt.close(fig)


def plot_ratio_hist(df: pd.DataFrame, output_path: Path) -> None:
    apply_plot_style()
    fig, ax = plt.subplots(figsize=(7.5, 5))
    ax.hist(df["nc_area_ratio"], bins=40, color="#b36f42", alpha=0.85, edgecolor="none")
    med = float(df["nc_area_ratio"].median()) if len(df) else np.nan
    if np.isfinite(med):
        ax.axvline(med, color="#7b1f1f", linestyle="--", linewidth=1.2)
    ax.set_xlabel("Nucleus-to-cell area ratio")
    ax.set_ylabel("Strict-core matched pairs")
    ax.set_title("Distribution of nucleus/cell area ratio")
    ax.grid(axis="y", linestyle=":", linewidth=0.9)
    fig.tight_layout()
    fig.savefig(output_path, dpi=180, bbox_inches="tight")
    plt.close(fig)


def write_methods_summary(output_dir: Path, summary: dict) -> None:
    lines = [
        "# Cell-Nucleus Linkage Methods",
        "",
        "- Matching unit: brightfield nucleus rows linked to saved cell label masks within the same image tile.",
        "- Primary linkage rule: if a saved raster nucleus mask is available, require the nucleus mask to overlap exactly one positive cell label across the full nucleus footprint (`nucleus_mask_overlap_unique`).",
        "- Secondary linkage rule: if a nucleus mask is unavailable or non-unique but an ImageJ ROI artifact exists, require the ROI polygon to overlap exactly one positive cell label (`roi_overlap_unique`).",
        "- Final fallback: if overlap-based linkage is unavailable or non-unique, sample the nucleus centroid inside the saved cell mask (`mask_centroid`) or, if the rounded centroid misses the mask by <= 1 pixel, within a 3x3 window containing exactly one positive cell label (`mask_window_unique`).",
        "- Interpretation note: overlap-based links are preferred because they use the full nucleus footprint rather than a single centroid sample.",
        "- Ambiguous cells with more than one linked nucleus are excluded from paired analyses rather than forced into a single match.",
        "- Physical validity rule: paired cell area must exceed nucleus area.",
        "- Trimming is performed within each image-species group to avoid slide-to-slide staining confounds.",
        "- `trim_5_95`: keep pairs within the 5th-95th percentile for cell area, nucleus area, nucleus IOD, and nucleus/cell area ratio.",
        "- `strict_core`: trimmed pairs with robust MAD-based |z| <= 2.5 for log(cell area), log(nucleus area), log(nucleus IOD), and log(nucleus/cell area ratio), with combined core distance <= 4.5.",
        "- `ultra_core`: same variables with |z| <= 2.0 and combined core distance <= 3.5.",
        "- Image-level advisory QC: mark `analysis_ready_image = True` when an image is in the upper quartile of strict-pair yield per cell and has within-image Spearman cell-vs-nucleus-area correlation >= 0.10.",
        "- Interpretation note: pooled cross-species correlations are descriptive only and can be distorted by species composition. For biological scaling claims, prioritize the within-image summaries across analysis-ready images.",
        "",
        "## Current Snapshot",
        "",
        f"- Cells loaded: `{summary['n_cells_total']}`",
        f"- Nuclei loaded: `{summary['n_nuclei_total']}`",
        f"- Nuclei with saved mask artifacts: `{summary['n_nuclei_with_mask_artifact']}`",
        f"- Nuclei with ROI artifacts: `{summary['n_nuclei_with_roi_artifact']}`",
        f"- Unique nucleus-mask links: `{summary['n_nucleus_mask_unique_hits']}`",
        f"- Unique ROI-overlap links: `{summary['n_roi_overlap_unique_hits']}`",
        f"- Centroid links: `{summary['n_direct_mask_hits']}`",
        f"- Unique-window rescue links: `{summary['n_window_mask_hits']}`",
        f"- Matched nuclei: `{summary['n_matched_nuclei']}`",
        f"- Strict-core matched pairs: `{summary['n_pairs_strict_core']}`",
        f"- Analysis-ready images: `{summary['n_analysis_ready_images']}`",
        f"- Analysis-ready strict pairs: `{summary['n_analysis_ready_pairs']}`",
        f"- Pearson cell vs nucleus area (strict): `{summary['pearson_cell_vs_nuc_area_strict']:.4f}`" if np.isfinite(summary["pearson_cell_vs_nuc_area_strict"]) else "- Pearson cell vs nucleus area (strict): `NA`",
        f"- Pearson cell vs nucleus area (analysis-ready): `{summary['pearson_cell_vs_nuc_area_analysis_ready']:.4f}`" if np.isfinite(summary["pearson_cell_vs_nuc_area_analysis_ready"]) else "- Pearson cell vs nucleus area (analysis-ready): `NA`",
        f"- Median image Spearman cell vs nucleus area (analysis-ready): `{summary['median_image_spearman_cell_vs_nuc_area_analysis_ready']:.4f}`" if np.isfinite(summary["median_image_spearman_cell_vs_nuc_area_analysis_ready"]) else "- Median image Spearman cell vs nucleus area (analysis-ready): `NA`",
        f"- Weighted image Spearman cell vs nucleus area (analysis-ready): `{summary['weighted_image_spearman_cell_vs_nuc_area_analysis_ready']:.4f}`" if np.isfinite(summary["weighted_image_spearman_cell_vs_nuc_area_analysis_ready"]) else "- Weighted image Spearman cell vs nucleus area (analysis-ready): `NA`",
        f"- Pearson cell vs nucleus IOD (strict): `{summary['pearson_cell_vs_nuc_iod_strict']:.4f}`" if np.isfinite(summary["pearson_cell_vs_nuc_iod_strict"]) else "- Pearson cell vs nucleus IOD (strict): `NA`",
        f"- Pearson cell vs nucleus IOD (analysis-ready): `{summary['pearson_cell_vs_nuc_iod_analysis_ready']:.4f}`" if np.isfinite(summary["pearson_cell_vs_nuc_iod_analysis_ready"]) else "- Pearson cell vs nucleus IOD (analysis-ready): `NA`",
        f"- Median image Spearman cell vs nucleus IOD (analysis-ready): `{summary['median_image_spearman_cell_vs_nuc_iod_analysis_ready']:.4f}`" if np.isfinite(summary["median_image_spearman_cell_vs_nuc_iod_analysis_ready"]) else "- Median image Spearman cell vs nucleus IOD (analysis-ready): `NA`",
        f"- Weighted image Spearman cell vs nucleus IOD (analysis-ready): `{summary['weighted_image_spearman_cell_vs_nuc_iod_analysis_ready']:.4f}`" if np.isfinite(summary["weighted_image_spearman_cell_vs_nuc_iod_analysis_ready"]) else "- Weighted image Spearman cell vs nucleus IOD (analysis-ready): `NA`",
        f"- Pearson nucleus area vs IOD (strict): `{summary['pearson_nuc_area_vs_iod_strict']:.4f}`" if np.isfinite(summary["pearson_nuc_area_vs_iod_strict"]) else "- Pearson nucleus area vs IOD (strict): `NA`",
    ]
    (output_dir / "methods_summary.md").write_text("\n".join(lines))


def build_html(output_dir: Path, summary: dict, image_summary: pd.DataFrame, species_summary: pd.DataFrame) -> str:
    def fmt_metric(value: float, digits: int = 3) -> str:
        return f"{value:.{digits}f}" if np.isfinite(value) else "NA"

    figures = [
        ("Retention by tier", "retention_by_tier.png"),
        ("Cell vs nucleus area", "cell_vs_nucleus_area.png"),
        ("Cell vs nucleus IOD", "cell_vs_nucleus_iod.png"),
        ("NC ratio histogram", "nc_ratio_hist.png"),
    ]
    figure_cards = []
    for label, name in figures:
        if not (output_dir / name).exists():
            continue
        figure_cards.append(
            "<div class='card figure'>"
            f"<a href='{html.escape(name)}'><img src='{html.escape(name)}' alt='{html.escape(label)}'></a>"
            f"<p><a href='{html.escape(name)}'>{html.escape(label)}</a></p>"
            "</div>"
        )

    top_species = species_summary.sort_values("n_pairs_strict_core", ascending=False).head(12)
    species_rows = []
    for _, row in top_species.iterrows():
        species_rows.append(
            "<tr>"
            f"<td>{html.escape(str(row['species']))}</td>"
            f"<td>{int(row['n_images'])}</td>"
            f"<td>{int(row['n_analysis_ready_images'])}</td>"
            f"<td>{int(row['n_pairs_strict_core'])}</td>"
            f"<td>{int(row['n_pairs_strict_core_analysis_ready'])}</td>"
            f"<td>{row['mean_image_median_cell_area_um2']:.3f}</td>"
            f"<td>{row['mean_image_median_nuc_area_um2']:.3f}</td>"
            f"<td>{row['mean_image_median_nuc_iod']:.3f}</td>"
            "</tr>"
        )

    image_rows = []
    top_images = image_summary.sort_values("n_pairs_strict_core", ascending=False).head(12)
    for _, row in top_images.iterrows():
        image_rows.append(
            "<tr>"
            f"<td>{html.escape(str(row['filename']))}</td>"
            f"<td>{html.escape(str(row['species']))}</td>"
            f"<td>{html.escape('yes' if bool(row['analysis_ready_image']) else 'no')}</td>"
            f"<td>{int(row['n_pairs_strict_core'])}</td>"
            f"<td>{row['strict_pair_rate_from_cells']:.3f}</td>"
            f"<td>{row['median_cell_area_um2_strict']:.3f}</td>"
            f"<td>{row['median_nuc_area_um2_strict']:.3f}</td>"
            f"<td>{row['median_nuc_iod_strict']:.3f}</td>"
            f"<td>{row['spearman_cell_vs_nuc_area_strict']:.3f}</td>"
            "</tr>"
        )

    return f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>Cell-Nucleus Linkage Report</title>
  <style>
    :root {{
      --bg: #f4ede5;
      --panel: #fffaf5;
      --ink: #2b241c;
      --muted: #6f5c49;
      --line: #d8cab9;
      --accent: #8c5632;
    }}
    body {{
      margin: 0;
      color: var(--ink);
      background: radial-gradient(circle at top left, rgba(140, 86, 50, 0.10), transparent 28%), var(--bg);
      font: 15px/1.5 Georgia, serif;
    }}
    .wrap {{
      max-width: 1480px;
      margin: 0 auto;
      padding: 24px;
    }}
    .hero, .card, table {{
      background: var(--panel);
      border: 1px solid var(--line);
      border-radius: 18px;
      box-shadow: 0 10px 28px rgba(35, 24, 14, 0.06);
    }}
    .hero, .card {{ padding: 18px; }}
    .grid {{
      display: grid;
      grid-template-columns: repeat(auto-fit, minmax(220px, 1fr));
      gap: 14px;
      margin: 18px 0 24px;
    }}
    .metric-label {{
      color: var(--muted);
      text-transform: uppercase;
      letter-spacing: 0.08em;
      font-size: 12px;
    }}
    .metric-value {{
      margin-top: 10px;
      font-size: 28px;
      color: var(--accent);
    }}
    .fig-grid {{
      display: grid;
      grid-template-columns: repeat(auto-fit, minmax(280px, 1fr));
      gap: 14px;
      margin-bottom: 24px;
    }}
    img {{
      width: 100%;
      height: auto;
      display: block;
      border-radius: 10px;
      border: 1px solid #e1d6ca;
      background: #faf4ec;
    }}
    table {{
      width: 100%;
      border-collapse: collapse;
      overflow: hidden;
      margin-bottom: 24px;
    }}
    th, td {{
      padding: 8px 10px;
      border-bottom: 1px solid var(--line);
      text-align: left;
      font-size: 13px;
    }}
    th {{ background: #f0e2d5; }}
    .links a {{
      display: inline-block;
      margin: 0 10px 10px 0;
      padding: 8px 12px;
      border-radius: 999px;
      background: #efe1d2;
      border: 1px solid #dfc9b5;
      color: var(--ink);
      text-decoration: none;
    }}
  </style>
</head>
<body>
  <div class="wrap">
    <section class="hero">
      <h1>Cell-Nucleus Linkage Report</h1>
      <p>Mask-aware pairing of brightfield cell-size estimates to nucleus size and IOD measurements. Ambiguous cells with more than one linked nucleus are excluded rather than force-matched.</p>
      <p>This build links nuclei to cells by full nucleus-mask overlap first, then falls back to ROI overlap, then centroid-inside-mask matching only when overlap evidence is unavailable. Treat overlap-linked pairs as the preferred combined output because they use the whole nucleus footprint instead of a single sampled point.</p>
      <div class="links">
        <a href="linked_nucleus_pairs.csv.gz">linked nucleus pairs</a>
        <a href="cell_linkage_summary.csv.gz">cell linkage summary</a>
        <a href="image_summary.csv">image summary</a>
        <a href="species_summary.csv">species summary</a>
        <a href="summary.json">summary json</a>
        <a href="methods_summary.md">methods summary</a>
      </div>
    </section>

    <section class="grid">
      <div class="card"><div class="metric-label">Cells</div><div class="metric-value">{summary['n_cells_total']:,}</div></div>
      <div class="card"><div class="metric-label">Nuclei</div><div class="metric-value">{summary['n_nuclei_total']:,}</div></div>
      <div class="card"><div class="metric-label">Mask Unique Links</div><div class="metric-value">{summary['n_nucleus_mask_unique_hits']:,}</div></div>
      <div class="card"><div class="metric-label">Matched Nuclei</div><div class="metric-value">{summary['n_matched_nuclei']:,}</div></div>
      <div class="card"><div class="metric-label">Strict Pairs</div><div class="metric-value">{summary['n_pairs_strict_core']:,}</div></div>
      <div class="card"><div class="metric-label">Analysis-Ready Images</div><div class="metric-value">{summary['n_analysis_ready_images']:,}</div></div>
      <div class="card"><div class="metric-label">Median Ready Image rho Area</div><div class="metric-value">{fmt_metric(summary['median_image_spearman_cell_vs_nuc_area_analysis_ready'])}</div></div>
      <div class="card"><div class="metric-label">Median Ready Image rho IOD</div><div class="metric-value">{fmt_metric(summary['median_image_spearman_cell_vs_nuc_iod_analysis_ready'])}</div></div>
    </section>

    <section class="fig-grid">
      {''.join(figure_cards)}
    </section>

    <section>
      <h2>Top Species by Strict-Core Pairs</h2>
      <table>
        <tr>
          <th>Species</th>
          <th>Images</th>
          <th>Ready images</th>
          <th>Strict pairs</th>
          <th>Ready strict pairs</th>
          <th>Mean image median cell area</th>
          <th>Mean image median nucleus area</th>
          <th>Mean image median nucleus IOD</th>
        </tr>
        {''.join(species_rows)}
      </table>
    </section>

    <section>
      <h2>Top Images by Strict-Core Pairs</h2>
      <table>
        <tr>
          <th>Filename</th>
          <th>Species</th>
          <th>Ready</th>
          <th>Strict pairs</th>
          <th>Strict pairs / cells</th>
          <th>Median cell area</th>
          <th>Median nucleus area</th>
          <th>Median nucleus IOD</th>
          <th>Spearman cell vs nucleus area</th>
        </tr>
        {''.join(image_rows)}
      </table>
    </section>
  </div>
</body>
</html>"""


def main() -> None:
    args = parse_args()
    require_exists(args.cell_csv)
    require_exists(args.nucleus_csv)
    args.output_dir.mkdir(parents=True, exist_ok=True)

    cells = normalize_cells(load_measurements(args.cell_csv), args.image_type)
    nuclei = normalize_nuclei(load_measurements(args.nucleus_csv), args.image_type)

    shared_files = sorted(set(cells["filename"]) & set(nuclei["filename"]))
    cells = cells[cells["filename"].isin(shared_files)].copy()
    nuclei = nuclei[nuclei["filename"].isin(shared_files)].copy()

    linked, cell_linkage, link_summary = link_nuclei_to_cells(cells, nuclei, args.window_radius)
    linked = annotate_pair_qc(linked, cell_linkage)
    image_summary = summarize_images(linked, cells, nuclei)
    image_summary, image_qc_config = annotate_image_qc(image_summary)
    species_summary = summarize_species(image_summary, linked)
    summary = global_summary(linked, image_summary, cell_linkage, link_summary, image_qc_config)
    summary["shared_images"] = len(shared_files)
    summary["shared_species"] = int(cells["species"].nunique())
    summary["cell_csv"] = str(args.cell_csv.resolve())
    summary["nucleus_csv"] = str(args.nucleus_csv.resolve())

    linked.to_csv(args.output_dir / "linked_nucleus_pairs.csv.gz", index=False)
    cell_linkage.to_csv(args.output_dir / "cell_linkage_summary.csv.gz", index=False)
    image_summary.to_csv(args.output_dir / "image_summary.csv", index=False, float_format="%.6f")
    species_summary.to_csv(args.output_dir / "species_summary.csv", index=False, float_format="%.6f")
    (args.output_dir / "summary.json").write_text(json.dumps(summary, indent=2))
    write_methods_summary(args.output_dir, summary)

    strict = linked[linked["keep_strict_core"]].copy()
    plot_retention(summary, args.output_dir / "retention_by_tier.png")
    if len(strict):
        plot_scatter(
            strict,
            "cell_area_um2",
            "nuc_area_um2",
            "Cell area (um^2)",
            "Nucleus area (um^2)",
            "Strict-core cell area vs nucleus area",
            args.output_dir / "cell_vs_nucleus_area.png",
        )
        plot_scatter(
            strict,
            "cell_area_um2",
            "nuc_iod",
            "Cell area (um^2)",
            "Nucleus IOD",
            "Strict-core cell area vs nucleus IOD",
            args.output_dir / "cell_vs_nucleus_iod.png",
        )
        plot_ratio_hist(strict, args.output_dir / "nc_ratio_hist.png")

    html_text = build_html(args.output_dir, summary, image_summary, species_summary)
    (args.output_dir / "index.html").write_text(html_text)

    print(f"Shared images: {len(shared_files)}")
    print(f"Linked nuclei with cells: {summary['n_matched_nuclei']}")
    print(f"Strict-core linked pairs: {summary['n_pairs_strict_core']}")
    print(f"Wrote report: {args.output_dir / 'index.html'}")


if __name__ == "__main__":
    main()
