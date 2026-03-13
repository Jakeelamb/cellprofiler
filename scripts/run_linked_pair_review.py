#!/usr/bin/env python3
"""Serve a local cell-by-cell linked-pair review app for report validation."""

from __future__ import annotations

import argparse
import csv
import json
import math
import threading
from dataclasses import dataclass
from datetime import datetime, timezone
from functools import lru_cache
from http import HTTPStatus
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from io import BytesIO
from pathlib import Path
from typing import Any
from urllib.parse import parse_qs, urlparse

import numpy as np
import pandas as pd
import tifffile
from PIL import Image
from scipy import ndimage as ndi
from skimage.measure import perimeter as binary_perimeter
from skimage.morphology import convex_hull_image
from skimage.segmentation import find_boundaries


PROJECT = Path(__file__).resolve().parent.parent
DEFAULT_PAIRS_CSV = (
    PROJECT
    / "output"
    / "runs"
    / "mixed_cellpose_yolo_full_dataset_v1_bgclean"
    / "linkage"
    / "linked_nucleus_pairs.csv.gz"
)
DEFAULT_OUTPUT_DIR = (
    PROJECT
    / "output"
    / "runs"
    / "mixed_cellpose_yolo_full_dataset_v1_bgclean"
    / "pair_review"
)
DEFAULT_PORT = 8765
SHAPE_QC_THRESHOLDS = {
    "cell": {
        "smoothness_min": 0.90,
        "ellipse_iou_min": 0.78,
        "solidity_min": 0.93,
    },
    "nucleus": {
        "smoothness_min": 0.92,
        "ellipse_iou_min": 0.82,
        "solidity_min": 0.95,
    },
}
REPAIR_PREVIEW_PARAMS = {
    "cell": {
        "sigma": 1.35,
        "ellipse_weight_force": 0.46,
        "ellipse_weight_jagged": 0.78,
        "ellipse_weight_ellipse": 0.88,
        "ellipse_weight_concave": 0.34,
        "hull_weight_concave": 0.22,
    },
    "nucleus": {
        "sigma": 0.95,
        "ellipse_weight_force": 0.58,
        "ellipse_weight_jagged": 0.70,
        "ellipse_weight_ellipse": 0.80,
        "ellipse_weight_concave": 0.25,
        "hull_weight_concave": 0.18,
    },
}
SHAPE_QC_FIELDNAMES = [
    "review_key",
    "cell_shape_area_px",
    "cell_shape_perimeter_px",
    "cell_shape_smoothness",
    "cell_shape_roughness",
    "cell_shape_ellipse_iou",
    "cell_shape_ellipse_deviation",
    "cell_shape_solidity",
    "cell_shape_aspect_ratio",
    "cell_shape_eccentricity",
    "nucleus_shape_area_px",
    "nucleus_shape_perimeter_px",
    "nucleus_shape_smoothness",
    "nucleus_shape_roughness",
    "nucleus_shape_ellipse_iou",
    "nucleus_shape_ellipse_deviation",
    "nucleus_shape_solidity",
    "nucleus_shape_aspect_ratio",
    "nucleus_shape_eccentricity",
    "shape_has_jagged_edge",
    "shape_has_ellipse_mismatch",
    "shape_has_concavity",
    "shape_is_suspect",
    "shape_qc_status",
    "shape_qc_reasons",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pairs-csv", type=Path, default=DEFAULT_PAIRS_CSV)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument(
        "--selection",
        choices=["strict_core", "mask_pair", "matched_all"],
        default="strict_core",
        help="Subset to review. Default is the current report-included strict-core set.",
    )
    parser.add_argument("--species", default="", help="Optional species filter at startup, e.g. 'D. monticola'")
    parser.add_argument("--max-pairs", type=int, default=0, help="Optional cap for smoke tests")
    parser.add_argument("--host", default="127.0.0.1")
    parser.add_argument("--port", type=int, default=DEFAULT_PORT)
    parser.add_argument("--min-half-size", type=int, default=112)
    parser.add_argument("--max-half-size", type=int, default=256)
    parser.add_argument(
        "--crop-scale",
        type=float,
        default=2.75,
        help="Crop half-width is scaled from cell equivalent radius, then clamped.",
    )
    return parser.parse_args()


def require_exists(path: Path) -> None:
    if not path.exists():
        raise FileNotFoundError(path)


def as_float(value: Any) -> float:
    try:
        out = float(value)
    except (TypeError, ValueError):
        return float("nan")
    return out


def as_int(value: Any, default: int = 0) -> int:
    try:
        if pd.isna(value):
            return default
    except Exception:
        pass
    try:
        return int(round(float(value)))
    except (TypeError, ValueError):
        return default


def as_bool(value: Any) -> bool:
    if isinstance(value, (bool, np.bool_)):
        return bool(value)
    text = str(value).strip().lower()
    if text in {"true", "t", "1", "yes", "y"}:
        return True
    if text in {"false", "f", "0", "no", "n", "", "nan", "none"}:
        return False
    return bool(value)


def clean_text(value: Any) -> str:
    if value is None:
        return ""
    try:
        if pd.isna(value):
            return ""
    except Exception:
        pass
    return str(value)


def normalize_repair_mode(value: Any) -> str:
    text = clean_text(value).strip().lower()
    if text in {"candidate_cell", "cell"}:
        return "cell"
    if text in {"candidate_nucleus", "nucleus"}:
        return "nucleus"
    if text in {"candidate_both", "both", "overlay"}:
        return "both"
    return ""


def normalize_species(value: str) -> str:
    text = clean_text(value).strip()
    if text.startswith("Desmognathus "):
        text = "D. " + text.split(" ", 1)[1]
    return text


def review_key_for_row(row: dict[str, Any]) -> str:
    nucleus_object_id = clean_text(row.get("nucleus_object_id"))
    if nucleus_object_id:
        return nucleus_object_id
    filename = clean_text(row.get("filename"))
    tile_name = clean_text(row.get("tile_name"))
    nucleus_label = as_int(row.get("nucleus_label"), default=-1)
    cell_label = as_int(row.get("mask_label_id"), default=-1)
    return f"{filename}::{tile_name}::nucleus:{nucleus_label}::cell:{cell_label}"


def now_iso() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def image_to_bytes(arr: np.ndarray, fmt: str, **save_kwargs: Any) -> bytes:
    handle = BytesIO()
    Image.fromarray(arr).save(handle, format=fmt, **save_kwargs)
    return handle.getvalue()


def safe_div(numerator: float, denominator: float) -> float:
    if not math.isfinite(numerator) or not math.isfinite(denominator) or denominator == 0:
        return float("nan")
    return numerator / denominator


def round_or_nan(value: Any, digits: int) -> float:
    value = as_float(value)
    return round(value, digits) if math.isfinite(value) else float("nan")


def crop_global(arr: np.ndarray, origin_x: int, origin_y: int, center_x: int, center_y: int, half: int) -> np.ndarray:
    left = center_x - half
    top = center_y - half
    right = center_x + half
    bottom = center_y + half
    src_left = max(left - origin_x, 0)
    src_top = max(top - origin_y, 0)
    src_right = min(right - origin_x, arr.shape[1])
    src_bottom = min(bottom - origin_y, arr.shape[0])

    out_shape = (half * 2, half * 2) + (() if arr.ndim == 2 else (arr.shape[2],))
    out = np.zeros(out_shape, dtype=arr.dtype)
    dst_left = src_left - (left - origin_x)
    dst_top = src_top - (top - origin_y)
    dst_right = dst_left + (src_right - src_left)
    dst_bottom = dst_top + (src_bottom - src_top)
    if src_right > src_left and src_bottom > src_top:
        out[dst_top:dst_bottom, dst_left:dst_right] = arr[src_top:src_bottom, src_left:src_right]
    return out


def normalize_raw_crop(arr: np.ndarray) -> np.ndarray:
    arr = np.asarray(arr)
    if arr.ndim == 3 and arr.shape[2] in {3, 4}:
        arr = arr[:, :, :3]
        if arr.dtype == np.uint8:
            return arr
        out = arr.astype(np.float32)
        lo, hi = np.percentile(out, [1, 99])
        hi = max(hi, lo + 1e-6)
        out = np.clip((out - lo) / (hi - lo), 0, 1)
        return np.rint(out * 255).astype(np.uint8)
    if arr.ndim > 2:
        arr = np.squeeze(arr)
    out = arr.astype(np.float32)
    valid = np.isfinite(out)
    if not valid.any():
        return np.zeros(out.shape + (3,), dtype=np.uint8)
    lo, hi = np.percentile(out[valid], [1, 99])
    hi = max(hi, lo + 1e-6)
    out = np.clip((out - lo) / (hi - lo), 0, 1)
    gray = np.rint(out * 255).astype(np.uint8)
    return np.repeat(gray[:, :, None], 3, axis=2)


def overlay_mask(rgb: np.ndarray, mask: np.ndarray, fill_rgb: tuple[int, int, int], edge_rgb: tuple[int, int, int], alpha: float) -> np.ndarray:
    if mask.shape != rgb.shape[:2]:
        raise ValueError("mask/rgb shape mismatch")
    if not np.any(mask):
        return rgb
    out = rgb.astype(np.float32)
    fill = np.array(fill_rgb, dtype=np.float32)
    edge = np.array(edge_rgb, dtype=np.float32)
    positive = mask.astype(bool)
    out[positive] = (1.0 - alpha) * out[positive] + alpha * fill
    boundary = find_boundaries(positive, mode="outer")
    out[boundary] = edge
    return np.clip(out, 0, 255).astype(np.uint8)


def add_crosshair(rgb: np.ndarray, x: int, y: int) -> np.ndarray:
    out = rgb.copy()
    h, w = out.shape[:2]
    if not (0 <= x < w and 0 <= y < h):
        return out
    color = np.array([255, 255, 255], dtype=np.uint8)
    border = np.array([20, 20, 20], dtype=np.uint8)
    for delta in range(-10, 11):
        if 0 <= x + delta < w:
            out[y, x + delta] = color
            if y - 1 >= 0:
                out[y - 1, x + delta] = border
            if y + 1 < h:
                out[y + 1, x + delta] = border
        if 0 <= y + delta < h:
            out[y + delta, x] = color
            if x - 1 >= 0:
                out[y + delta, x - 1] = border
            if x + 1 < w:
                out[y + delta, x + 1] = border
    return out


@lru_cache(maxsize=32)
def load_raw_array(path_text: str) -> np.ndarray:
    path = Path(path_text)
    if path.suffix.lower() in {".png", ".jpg", ".jpeg"}:
        return np.asarray(Image.open(path))
    return np.asarray(tifffile.imread(path))


@lru_cache(maxsize=64)
def load_label_array(path_text: str) -> np.ndarray:
    return np.asarray(tifffile.imread(path_text))


@lru_cache(maxsize=64)
def load_label_slices(path_text: str) -> tuple[slice | tuple[slice, ...] | None, ...]:
    return tuple(ndi.find_objects(load_label_array(path_text)))


def extract_object_binary(path_text: str, label: int) -> np.ndarray:
    arr = load_label_array(path_text)
    if label > 0:
        slices = load_label_slices(path_text)
        if label - 1 < len(slices):
            obj_slice = slices[label - 1]
            if obj_slice is not None:
                return np.asarray(arr[obj_slice] == label, dtype=bool)
        return np.asarray(arr == label, dtype=bool)
    return np.zeros((0, 0), dtype=bool)


def fit_moment_ellipse_params(mask: np.ndarray) -> tuple[float, float, np.ndarray, float, float] | None:
    coords = np.argwhere(mask)
    if coords.shape[0] < 8:
        return None
    cy, cx = coords.mean(axis=0)
    centered = np.stack([coords[:, 1] - cx, coords[:, 0] - cy], axis=1).astype(np.float64)
    cov = np.cov(centered, rowvar=False, bias=True)
    if cov.shape != (2, 2) or not np.isfinite(cov).all():
        return None
    eigvals, eigvecs = np.linalg.eigh(cov)
    order = np.argsort(eigvals)[::-1]
    eigvals = np.maximum(eigvals[order], 1e-6)
    eigvecs = eigvecs[:, order]
    semi_major = max(1.0, 2.0 * math.sqrt(float(eigvals[0])))
    semi_minor = max(1.0, 2.0 * math.sqrt(float(eigvals[1])))
    return cx, cy, eigvecs, semi_major, semi_minor


def render_ellipse_mask(
    shape: tuple[int, int],
    cx: float,
    cy: float,
    eigvecs: np.ndarray,
    semi_major: float,
    semi_minor: float,
) -> np.ndarray:
    yy, xx = np.indices(shape, dtype=np.float32)
    centered_grid = np.stack([xx - cx, yy - cy], axis=-1)
    rotated = centered_grid @ eigvecs
    return ((rotated[..., 0] / semi_major) ** 2 + (rotated[..., 1] / semi_minor) ** 2) <= 1.0


def fit_moment_ellipse(mask: np.ndarray) -> tuple[np.ndarray, float, float] | None:
    params = fit_moment_ellipse_params(mask)
    if params is None:
        return None
    cx, cy, eigvecs, semi_major, semi_minor = params

    ellipse = render_ellipse_mask(mask.shape, cx, cy, eigvecs, semi_major, semi_minor)
    return ellipse, semi_major, semi_minor


def compute_binary_shape_metrics(mask: np.ndarray) -> dict[str, Any]:
    mask = np.asarray(mask, dtype=bool)
    area_px = int(mask.sum())
    if area_px == 0:
        return {
            "area_px": 0.0,
            "perimeter_px": float("nan"),
            "smoothness": float("nan"),
            "roughness": float("nan"),
            "ellipse_iou": float("nan"),
            "ellipse_deviation": float("nan"),
            "solidity": float("nan"),
            "aspect_ratio": float("nan"),
            "eccentricity": float("nan"),
        }

    perimeter_px = float(binary_perimeter(mask, neighborhood=8))
    convex = convex_hull_image(mask)
    convex_area_px = float(convex.sum())
    solidity = safe_div(float(area_px), convex_area_px)

    ellipse_fit = fit_moment_ellipse(mask)
    if ellipse_fit is None:
        ellipse_iou = float("nan")
        ellipse_deviation = float("nan")
        smoothness = float("nan")
        roughness = float("nan")
        aspect_ratio = float("nan")
        eccentricity = float("nan")
    else:
        ellipse_mask, semi_major, semi_minor = ellipse_fit
        intersection = float(np.logical_and(mask, ellipse_mask).sum())
        union = float(np.logical_or(mask, ellipse_mask).sum())
        ellipse_iou = safe_div(intersection, union)
        ellipse_deviation = 1.0 - ellipse_iou if math.isfinite(ellipse_iou) else float("nan")
        ellipse_perimeter_px = float(binary_perimeter(ellipse_mask, neighborhood=8))
        roughness = safe_div(perimeter_px, ellipse_perimeter_px)
        smoothness_ratio = safe_div(ellipse_perimeter_px, perimeter_px)
        smoothness = min(1.0, smoothness_ratio) if math.isfinite(smoothness_ratio) else float("nan")
        aspect_ratio = safe_div(semi_major, semi_minor)
        ratio_sq = safe_div(semi_minor * semi_minor, semi_major * semi_major)
        eccentricity = math.sqrt(max(0.0, 1.0 - ratio_sq)) if math.isfinite(ratio_sq) else float("nan")

    return {
        "area_px": float(area_px),
        "perimeter_px": perimeter_px,
        "smoothness": smoothness,
        "roughness": roughness,
        "ellipse_iou": ellipse_iou,
        "ellipse_deviation": ellipse_deviation,
        "solidity": solidity,
        "aspect_ratio": aspect_ratio,
        "eccentricity": eccentricity,
    }


def prefix_shape_metrics(prefix: str, metrics: dict[str, Any]) -> dict[str, Any]:
    return {f"{prefix}_{key}": value for key, value in metrics.items()}


def largest_component(mask: np.ndarray) -> np.ndarray:
    mask = np.asarray(mask, dtype=bool)
    if not mask.any():
        return mask
    labels, count = ndi.label(mask)
    if count <= 1:
        return mask
    sizes = np.bincount(labels.ravel())
    sizes[0] = 0
    largest = int(np.argmax(sizes))
    return labels == largest


def ellipse_from_moments(mask: np.ndarray) -> np.ndarray | None:
    params = fit_moment_ellipse_params(mask)
    if params is None:
        return None
    cx, cy, eigvecs, semi_major, semi_minor = params
    target_area = float(np.asarray(mask, dtype=bool).sum())
    ellipse_area = math.pi * semi_major * semi_minor
    if ellipse_area <= 0:
        return None
    scale = math.sqrt(max(target_area, 1.0) / ellipse_area)
    return render_ellipse_mask(mask.shape, cx, cy, eigvecs, semi_major * scale, semi_minor * scale)


def signed_distance_field(mask: np.ndarray) -> np.ndarray:
    mask = np.asarray(mask, dtype=bool)
    inside = ndi.distance_transform_edt(mask)
    outside = ndi.distance_transform_edt(~mask)
    return inside - outside


def threshold_field_to_area(field: np.ndarray, target_area: int) -> np.ndarray:
    field = np.asarray(field, dtype=np.float32)
    if target_area <= 0 or field.size == 0:
        return np.zeros(field.shape, dtype=bool)
    flat = field.ravel()
    target_area = max(1, min(int(target_area), flat.size))
    kth = flat.size - target_area
    threshold = np.partition(flat, kth)[kth]
    mask = field >= threshold
    current_area = int(mask.sum())
    if current_area > target_area:
        ties = np.argwhere(field == threshold)
        drop = current_area - target_area
        if drop > 0 and drop <= len(ties):
            mask = mask.copy()
            for y, x in ties[:drop]:
                mask[y, x] = False
    return mask


def smooth_mask_to_target_area(mask: np.ndarray, sigma: float) -> np.ndarray:
    field = signed_distance_field(mask)
    smoothed = ndi.gaussian_filter(field, sigma=sigma)
    return threshold_field_to_area(smoothed, int(np.asarray(mask, dtype=bool).sum()))


def row_has_shape_reason(row: dict[str, Any], prefix: str, suffix: str) -> bool:
    reasons = {item.strip() for item in clean_text(row.get("shape_qc_reasons")).split(",") if item.strip()}
    return f"{prefix}_{suffix}" in reasons


def build_repair_candidate(mask: np.ndarray, prefix: str, row: dict[str, Any], *, force: bool = False) -> np.ndarray:
    mask = np.asarray(mask, dtype=bool)
    if not mask.any():
        return mask

    smoothness = as_float(row.get(f"{prefix}_shape_smoothness"))
    ellipse_iou = as_float(row.get(f"{prefix}_shape_ellipse_iou"))
    solidity = as_float(row.get(f"{prefix}_shape_solidity"))
    thresholds = SHAPE_QC_THRESHOLDS[prefix]

    needs_shape_repair = any(
        (
            math.isfinite(smoothness) and smoothness < thresholds["smoothness_min"],
            math.isfinite(ellipse_iou) and ellipse_iou < thresholds["ellipse_iou_min"],
            math.isfinite(solidity) and solidity < thresholds["solidity_min"],
        )
    )
    if not needs_shape_repair and not force:
        return mask

    params = REPAIR_PREVIEW_PARAMS[prefix]
    target_area = int(mask.sum())
    candidate = smooth_mask_to_target_area(mask, sigma=params["sigma"])
    candidate = ndi.binary_fill_holes(candidate)
    candidate = largest_component(candidate)

    ellipse_mask = ellipse_from_moments(mask)
    if ellipse_mask is not None:
        candidate_field = signed_distance_field(candidate)
        ellipse_field = signed_distance_field(ellipse_mask)
        ellipse_weight = params["ellipse_weight_force"] if force else 0.0
        if row_has_shape_reason(row, prefix, "ellipse"):
            ellipse_weight = max(ellipse_weight, params["ellipse_weight_ellipse"])
        if row_has_shape_reason(row, prefix, "jagged"):
            ellipse_weight = max(ellipse_weight, params["ellipse_weight_jagged"])
        if row_has_shape_reason(row, prefix, "concave"):
            ellipse_weight = max(ellipse_weight, params["ellipse_weight_concave"])

        if ellipse_weight > 0:
            blended_field = (1.0 - ellipse_weight) * candidate_field + ellipse_weight * ellipse_field
        else:
            blended_field = candidate_field

        if row_has_shape_reason(row, prefix, "concave"):
            hull = convex_hull_image(mask)
            hull_field = signed_distance_field(hull)
            hull_weight = params["hull_weight_concave"]
            blended_field = (1.0 - hull_weight) * blended_field + hull_weight * hull_field

        candidate = threshold_field_to_area(blended_field, target_area)
        candidate = smooth_mask_to_target_area(candidate, sigma=max(0.6, params["sigma"] * 0.7))
        candidate = ndi.binary_fill_holes(candidate)
        candidate = largest_component(candidate)

    return candidate


def shape_flag_reasons(prefix: str, metrics: dict[str, Any]) -> list[str]:
    thresholds = SHAPE_QC_THRESHOLDS[prefix]
    reasons: list[str] = []
    smoothness = as_float(metrics.get(f"{prefix}_shape_smoothness"))
    ellipse_iou = as_float(metrics.get(f"{prefix}_shape_ellipse_iou"))
    solidity = as_float(metrics.get(f"{prefix}_shape_solidity"))
    if math.isfinite(smoothness) and smoothness < thresholds["smoothness_min"]:
        reasons.append(f"{prefix}_jagged")
    if math.isfinite(ellipse_iou) and ellipse_iou < thresholds["ellipse_iou_min"]:
        reasons.append(f"{prefix}_ellipse")
    if math.isfinite(solidity) and solidity < thresholds["solidity_min"]:
        reasons.append(f"{prefix}_concave")
    return reasons


def compute_shape_qc_for_row(row: dict[str, Any]) -> dict[str, Any]:
    out = {name: float("nan") for name in SHAPE_QC_FIELDNAMES if name not in {"review_key", "shape_has_jagged_edge", "shape_has_ellipse_mismatch", "shape_has_concavity", "shape_is_suspect", "shape_qc_status", "shape_qc_reasons"}}

    cell_path = clean_text(row.get("cell_mask_path"))
    cell_label = as_int(row.get("mask_label_id"), default=0)
    if cell_path and Path(cell_path).exists():
        cell_binary = extract_object_binary(cell_path, cell_label)
        out.update(prefix_shape_metrics("cell_shape", compute_binary_shape_metrics(cell_binary)))

    nucleus_path = clean_text(row.get("nucleus_mask_path"))
    nucleus_label = as_int(row.get("nucleus_label"), default=0)
    if nucleus_path and Path(nucleus_path).exists():
        nucleus_binary = extract_object_binary(nucleus_path, nucleus_label)
        out.update(prefix_shape_metrics("nucleus_shape", compute_binary_shape_metrics(nucleus_binary)))

    cell_flags = shape_flag_reasons("cell", out)
    nucleus_flags = shape_flag_reasons("nucleus", out)
    reasons = cell_flags + nucleus_flags
    jagged = any(reason.endswith("_jagged") for reason in reasons)
    ellipse = any(reason.endswith("_ellipse") for reason in reasons)
    concave = any(reason.endswith("_concave") for reason in reasons)
    out.update(
        {
            "shape_has_jagged_edge": jagged,
            "shape_has_ellipse_mismatch": ellipse,
            "shape_has_concavity": concave,
            "shape_is_suspect": bool(reasons),
            "shape_qc_status": "suspect" if reasons else "clean",
            "shape_qc_reasons": ", ".join(reasons),
        }
    )
    return out


@dataclass
class ReviewRecord:
    review_index: int
    review_key: str
    row: dict[str, Any]

    @property
    def decision(self) -> str:
        return clean_text(self.row.get("review_decision")) or "unlabeled"


class DecisionStore:
    def __init__(self, path: Path) -> None:
        self.path = path
        self.path.parent.mkdir(parents=True, exist_ok=True)
        self._lock = threading.Lock()
        self._rows: dict[str, dict[str, str]] = {}
        if self.path.exists():
            with self.path.open(newline="") as handle:
                for row in csv.DictReader(handle):
                    key = clean_text(row.get("review_key"))
                    if key:
                        self._rows[key] = {k: clean_text(v) for k, v in row.items()}

    def get(self, key: str) -> dict[str, str] | None:
        return self._rows.get(key)

    def set(self, row: dict[str, str]) -> None:
        with self._lock:
            self._rows[row["review_key"]] = row
            self._flush()

    def _flush(self) -> None:
        temp = self.path.with_suffix(".tmp")
        fieldnames = [
            "review_key",
            "decision",
            "repair_mode",
            "note",
            "updated_at",
            "review_index",
            "species",
            "filename",
            "tile_name",
            "nucleus_label",
            "cell_label",
            "current_report_keep_strict_core",
        ]
        with temp.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fieldnames)
            writer.writeheader()
            for key in sorted(self._rows):
                writer.writerow({name: clean_text(self._rows[key].get(name)) for name in fieldnames})
        temp.replace(self.path)


class PairReviewApp:
    def __init__(self, args: argparse.Namespace) -> None:
        self.args = args
        require_exists(args.pairs_csv)
        args.output_dir.mkdir(parents=True, exist_ok=True)
        self.decisions = DecisionStore(args.output_dir / "decisions.csv")
        self.records = self._load_records()
        self._attach_shape_qc()
        self._write_manifest()
        self._write_readme()

    def _load_records(self) -> list[ReviewRecord]:
        df = pd.read_csv(self.args.pairs_csv, low_memory=False)
        for col in ["keep_strict_core", "keep_mask_pair", "has_cell_match", "physical_pair_ok", "one_to_one_cell", "keep_trim_5_95", "keep_ultra_core"]:
            if col in df.columns:
                df[col] = df[col].map(as_bool)

        if self.args.selection == "strict_core":
            df = df[df["keep_strict_core"].map(as_bool)]
        elif self.args.selection == "mask_pair":
            df = df[df["keep_mask_pair"].map(as_bool)]
        elif self.args.selection == "matched_all":
            df = df[df["has_cell_match"].map(as_bool)]

        if self.args.species:
            want = normalize_species(self.args.species).lower()
            df = df[df["species"].astype(str).str.lower() == want]

        if self.args.max_pairs > 0:
            df = df.head(self.args.max_pairs)

        if df.empty:
            raise SystemExit("No linked pairs matched the requested selection.")

        records: list[ReviewRecord] = []
        for review_index, row in enumerate(df.to_dict(orient="records")):
            key = review_key_for_row(row)
            decision_row = self.decisions.get(key)
            row["review_key"] = key
            row["review_decision"] = decision_row["decision"] if decision_row else ""
            row["review_repair_mode"] = normalize_repair_mode(decision_row.get("repair_mode")) if decision_row else ""
            row["review_note"] = decision_row["note"] if decision_row else ""
            row["review_updated_at"] = decision_row["updated_at"] if decision_row else ""
            records.append(ReviewRecord(review_index=review_index, review_key=key, row=row))
        return records

    def _shape_metrics_path(self) -> Path:
        return self.args.output_dir / "shape_qc_metrics.csv"

    def _load_shape_cache(self) -> dict[str, dict[str, Any]]:
        path = self._shape_metrics_path()
        if not path.exists():
            return {}
        cache: dict[str, dict[str, Any]] = {}
        with path.open(newline="") as handle:
            for row in csv.DictReader(handle):
                key = clean_text(row.get("review_key"))
                if not key:
                    continue
                parsed: dict[str, Any] = {"review_key": key}
                for name in SHAPE_QC_FIELDNAMES:
                    if name == "review_key":
                        continue
                    value = row.get(name)
                    if name.startswith("shape_has_") or name == "shape_is_suspect":
                        parsed[name] = as_bool(value)
                    elif name in {"shape_qc_status", "shape_qc_reasons"}:
                        parsed[name] = clean_text(value)
                    else:
                        parsed[name] = as_float(value)
                cache[key] = parsed
        return cache

    def _write_shape_cache(self, cache: dict[str, dict[str, Any]]) -> None:
        path = self._shape_metrics_path()
        temp = path.with_suffix(".tmp")
        with temp.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=SHAPE_QC_FIELDNAMES)
            writer.writeheader()
            for key in sorted(cache):
                row = cache[key]
                writer.writerow({name: row.get(name, "") for name in SHAPE_QC_FIELDNAMES})
        temp.replace(path)

    def _attach_shape_qc(self) -> None:
        cache = self._load_shape_cache()
        missing: list[ReviewRecord] = []
        for record in self.records:
            cached = cache.get(record.review_key)
            if cached is None:
                missing.append(record)
                continue
            record.row.update(cached)

        if missing:
            total = len(missing)
            print(f"Computing shape QC metrics for {total} linked pairs...")
            for offset, record in enumerate(missing, start=1):
                metrics = compute_shape_qc_for_row(record.row)
                metrics["review_key"] = record.review_key
                record.row.update(metrics)
                cache[record.review_key] = metrics
                if offset == total or offset % 250 == 0:
                    print(f"  shape QC {offset}/{total}")
            self._write_shape_cache(cache)

    def _write_manifest(self) -> None:
        rows = []
        for record in self.records:
            row = self.serialize_record(record)
            rows.append(row)
        manifest_path = self.args.output_dir / "review_manifest.csv"
        with manifest_path.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
            writer.writeheader()
            writer.writerows(rows)

    def _write_readme(self) -> None:
        readme = self.args.output_dir / "README.md"
        lines = [
            "# Linked Pair Review",
            "",
            "This review app is for validating whether linked cell+nucleus pairs should be counted in downstream reports.",
            "",
            "## Decisions",
            "",
            f"- Decisions are saved to `{self.args.output_dir / 'decisions.csv'}`",
            f"- Review manifest: `{self.args.output_dir / 'review_manifest.csv'}`",
            "",
            "## Hotkeys",
            "",
            "- `K`: keep current pair and advance",
            "- `D`: discard current pair and advance",
            "- `M`: mark current pair maybe and advance",
            "- `R`: mark current pair as repair candidate and advance",
            "- `U`: clear decision",
            "- `ArrowRight` / `L`: next pair",
            "- `ArrowLeft` / `H`: previous pair",
            "- `]`: next unlabeled pair",
            "- `[`: previous unlabeled pair",
            "- `J`: next shape suspect pair",
            "- `O`: show original overlay",
            "- `C`: show cell-repair-only overlay",
            "- `N`: show nucleus-repair-only overlay",
            "- `B` or `V`: show both-repaired overlay",
            "",
            "## Launch",
            "",
            f"```bash\npython3 scripts/run_linked_pair_review.py --pairs-csv {self.args.pairs_csv} --selection {self.args.selection}\n```",
            "",
        ]
        readme.write_text("\n".join(lines))

    def serialize_record(self, record: ReviewRecord) -> dict[str, Any]:
        row = record.row
        decision = clean_text(row.get("review_decision")) or "unlabeled"
        repair_mode = normalize_repair_mode(row.get("review_repair_mode"))
        if decision == "repair" and not repair_mode:
            repair_mode = "both"
        return {
            "review_index": record.review_index,
            "review_key": record.review_key,
            "species": normalize_species(row.get("species")),
            "filename": clean_text(row.get("filename")),
            "tile_name": clean_text(row.get("tile_name")),
            "slide_id": as_int(row.get("slide_id"), default=0),
            "specimen_id": as_int(row.get("specimen_id"), default=0),
            "decision": decision,
            "repair_mode": repair_mode,
            "note": clean_text(row.get("review_note")),
            "updated_at": clean_text(row.get("review_updated_at")),
            "keep_strict_core": as_bool(row.get("keep_strict_core")),
            "keep_mask_pair": as_bool(row.get("keep_mask_pair")),
            "keep_trim_5_95": as_bool(row.get("keep_trim_5_95")),
            "keep_ultra_core": as_bool(row.get("keep_ultra_core")),
            "link_method": clean_text(row.get("link_method")),
            "flag_summary": clean_text(row.get("flag_summary")),
            "cell_area_um2": round(as_float(row.get("cell_area_um2")), 4),
            "nuc_area_um2": round(as_float(row.get("nuc_area_um2")), 4),
            "nc_area_ratio": round(as_float(row.get("nc_area_ratio")), 6),
            "nuc_iod": round(as_float(row.get("nuc_iod")), 6),
            "centroid_distance_px": round(as_float(row.get("centroid_distance_px")), 4),
            "distance_over_cell_radius": round(as_float(row.get("distance_over_cell_radius")), 4),
            "pair_core_distance": round(as_float(row.get("pair_core_distance")), 4),
            "shape_qc_status": clean_text(row.get("shape_qc_status")) or "clean",
            "shape_qc_reasons": clean_text(row.get("shape_qc_reasons")),
            "shape_is_suspect": as_bool(row.get("shape_is_suspect")),
            "shape_has_jagged_edge": as_bool(row.get("shape_has_jagged_edge")),
            "shape_has_ellipse_mismatch": as_bool(row.get("shape_has_ellipse_mismatch")),
            "shape_has_concavity": as_bool(row.get("shape_has_concavity")),
            "cell_shape_smoothness": round_or_nan(row.get("cell_shape_smoothness"), 4),
            "cell_shape_ellipse_iou": round_or_nan(row.get("cell_shape_ellipse_iou"), 4),
            "cell_shape_solidity": round_or_nan(row.get("cell_shape_solidity"), 4),
            "cell_shape_aspect_ratio": round_or_nan(row.get("cell_shape_aspect_ratio"), 4),
            "nucleus_shape_smoothness": round_or_nan(row.get("nucleus_shape_smoothness"), 4),
            "nucleus_shape_ellipse_iou": round_or_nan(row.get("nucleus_shape_ellipse_iou"), 4),
            "nucleus_shape_solidity": round_or_nan(row.get("nucleus_shape_solidity"), 4),
            "nucleus_shape_aspect_ratio": round_or_nan(row.get("nucleus_shape_aspect_ratio"), 4),
            "raw_url": f"/api/image/{record.review_index}?mode=raw",
            "overlay_url": f"/api/image/{record.review_index}?mode=overlay",
            "candidate_cell_overlay_url": f"/api/image/{record.review_index}?mode=candidate_cell_overlay",
            "candidate_nucleus_overlay_url": f"/api/image/{record.review_index}?mode=candidate_nucleus_overlay",
            "candidate_both_overlay_url": f"/api/image/{record.review_index}?mode=candidate_both_overlay",
        }

    def records_payload(self) -> list[dict[str, Any]]:
        return [self.serialize_record(record) for record in self.records]

    def update_decision(self, review_index: int, decision: str, note: str = "", repair_mode: str = "") -> dict[str, Any]:
        if decision not in {"keep", "discard", "maybe", "repair", "unlabeled"}:
            raise ValueError(f"Unsupported decision: {decision}")
        record = self.records[review_index]
        row = record.row
        row["review_decision"] = "" if decision == "unlabeled" else decision
        normalized_mode = normalize_repair_mode(repair_mode)
        row["review_repair_mode"] = (normalized_mode or "both") if decision == "repair" else ""
        row["review_note"] = note.strip()
        row["review_updated_at"] = now_iso()
        decision_row = {
            "review_key": record.review_key,
            "decision": clean_text(row["review_decision"]) or "unlabeled",
            "repair_mode": clean_text(row["review_repair_mode"]),
            "note": clean_text(row["review_note"]),
            "updated_at": clean_text(row["review_updated_at"]),
            "review_index": str(review_index),
            "species": normalize_species(row.get("species")),
            "filename": clean_text(row.get("filename")),
            "tile_name": clean_text(row.get("tile_name")),
            "nucleus_label": str(as_int(row.get("nucleus_label"), default=-1)),
            "cell_label": str(as_int(row.get("mask_label_id"), default=-1)),
            "current_report_keep_strict_core": "1" if as_bool(row.get("keep_strict_core")) else "0",
        }
        self.decisions.set(decision_row)
        return self.serialize_record(record)

    def review_summary(self) -> dict[str, Any]:
        counts = {"keep": 0, "discard": 0, "maybe": 0, "repair": 0, "unlabeled": 0}
        suspect = 0
        for record in self.records:
            counts[record.decision] += 1
            if as_bool(record.row.get("shape_is_suspect")):
                suspect += 1
        return {
            "n_records": len(self.records),
            "counts": counts,
            "shape_suspect_count": suspect,
            "selection": self.args.selection,
            "decisions_path": str((self.args.output_dir / "decisions.csv").resolve()),
        }

    def render_pair_image(self, review_index: int, mode: str) -> tuple[bytes, str]:
        record = self.records[review_index]
        row = record.row
        center_x = as_int(row.get("cell_centroid_x"), default=0)
        center_y = as_int(row.get("cell_centroid_y"), default=0)
        if not math.isfinite(center_x) or center_x == 0:
            center_x = as_int(row.get("nuc_centroid_x"), default=0)
        if not math.isfinite(center_y) or center_y == 0:
            center_y = as_int(row.get("nuc_centroid_y"), default=0)

        cell_radius = as_float(row.get("cell_equiv_radius_px"))
        nuc_area_px = max(as_float(row.get("nuc_area_px")), 0.0)
        nuc_radius = math.sqrt(nuc_area_px / math.pi) if nuc_area_px > 0 else 0.0
        half = int(round(max(self.args.min_half_size, cell_radius * self.args.crop_scale, nuc_radius * 6.0)))
        half = max(self.args.min_half_size, min(self.args.max_half_size, half))

        raw_path = Path(clean_text(row.get("nucleus_source_image_path")))
        raw_origin_x = as_int(row.get("tile_x0"), default=0)
        raw_origin_y = as_int(row.get("tile_y0"), default=0)
        if not raw_path.exists():
            raw_path = Path(clean_text(row.get("cell_source_image_path")))
            raw_origin_x = 0
            raw_origin_y = 0
        raw_arr = load_raw_array(str(raw_path))
        raw_crop = crop_global(raw_arr, raw_origin_x, raw_origin_y, center_x, center_y, half)
        raw_rgb = normalize_raw_crop(raw_crop)

        if mode == "raw":
            return image_to_bytes(raw_rgb, "JPEG", quality=90), "image/jpeg"

        cell_mask_path = Path(clean_text(row.get("cell_mask_path")))
        cell_origin_x = as_int(row.get("cell_tile_x0"), default=raw_origin_x)
        cell_origin_y = as_int(row.get("cell_tile_y0"), default=raw_origin_y)
        cell_binary = np.zeros((half * 2, half * 2), dtype=bool)
        if cell_mask_path.exists():
            cell_mask = load_label_array(str(cell_mask_path))
            cell_label = as_int(row.get("mask_label_id"), default=0)
            cell_crop = crop_global(cell_mask, cell_origin_x, cell_origin_y, center_x, center_y, half)
            cell_binary = cell_crop == cell_label if cell_label > 0 else cell_crop > 0

        nucleus_mask_path = Path(clean_text(row.get("nucleus_mask_path")))
        nucleus_origin_x = as_int(row.get("tile_x0"), default=raw_origin_x)
        nucleus_origin_y = as_int(row.get("tile_y0"), default=raw_origin_y)
        nucleus_binary = np.zeros((half * 2, half * 2), dtype=bool)
        if nucleus_mask_path.exists():
            nucleus_mask = load_label_array(str(nucleus_mask_path))
            nucleus_label = as_int(row.get("nucleus_label"), default=0)
            nucleus_crop = crop_global(nucleus_mask, nucleus_origin_x, nucleus_origin_y, center_x, center_y, half)
            nucleus_binary = nucleus_crop == nucleus_label if nucleus_label > 0 else nucleus_crop > 0

        if mode in {"candidate_cell_overlay", "candidate_both_overlay"}:
            cell_binary = build_repair_candidate(cell_binary, "cell", row, force=True)
        if mode in {"candidate_nucleus_overlay", "candidate_both_overlay"}:
            nucleus_binary = build_repair_candidate(nucleus_binary, "nucleus", row, force=True)

        overlay = raw_rgb
        overlay = overlay_mask(overlay, cell_binary, fill_rgb=(55, 168, 255), edge_rgb=(0, 255, 255), alpha=0.24)
        overlay = overlay_mask(overlay, nucleus_binary, fill_rgb=(255, 151, 59), edge_rgb=(255, 196, 0), alpha=0.45)
        return image_to_bytes(overlay, "PNG", compress_level=1), "image/png"


def build_index_html(app: PairReviewApp) -> str:
    payload = json.dumps(
        {
            "records": app.records_payload(),
            "summary": app.review_summary(),
            "selection": app.args.selection,
            "speciesFilter": normalize_species(app.args.species),
        },
        separators=(",", ":"),
    )
    return f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>Linked Pair Review</title>
  <style>
    * {{ box-sizing: border-box; }}
    body {{
      margin: 0;
      background: #0f1514;
      color: #ecf4f1;
      font: 15px/1.45 Georgia, "Times New Roman", serif;
    }}
    .app {{
      min-height: 100vh;
      display: grid;
      grid-template-rows: auto 1fr;
    }}
    .toolbar {{
      position: sticky;
      top: 0;
      z-index: 20;
      display: grid;
      gap: 10px;
      padding: 12px 16px;
      background: rgba(13, 20, 19, 0.97);
      border-bottom: 1px solid #273835;
    }}
    .toolbar-row {{
      display: flex;
      flex-wrap: wrap;
      gap: 10px 14px;
      align-items: center;
    }}
    .title {{
      font-size: 18px;
      font-weight: 600;
    }}
    .subtle {{
      color: #9eb5ae;
      font-size: 13px;
    }}
    select, button {{
      border: 1px solid #43625b;
      background: #17302a;
      color: #ecf4f1;
      border-radius: 10px;
      padding: 8px 12px;
      font: inherit;
    }}
    button {{
      cursor: pointer;
      min-width: 88px;
    }}
    button.keep {{ background: #1d4b2f; border-color: #3b7b55; }}
    button.discard {{ background: #4a211d; border-color: #8b5248; }}
    button.maybe {{ background: #544319; border-color: #90763e; }}
    button.repair {{ background: #214543; border-color: #4f8a86; }}
    button.secondary {{ background: #1c2624; }}
    .status-pill {{
      display: inline-flex;
      align-items: center;
      gap: 8px;
      padding: 7px 11px;
      border-radius: 999px;
      border: 1px solid #456159;
      background: #17211f;
      color: #d7e7e1;
    }}
    .pill-keep {{ background: #153223; border-color: #356f50; }}
    .pill-discard {{ background: #371a17; border-color: #7f4b43; }}
    .pill-maybe {{ background: #3a3117; border-color: #8e7642; }}
    .pill-suspect {{ background: #332215; border-color: #8f6238; color: #f0dbc6; }}
    .pill-clean {{ background: #132922; border-color: #356b5d; color: #d4eee6; }}
    .main {{
      display: grid;
      grid-template-columns: minmax(0, 1.25fr) minmax(320px, 0.75fr);
      gap: 14px;
      padding: 14px;
    }}
    .viewer {{
      display: grid;
      grid-template-columns: repeat(2, minmax(0, 1fr));
      gap: 14px;
    }}
    .panel {{
      background: #16211f;
      border: 1px solid #2e4640;
      border-radius: 16px;
      overflow: hidden;
      min-height: 0;
    }}
    .panel-head {{
      padding: 10px 12px;
      border-bottom: 1px solid #243632;
      color: #b8ccc5;
      font-size: 13px;
      letter-spacing: 0.02em;
      text-transform: uppercase;
    }}
    .panel img {{
      width: 100%;
      height: calc(100vh - 250px);
      object-fit: contain;
      display: block;
      background: #000;
      cursor: zoom-in;
    }}
    .meta {{
      padding: 14px;
      display: grid;
      gap: 12px;
      align-content: start;
    }}
    .meta-grid {{
      display: grid;
      grid-template-columns: repeat(2, minmax(0, 1fr));
      gap: 8px 12px;
      font-size: 13px;
    }}
    .meta-grid div {{
      padding: 8px 10px;
      background: #101918;
      border-radius: 10px;
      border: 1px solid #22312e;
    }}
    .meta-grid strong {{
      display: block;
      color: #91aea5;
      font-size: 11px;
      text-transform: uppercase;
      letter-spacing: 0.03em;
      margin-bottom: 4px;
    }}
    .hotkeys {{
      font-size: 13px;
      color: #bdd0c9;
      border-top: 1px solid #22312e;
      padding-top: 12px;
    }}
    .progress {{
      font-variant-numeric: tabular-nums;
      font-size: 14px;
    }}
    .lightbox {{
      position: fixed;
      inset: 0;
      display: none;
      align-items: center;
      justify-content: center;
      padding: 20px;
      background: rgba(0,0,0,0.92);
      z-index: 40;
    }}
    .lightbox.active {{ display: flex; }}
    .lightbox img {{
      max-width: 96vw;
      max-height: 95vh;
      background: #000;
      border-radius: 12px;
    }}
    @media (max-width: 1120px) {{
      .main {{ grid-template-columns: 1fr; }}
      .panel img {{ height: 42vh; }}
    }}
  </style>
</head>
<body>
  <div class="app">
    <div class="toolbar">
      <div class="toolbar-row">
        <div class="title">Linked Pair Review</div>
        <div class="subtle">Selection: {app.args.selection}</div>
        <div class="subtle" id="savePath"></div>
      </div>
      <div class="toolbar-row">
        <select id="speciesFilter"></select>
        <select id="decisionFilter"></select>
        <select id="shapeFilter"></select>
        <button class="secondary" onclick="jumpPrevUnlabeled()">Prev Unlabeled [</button>
        <button class="secondary" onclick="jumpNextUnlabeled()">Next Unlabeled ]</button>
        <button class="secondary" onclick="jumpNextSuspect()">Next Suspect J</button>
        <div class="status-pill progress" id="progressText"></div>
        <div class="status-pill" id="summaryText"></div>
        <div class="status-pill" id="shapeSummaryText"></div>
      </div>
      <div class="toolbar-row">
        <button class="keep" onclick="setDecision('keep')">Keep K</button>
        <button class="discard" onclick="setDecision('discard')">Discard D</button>
        <button class="maybe" onclick="setDecision('maybe')">Maybe M</button>
        <button class="repair" onclick="setDecision('repair')">Repair R</button>
        <button class="secondary" onclick="setDecision('unlabeled')">Clear U</button>
        <button class="secondary" onclick="setOverlayMode('overlay')">Original O</button>
        <button class="secondary" onclick="setOverlayMode('candidate_cell')">Cell C</button>
        <button class="secondary" onclick="setOverlayMode('candidate_nucleus')">Nucleus N</button>
        <button class="secondary" onclick="setOverlayMode('candidate_both')">Both B</button>
        <button class="secondary" onclick="move(-1)">Prev ←</button>
        <button class="secondary" onclick="move(1)">Next →</button>
        <div class="status-pill" id="decisionPill"></div>
      </div>
    </div>
    <div class="main">
      <div class="viewer">
        <section class="panel">
          <div class="panel-head">Raw Tile Crop</div>
          <img id="rawImg" alt="raw pair crop" onclick="openLightbox(this.src)">
        </section>
        <section class="panel">
          <div class="panel-head" id="overlayPanelHead">Overlay Crop: Cell + Nucleus Masks</div>
          <img id="overlayImg" alt="overlay pair crop" onclick="openLightbox(this.src)">
        </section>
      </div>
      <aside class="panel">
        <div class="panel-head">Pair Metadata</div>
        <div class="meta">
          <div id="headline"></div>
          <div class="meta-grid" id="metaGrid"></div>
          <div class="hotkeys">
            <div><strong>Hotkeys</strong></div>
            <div>`K` keep, `D` discard, `M` maybe, `R` repair, `U` clear</div>
            <div>`←/H` previous, `→/L` next</div>
            <div>`[` previous unlabeled, `]` next unlabeled</div>
            <div>`J` next shape suspect</div>
            <div>`O` original, `C` cell repair, `N` nucleus repair, `B/V` both repair</div>
            <div>`Esc` closes zoom</div>
          </div>
        </div>
      </aside>
    </div>
  </div>
  <div class="lightbox" id="lightbox" onclick="closeLightbox()">
    <img id="lightboxImg" alt="zoomed review crop">
  </div>
  <script>
    const BOOT = {payload};
    const RECORDS = BOOT.records;
    const DECISION_ORDER = ['unlabeled', 'keep', 'discard', 'maybe', 'repair'];
    const SHAPE_FILTERS = [
      ['__all__', 'All shapes'],
      ['suspect_any', 'Suspect only'],
      ['jagged_any', 'Jagged edge'],
      ['ellipse_any', 'Ellipse mismatch'],
      ['concave_any', 'Concave / low solidity'],
      ['clean_only', 'Likely clean only']
    ];
    let speciesFilter = BOOT.speciesFilter || '__all__';
    let decisionFilter = '__all__';
    let shapeFilter = '__all__';
    let overlayMode = 'overlay';
    let visible = [];
    let visiblePos = 0;

    function counts() {{
      const out = {{ keep: 0, discard: 0, maybe: 0, repair: 0, unlabeled: 0, suspect: 0 }};
      for (const row of RECORDS) {{
        out[row.decision || 'unlabeled'] += 1;
        if (row.shape_is_suspect) out.suspect += 1;
      }}
      return out;
    }}

    function current() {{
      if (!visible.length) return null;
      return RECORDS[visible[visiblePos]];
    }}

    function applyFilters() {{
      visible = [];
      for (let i = 0; i < RECORDS.length; i += 1) {{
        const row = RECORDS[i];
        if (speciesFilter !== '__all__' && row.species !== speciesFilter) continue;
        if (decisionFilter !== '__all__' && row.decision !== decisionFilter) continue;
        if (shapeFilter === 'suspect_any' && !row.shape_is_suspect) continue;
        if (shapeFilter === 'jagged_any' && !row.shape_has_jagged_edge) continue;
        if (shapeFilter === 'ellipse_any' && !row.shape_has_ellipse_mismatch) continue;
        if (shapeFilter === 'concave_any' && !row.shape_has_concavity) continue;
        if (shapeFilter === 'clean_only' && row.shape_is_suspect) continue;
        visible.push(i);
      }}
      if (!visible.length) {{
        visiblePos = 0;
      }} else {{
        visiblePos = Math.max(0, Math.min(visiblePos, visible.length - 1));
      }}
    }}

    function populateFilters() {{
      const species = ['__all__', ...Array.from(new Set(RECORDS.map((row) => row.species))).sort()];
      const speciesSelect = document.getElementById('speciesFilter');
      speciesSelect.innerHTML = species.map((value) => {{
        const label = value === '__all__' ? 'All species' : value;
        const selected = value === speciesFilter ? ' selected' : '';
        return `<option value="${{value}}"${{selected}}>${{label}}</option>`;
      }}).join('');
      const decisionSelect = document.getElementById('decisionFilter');
      const decisionOptions = ['__all__', ...DECISION_ORDER];
      decisionSelect.innerHTML = decisionOptions.map((value) => {{
        const label = value === '__all__' ? 'All decisions' : value;
        const selected = value === decisionFilter ? ' selected' : '';
        return `<option value="${{value}}"${{selected}}>${{label}}</option>`;
      }}).join('');
      const shapeSelect = document.getElementById('shapeFilter');
      shapeSelect.innerHTML = SHAPE_FILTERS.map(([value, label]) => {{
        const selected = value === shapeFilter ? ' selected' : '';
        return `<option value="${{value}}"${{selected}}>${{label}}</option>`;
      }}).join('');
    }}

    function render() {{
      populateFilters();
      applyFilters();
      document.getElementById('savePath').textContent = BOOT.summary.decisions_path;
      const summary = counts();
      document.getElementById('summaryText').textContent =
        `keep ${{summary.keep}} | repair ${{summary.repair}} | discard ${{summary.discard}} | maybe ${{summary.maybe}} | unlabeled ${{summary.unlabeled}}`;
      document.getElementById('shapeSummaryText').textContent =
        `shape suspects ${{summary.suspect}} / ${{RECORDS.length}}`;
      if (!visible.length) {{
        document.getElementById('progressText').textContent = '0 / 0';
        document.getElementById('headline').innerHTML = '<h2>No pairs match the current filters.</h2>';
        document.getElementById('metaGrid').innerHTML = '';
        document.getElementById('rawImg').removeAttribute('src');
        document.getElementById('overlayImg').removeAttribute('src');
        document.getElementById('decisionPill').textContent = 'No selection';
        document.getElementById('decisionPill').className = 'status-pill';
        return;
      }}
      const row = current();
      document.getElementById('progressText').textContent = `${{visiblePos + 1}} / ${{visible.length}} shown | global #${{row.review_index + 1}} / ${{RECORDS.length}}`;
      document.getElementById('headline').innerHTML =
        `<h2 style="margin:0 0 6px 0;">${{row.species}}</h2><div class="subtle">${{row.filename}} :: ${{row.tile_name}}</div>`;
      const metrics = [
        ['Decision', row.decision],
        ['Repair target', row.repair_mode || ''],
        ['Slide / specimen', `${{row.slide_id}} / ${{row.specimen_id}}`],
        ['Report strict-core', row.keep_strict_core ? 'yes' : 'no'],
        ['Mask pair', row.keep_mask_pair ? 'yes' : 'no'],
        ['Link method', row.link_method || ''],
        ['Flags', row.flag_summary || ''],
        ['Cell area um²', Number.isFinite(row.cell_area_um2) ? row.cell_area_um2.toFixed(4) : ''],
        ['Nucleus area um²', Number.isFinite(row.nuc_area_um2) ? row.nuc_area_um2.toFixed(4) : ''],
        ['N:C ratio', Number.isFinite(row.nc_area_ratio) ? row.nc_area_ratio.toFixed(6) : ''],
        ['Nucleus IOD', Number.isFinite(row.nuc_iod) ? row.nuc_iod.toFixed(6) : ''],
        ['Shape QC', row.shape_qc_status || ''],
        ['Shape reasons', row.shape_qc_reasons || ''],
        ['Cell smoothness', Number.isFinite(row.cell_shape_smoothness) ? row.cell_shape_smoothness.toFixed(4) : ''],
        ['Cell ellipse IoU', Number.isFinite(row.cell_shape_ellipse_iou) ? row.cell_shape_ellipse_iou.toFixed(4) : ''],
        ['Cell solidity', Number.isFinite(row.cell_shape_solidity) ? row.cell_shape_solidity.toFixed(4) : ''],
        ['Cell aspect ratio', Number.isFinite(row.cell_shape_aspect_ratio) ? row.cell_shape_aspect_ratio.toFixed(4) : ''],
        ['Nucleus smoothness', Number.isFinite(row.nucleus_shape_smoothness) ? row.nucleus_shape_smoothness.toFixed(4) : ''],
        ['Nucleus ellipse IoU', Number.isFinite(row.nucleus_shape_ellipse_iou) ? row.nucleus_shape_ellipse_iou.toFixed(4) : ''],
        ['Nucleus solidity', Number.isFinite(row.nucleus_shape_solidity) ? row.nucleus_shape_solidity.toFixed(4) : ''],
        ['Nucleus aspect ratio', Number.isFinite(row.nucleus_shape_aspect_ratio) ? row.nucleus_shape_aspect_ratio.toFixed(4) : ''],
        ['Centroid distance px', Number.isFinite(row.centroid_distance_px) ? row.centroid_distance_px.toFixed(4) : ''],
        ['Distance / cell radius', Number.isFinite(row.distance_over_cell_radius) ? row.distance_over_cell_radius.toFixed(4) : ''],
        ['Core distance', Number.isFinite(row.pair_core_distance) ? row.pair_core_distance.toFixed(4) : ''],
        ['Updated', row.updated_at || '']
      ];
      document.getElementById('metaGrid').innerHTML = metrics.map(([label, value]) =>
        `<div><strong>${{label}}</strong>${{value || '&nbsp;'}}</div>`
      ).join('');
      const cacheBust = row.updated_at ? `&t=${{encodeURIComponent(row.updated_at)}}` : '';
      document.getElementById('rawImg').src = row.raw_url + cacheBust;
      const overlayMap = {{
        overlay: row.overlay_url,
        candidate_cell: row.candidate_cell_overlay_url,
        candidate_nucleus: row.candidate_nucleus_overlay_url,
        candidate_both: row.candidate_both_overlay_url
      }};
      const overlayTitles = {{
        overlay: 'Overlay Crop: Cell + Nucleus Masks',
        candidate_cell: 'Overlay Crop: Cell Repair Only',
        candidate_nucleus: 'Overlay Crop: Nucleus Repair Only',
        candidate_both: 'Overlay Crop: Cell + Nucleus Repair'
      }};
      const overlaySrc = overlayMap[overlayMode] || row.overlay_url;
      document.getElementById('overlayImg').src = overlaySrc + cacheBust;
      document.getElementById('overlayPanelHead').textContent = overlayTitles[overlayMode] || overlayTitles.overlay;
      const pill = document.getElementById('decisionPill');
      pill.textContent = `Current: ${{row.decision}}`;
      pill.className = 'status-pill';
      if (row.decision === 'keep') pill.classList.add('pill-keep');
      if (row.decision === 'discard') pill.classList.add('pill-discard');
      if (row.decision === 'maybe') pill.classList.add('pill-maybe');
      if (row.decision === 'repair') pill.classList.add('pill-clean');
      if (row.shape_is_suspect) pill.classList.add('pill-suspect');
    }}

    function move(step) {{
      if (!visible.length) return;
      visiblePos = Math.max(0, Math.min(visible.length - 1, visiblePos + step));
      render();
    }}

    function jumpToUnlabeled(direction) {{
      if (!visible.length) return;
      let pos = visiblePos + direction;
      while (pos >= 0 && pos < visible.length) {{
        const row = RECORDS[visible[pos]];
        if (row.decision === 'unlabeled') {{
          visiblePos = pos;
          render();
          return;
        }}
        pos += direction;
      }}
    }}

    function jumpPrevUnlabeled() {{ jumpToUnlabeled(-1); }}
    function jumpNextUnlabeled() {{ jumpToUnlabeled(1); }}

    function setOverlayMode(mode) {{
      overlayMode = mode;
      render();
    }}

    function jumpNextSuspect() {{
      if (!visible.length) return;
      let pos = visiblePos + 1;
      while (pos < visible.length) {{
        const row = RECORDS[visible[pos]];
        if (row.shape_is_suspect) {{
          visiblePos = pos;
          render();
          return;
        }}
        pos += 1;
      }}
    }}

    async function setDecision(decision) {{
      const row = current();
      if (!row) return;
      const repairModeMap = {{
        overlay: 'both',
        candidate_cell: 'cell',
        candidate_nucleus: 'nucleus',
        candidate_both: 'both'
      }};
      const repair_mode = decision === 'repair' ? (repairModeMap[overlayMode] || 'both') : '';
      const response = await fetch('/api/decision', {{
        method: 'POST',
        headers: {{ 'Content-Type': 'application/json' }},
        body: JSON.stringify({{ review_index: row.review_index, decision, repair_mode }})
      }});
      if (!response.ok) {{
        alert(`Failed to save decision: ${{response.status}}`);
        return;
      }}
      const updated = await response.json();
      Object.assign(RECORDS[row.review_index], updated);
      const step = decision === 'unlabeled' ? 0 : 1;
      render();
      if (step > 0 && visible.length && visiblePos < visible.length - 1) {{
        move(1);
      }}
    }}

    function openLightbox(src) {{
      document.getElementById('lightboxImg').src = src;
      document.getElementById('lightbox').classList.add('active');
    }}

    function closeLightbox() {{
      document.getElementById('lightbox').classList.remove('active');
    }}

    document.addEventListener('change', (event) => {{
      if (event.target.id === 'speciesFilter') {{
        speciesFilter = event.target.value;
        visiblePos = 0;
        render();
      }}
      if (event.target.id === 'decisionFilter') {{
        decisionFilter = event.target.value;
        visiblePos = 0;
        render();
      }}
      if (event.target.id === 'shapeFilter') {{
        shapeFilter = event.target.value;
        visiblePos = 0;
        render();
      }}
    }});

    document.addEventListener('keydown', (event) => {{
      if (event.key === 'Escape') {{
        closeLightbox();
        return;
      }}
      if (event.metaKey || event.ctrlKey || event.altKey) return;
      const key = event.key.toLowerCase();
      if (key === 'k') {{ event.preventDefault(); setDecision('keep'); }}
      else if (key === 'd') {{ event.preventDefault(); setDecision('discard'); }}
      else if (key === 'm') {{ event.preventDefault(); setDecision('maybe'); }}
      else if (key === 'r') {{ event.preventDefault(); setDecision('repair'); }}
      else if (key === 'u') {{ event.preventDefault(); setDecision('unlabeled'); }}
      else if (event.key === 'ArrowRight' || key === 'l') {{ event.preventDefault(); move(1); }}
      else if (event.key === 'ArrowLeft' || key === 'h') {{ event.preventDefault(); move(-1); }}
      else if (event.key === ']') {{ event.preventDefault(); jumpNextUnlabeled(); }}
      else if (event.key === '[') {{ event.preventDefault(); jumpPrevUnlabeled(); }}
      else if (key === 'j') {{ event.preventDefault(); jumpNextSuspect(); }}
      else if (key === 'o') {{ event.preventDefault(); setOverlayMode('overlay'); }}
      else if (key === 'c') {{ event.preventDefault(); setOverlayMode('candidate_cell'); }}
      else if (key === 'n') {{ event.preventDefault(); setOverlayMode('candidate_nucleus'); }}
      else if (key === 'b' || key === 'v') {{ event.preventDefault(); setOverlayMode('candidate_both'); }}
    }});

    render();
  </script>
</body>
</html>
"""


class ReviewHandler(BaseHTTPRequestHandler):
    server_version = "LinkedPairReview/1.0"

    @property
    def app(self) -> PairReviewApp:
        return self.server.app  # type: ignore[attr-defined]

    def do_GET(self) -> None:
        parsed = urlparse(self.path)
        if parsed.path == "/":
            body = build_index_html(self.app).encode("utf-8")
            self._send_bytes(body, "text/html; charset=utf-8")
            return
        if parsed.path == "/api/records":
            body = json.dumps({"records": self.app.records_payload(), "summary": self.app.review_summary()}).encode("utf-8")
            self._send_bytes(body, "application/json; charset=utf-8")
            return
        if parsed.path.startswith("/api/image/"):
            try:
                review_index = int(parsed.path.rsplit("/", 1)[1])
            except ValueError:
                self.send_error(HTTPStatus.BAD_REQUEST, "invalid review index")
                return
            params = parse_qs(parsed.query)
            mode = params.get("mode", ["raw"])[0]
            if mode not in {"raw", "overlay", "candidate_cell_overlay", "candidate_nucleus_overlay", "candidate_both_overlay"}:
                self.send_error(HTTPStatus.BAD_REQUEST, "invalid image mode")
                return
            try:
                body, content_type = self.app.render_pair_image(review_index, mode)
            except IndexError:
                self.send_error(HTTPStatus.NOT_FOUND, "review index out of range")
                return
            except Exception as exc:
                self.send_error(HTTPStatus.INTERNAL_SERVER_ERROR, str(exc))
                return
            self._send_bytes(body, content_type)
            return
        if parsed.path == "/decisions.csv":
            path = self.app.args.output_dir / "decisions.csv"
            if not path.exists():
                self.send_error(HTTPStatus.NOT_FOUND, "decisions file not found")
                return
            self._send_bytes(path.read_bytes(), "text/csv; charset=utf-8")
            return
        self.send_error(HTTPStatus.NOT_FOUND, "not found")

    def do_POST(self) -> None:
        parsed = urlparse(self.path)
        if parsed.path != "/api/decision":
            self.send_error(HTTPStatus.NOT_FOUND, "not found")
            return
        length = int(self.headers.get("Content-Length", "0"))
        try:
            payload = json.loads(self.rfile.read(length) or b"{}")
            review_index = int(payload["review_index"])
            decision = clean_text(payload.get("decision")) or "unlabeled"
            note = clean_text(payload.get("note"))
            repair_mode = normalize_repair_mode(payload.get("repair_mode"))
            row = self.app.update_decision(review_index, decision, note=note, repair_mode=repair_mode)
        except Exception as exc:
            self.send_error(HTTPStatus.BAD_REQUEST, str(exc))
            return
        self._send_bytes(json.dumps(row).encode("utf-8"), "application/json; charset=utf-8")

    def _send_bytes(self, body: bytes, content_type: str) -> None:
        self.send_response(HTTPStatus.OK)
        self.send_header("Content-Type", content_type)
        self.send_header("Cache-Control", "no-store")
        self.send_header("Content-Length", str(len(body)))
        self.end_headers()
        self.wfile.write(body)

    def log_message(self, fmt: str, *args: Any) -> None:
        print(f"[review] {self.address_string()} - {fmt % args}")


def main() -> None:
    args = parse_args()
    app = PairReviewApp(args)
    server = ThreadingHTTPServer((args.host, args.port), ReviewHandler)
    server.app = app  # type: ignore[attr-defined]
    url = f"http://{args.host}:{args.port}/"
    summary = {
        "url": url,
        "selection": args.selection,
        "n_records": len(app.records),
        "species_filter": normalize_species(args.species),
        "decisions_path": str((args.output_dir / "decisions.csv").resolve()),
    }
    (args.output_dir / "server_summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))
    print(f"Serving linked-pair review app at {url}")
    server.serve_forever()


if __name__ == "__main__":
    main()
