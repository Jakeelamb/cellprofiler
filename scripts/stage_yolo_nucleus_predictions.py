#!/usr/bin/env python3
"""Stage YOLO nucleus predictions for unlabeled tiles as a Cellpose correction bundle."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import cv2
import numpy as np
from PIL import Image

PROJECT = Path(__file__).resolve().parent.parent
DEFAULT_MANIFEST = PROJECT / "output" / "nucleus_label_manual_round1" / "manifest.csv"
DEFAULT_LABEL_DIR = PROJECT / "output" / "nucleus_label_manual_round1" / "correction_bundle"
DEFAULT_OUTPUT = PROJECT / "output" / "yolo_nucleus_label_round2"
DEFAULT_MODEL = PROJECT / "runs" / "segment" / "output" / "yolo_nucleus_training_round1" / "weights" / "best.pt"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--label-dir", type=Path, default=DEFAULT_LABEL_DIR)
    parser.add_argument("--extra-label-dir", type=Path, action="append", default=[])
    parser.add_argument("--model", type=Path, default=DEFAULT_MODEL)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--patch-size", type=int, default=1024)
    parser.add_argument("--stride", type=int, default=1024)
    parser.add_argument("--imgsz", type=int, default=1024)
    parser.add_argument("--batch", type=int, default=8)
    parser.add_argument("--conf", type=float, default=0.10)
    parser.add_argument("--iou", type=float, default=0.50)
    parser.add_argument("--max-det", type=int, default=512)
    parser.add_argument("--min-mask-area", type=int, default=16)
    parser.add_argument("--max-mask-area", type=int, default=0)
    parser.add_argument("--min-circularity", type=float, default=0.0)
    parser.add_argument("--min-solidity", type=float, default=0.0)
    parser.add_argument("--max-aspect-ratio", type=float, default=0.0)
    parser.add_argument("--device", default=None)
    return parser.parse_args()


def require_exists(path: Path) -> None:
    if not path.exists():
        raise FileNotFoundError(path)


def project_relative(path: Path) -> str:
    path = path.resolve() if path.is_absolute() else path
    try:
        return str(path.relative_to(PROJECT))
    except ValueError:
        return str(path)


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


def choose_device(requested: str | None) -> str:
    import torch

    if requested:
        return requested
    return "0" if torch.cuda.is_available() else "cpu"


def load_rows(manifest_path: Path) -> tuple[list[str], list[dict[str, str]]]:
    require_exists(manifest_path)
    with manifest_path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        rows = list(reader)
        if reader.fieldnames is None or not rows:
            raise SystemExit("manifest has no rows")
        return reader.fieldnames, rows


def first_present(row: dict[str, str], *keys: str) -> str:
    for key in keys:
        value = str(row.get(key, "")).strip()
        if value:
            return value
    return ""


def resolve_manifest_path(manifest_path: Path, value: str) -> Path:
    path = Path(value)
    if path.is_absolute():
        return path
    return (manifest_path.parent / path).resolve()


def row_image_name(manifest_path: Path, row: dict[str, str]) -> str:
    value = first_present(row, "correction_image_path", "annotation_image_path", "image_path_abs", "image_path")
    if not value:
        raise KeyError("row missing image path fields")
    return resolve_manifest_path(manifest_path, value).name


def row_source_image_path(manifest_path: Path, row: dict[str, str]) -> Path:
    value = first_present(
        row,
        "source_round_image_path",
        "annotation_image_path",
        "image_path_abs",
        "correction_image_path",
        "image_path",
    )
    if not value:
        raise KeyError("row missing source image path fields")
    return resolve_manifest_path(manifest_path, value)


def label_path_for_image(label_dir: Path, image_name: str) -> Path:
    return label_dir / f"{Path(image_name).stem}_seg.npy"


def has_label(image_name: str, label_dirs: list[Path]) -> bool:
    return any(label_path_for_image(label_dir, image_name).exists() for label_dir in label_dirs)


def iter_starts(length: int, patch_size: int, stride: int) -> list[int]:
    if length <= patch_size:
        return [0]
    starts = list(range(0, length - patch_size + 1, stride))
    last = length - patch_size
    if starts[-1] != last:
        starts.append(last)
    return starts


def pad_to_patch(arr: np.ndarray, patch_size: int, fill_value: int = 0) -> np.ndarray:
    h, w = arr.shape[:2]
    if arr.ndim == 2:
        out = np.full((patch_size, patch_size), fill_value, dtype=arr.dtype)
        out[:h, :w] = arr
        return out
    out = np.full((patch_size, patch_size, arr.shape[2]), fill_value, dtype=arr.dtype)
    out[:h, :w, :] = arr
    return out


def ensure_rgb(image: np.ndarray) -> np.ndarray:
    if image.ndim == 2:
        return np.repeat(image[..., None], 3, axis=2)
    if image.ndim == 3 and image.shape[2] == 1:
        return np.repeat(image, 3, axis=2)
    return image


def resize_binary_mask(mask: np.ndarray, shape: tuple[int, int]) -> np.ndarray:
    if mask.shape == shape:
        return mask.astype(bool, copy=False)
    resized = cv2.resize(mask.astype(np.uint8), (shape[1], shape[0]), interpolation=cv2.INTER_NEAREST)
    return resized.astype(bool, copy=False)


def shape_metrics(binary: np.ndarray) -> tuple[float, float, float]:
    contours, _ = cv2.findContours(binary.astype(np.uint8), cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_NONE)
    if not contours:
        return 0.0, 0.0, 0.0
    contour = max(contours, key=cv2.contourArea)
    area = float(cv2.contourArea(contour))
    perimeter = float(cv2.arcLength(contour, True))
    hull = cv2.convexHull(contour)
    hull_area = float(cv2.contourArea(hull))
    x, y, w, h = cv2.boundingRect(contour)

    circularity = 0.0 if perimeter <= 0 else float(4 * np.pi * area / (perimeter * perimeter))
    solidity = 0.0 if hull_area <= 0 else float(area / hull_area)
    aspect_ratio = float(max(w, h) / max(1, min(w, h)))
    return circularity, solidity, aspect_ratio


def filter_mask_instances(
    mask: np.ndarray,
    min_area: int,
    max_area: int,
    min_circularity: float,
    min_solidity: float,
    max_aspect_ratio: float,
) -> np.ndarray:
    filtered = np.zeros(mask.shape, dtype=np.uint16)
    next_label = 1
    values, counts = np.unique(mask, return_counts=True)
    for label, area in zip(values, counts, strict=True):
        if label <= 0:
            continue
        area = int(area)
        if area < min_area:
            continue
        if max_area > 0 and area > max_area:
            continue
        binary = (mask == label)
        circularity, solidity, aspect_ratio = shape_metrics(binary)
        if min_circularity > 0 and circularity < min_circularity:
            continue
        if min_solidity > 0 and solidity < min_solidity:
            continue
        if max_aspect_ratio > 0 and aspect_ratio > max_aspect_ratio:
            continue
        filtered[binary] = next_label
        next_label += 1
    return filtered


def tile_predictions(
    model,
    image: np.ndarray,
    patch_size: int,
    stride: int,
    imgsz: int,
    batch: int,
    conf: float,
    iou: float,
    max_det: int,
    device: str,
    min_mask_area: int,
    max_mask_area: int,
    min_circularity: float,
    min_solidity: float,
    max_aspect_ratio: float,
) -> tuple[np.ndarray, int]:
    image = ensure_rgb(image)
    full_mask = np.zeros(image.shape[:2], dtype=np.uint16)
    patches: list[np.ndarray] = []
    metadata: list[tuple[int, int, int, int]] = []

    for y0 in iter_starts(image.shape[0], patch_size, stride):
        for x0 in iter_starts(image.shape[1], patch_size, stride):
            crop = image[y0 : y0 + patch_size, x0 : x0 + patch_size]
            crop_h, crop_w = crop.shape[:2]
            patches.append(pad_to_patch(crop, patch_size, fill_value=0))
            metadata.append((y0, x0, crop_h, crop_w))

    results = model.predict(
        source=patches,
        imgsz=imgsz,
        conf=conf,
        iou=iou,
        max_det=max_det,
        batch=batch,
        device=device,
        retina_masks=True,
        verbose=False,
    )

    next_label = 1
    for result, (y0, x0, crop_h, crop_w) in zip(results, metadata, strict=True):
        if result.masks is None or result.boxes is None or len(result.boxes) == 0:
            continue

        mask_stack = result.masks.data.cpu().numpy()
        scores = result.boxes.conf.cpu().numpy()
        tile_view = full_mask[y0 : y0 + crop_h, x0 : x0 + crop_w]

        for idx in np.argsort(-scores):
            binary = resize_binary_mask(mask_stack[idx] > 0.5, (crop_h, crop_w))
            if int(binary.sum()) < min_mask_area:
                continue
            assign = binary & (tile_view == 0)
            if int(assign.sum()) < min_mask_area:
                continue
            tile_view[assign] = next_label
            next_label += 1

    full_mask = filter_mask_instances(
        full_mask,
        min_area=min_mask_area,
        max_area=max_mask_area,
        min_circularity=min_circularity,
        min_solidity=min_solidity,
        max_aspect_ratio=max_aspect_ratio,
    )
    return full_mask, int(np.unique(full_mask)[1:].size)


def main() -> None:
    args = parse_args()
    require_exists(args.label_dir)
    label_dirs = [args.label_dir, *args.extra_label_dir]
    for label_dir in label_dirs:
        require_exists(label_dir)
    require_exists(args.model)
    fieldnames, rows = load_rows(args.manifest)

    correction_dir = args.output_dir / "correction_bundle"
    args.output_dir.mkdir(parents=True, exist_ok=True)
    clear_dir(args.output_dir)
    correction_dir.mkdir(parents=True, exist_ok=True)

    from ultralytics import YOLO

    device = choose_device(args.device)
    model = YOLO(str(args.model))
    unlabeled_rows = []
    total_instances = 0
    total_positive_pixels = 0

    manifest_rows: list[dict[str, str]] = []
    extra_fields = ["seed_mask_path", "predicted_instances", "predicted_positive_pixels", "prediction_model"]

    for row in rows:
        image_name = row_image_name(args.manifest, row)
        if has_label(image_name, label_dirs):
            continue

        raw_path = row_source_image_path(args.manifest, row)
        require_exists(raw_path)
        image = np.asarray(Image.open(raw_path))
        full_mask, n_instances = tile_predictions(
            model=model,
            image=image,
            patch_size=args.patch_size,
            stride=args.stride,
            imgsz=args.imgsz,
            batch=args.batch,
            conf=args.conf,
            iou=args.iou,
            max_det=args.max_det,
            device=device,
            min_mask_area=args.min_mask_area,
            max_mask_area=args.max_mask_area,
            min_circularity=args.min_circularity,
            min_solidity=args.min_solidity,
            max_aspect_ratio=args.max_aspect_ratio,
        )
        positive_pixels = int((full_mask > 0).sum())
        total_instances += n_instances
        total_positive_pixels += positive_pixels

        correction_image_path = correction_dir / image_name
        seed_mask_path = correction_dir / f"{Path(image_name).stem}_masks.png"
        Image.fromarray(image).save(correction_image_path)
        Image.fromarray(full_mask).save(seed_mask_path)

        manifest_row = dict(row)
        manifest_row["correction_image_path"] = project_relative(correction_image_path)
        manifest_row["seed_mask_path"] = project_relative(seed_mask_path)
        manifest_row["predicted_instances"] = str(n_instances)
        manifest_row["predicted_positive_pixels"] = str(positive_pixels)
        manifest_row["prediction_model"] = str(args.model.resolve())
        manifest_rows.append(manifest_row)
        unlabeled_rows.append(row)

    if not manifest_rows:
        raise SystemExit("no unlabeled rows found to seed")

    manifest_out = args.output_dir / "manifest.csv"
    with manifest_out.open("w", newline="") as handle:
        output_fields = list(dict.fromkeys([*fieldnames, "correction_image_path", *extra_fields]))
        writer = csv.DictWriter(handle, fieldnames=output_fields)
        writer.writeheader()
        writer.writerows(manifest_rows)

    summary = {
        "source_manifest": str(args.manifest.resolve()),
        "label_dir": str(args.label_dir.resolve()),
        "extra_label_dirs": [str(path.resolve()) for path in args.extra_label_dir],
        "model": str(args.model.resolve()),
        "device": device,
        "output_dir": str(args.output_dir.resolve()),
        "correction_dir": str(correction_dir.resolve()),
        "n_seeded_tiles": len(manifest_rows),
        "n_total_predicted_instances": int(total_instances),
        "n_total_positive_pixels": int(total_positive_pixels),
        "patch_size": int(args.patch_size),
        "stride": int(args.stride),
        "imgsz": int(args.imgsz),
        "batch": int(args.batch),
        "conf": float(args.conf),
        "iou": float(args.iou),
        "max_det": int(args.max_det),
        "min_mask_area": int(args.min_mask_area),
        "max_mask_area": int(args.max_mask_area),
        "min_circularity": float(args.min_circularity),
        "min_solidity": float(args.min_solidity),
        "max_aspect_ratio": float(args.max_aspect_ratio),
    }
    (args.output_dir / "summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
