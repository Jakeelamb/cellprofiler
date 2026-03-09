#!/usr/bin/env python3
"""Convert manual nucleus tile labels into a patch-based YOLO segmentation dataset."""

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
DEFAULT_OUTPUT = PROJECT / "output" / "yolo_nucleus_dataset_round1"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--label-dir", type=Path, default=DEFAULT_LABEL_DIR)
    parser.add_argument("--extra-label-dir", type=Path, action="append", default=[])
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--class-name", default="nucleus")
    parser.add_argument("--patch-size", type=int, default=1024)
    parser.add_argument("--stride", type=int, default=1024)
    parser.add_argument("--min-instance-area", type=int, default=8)
    parser.add_argument("--max-instance-area", type=int, default=0)
    parser.add_argument("--min-circularity", type=float, default=0.0)
    parser.add_argument("--min-solidity", type=float, default=0.0)
    parser.add_argument("--max-aspect-ratio", type=float, default=0.0)
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


def load_rows(manifest_path: Path) -> list[dict[str, str]]:
    require_exists(manifest_path)
    rows = list(csv.DictReader(manifest_path.open(newline="")))
    if not rows:
        raise SystemExit("manifest has no rows")
    return rows


def label_path_for_image(label_dir: Path, image_name: str) -> Path:
    return label_dir / f"{Path(image_name).stem}_seg.npy"


def resolve_label_path(label_dirs: list[Path], image_name: str) -> Path | None:
    for label_dir in label_dirs:
        candidate = label_path_for_image(label_dir, image_name)
        if candidate.exists():
            return candidate
    return None


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


def row_split(row: dict[str, str]) -> str:
    split = str(row.get("split", "")).strip().lower()
    if split in {"train", "test", "val"}:
        return image_output_split(split)
    image_path = str(row.get("image_path", "")).strip().lower()
    if image_path.startswith("train/"):
        return "train"
    if image_path.startswith("test/") or image_path.startswith("val/"):
        return "val"
    raise KeyError("row missing split and image_path does not encode train/test")


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


def contour_lines(mask: np.ndarray, patch_size: int, min_instance_area: int) -> list[str]:
    lines: list[str] = []
    for label in np.unique(mask):
        if label <= 0:
            continue
        binary = (mask == label).astype(np.uint8)
        if int(binary.sum()) < min_instance_area:
            continue
        contours, _ = cv2.findContours(binary, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
        for contour in contours:
            area = cv2.contourArea(contour)
            if area < float(min_instance_area):
                continue
            points = contour[:, 0, :].astype(np.float32)
            if len(points) < 3:
                continue
            coords: list[str] = ["0"]
            for x, y in points:
                coords.append(f"{x / patch_size:.6f}")
                coords.append(f"{y / patch_size:.6f}")
            lines.append(" ".join(coords))
    return lines


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


def load_mask(
    seg_path: Path,
    min_area: int,
    max_area: int,
    min_circularity: float,
    min_solidity: float,
    max_aspect_ratio: float,
) -> np.ndarray:
    data = np.load(seg_path, allow_pickle=True).item()
    masks = data.get("masks")
    if masks is None:
        raise ValueError(f"{seg_path} missing 'masks'")
    masks = np.asarray(masks)
    if masks.ndim != 2:
        masks = np.squeeze(masks)
    masks = masks.astype(np.uint16, copy=False)
    return filter_mask_instances(
        masks,
        min_area=min_area,
        max_area=max_area,
        min_circularity=min_circularity,
        min_solidity=min_solidity,
        max_aspect_ratio=max_aspect_ratio,
    )


def image_output_split(split: str) -> str:
    return "val" if split == "test" else "train"


def write_dataset_yaml(output_dir: Path, class_name: str) -> None:
    text = "\n".join(
        [
            f"path: {output_dir.resolve()}",
            "train: images/train",
            "val: images/val",
            "",
            "names:",
            f"  0: {class_name}",
            "",
        ]
    )
    (output_dir / "dataset.yaml").write_text(text)


def main() -> None:
    args = parse_args()
    require_exists(args.label_dir)
    label_dirs = [args.label_dir, *args.extra_label_dir]
    for label_dir in label_dirs:
        require_exists(label_dir)
    rows = load_rows(args.manifest)

    args.output_dir.mkdir(parents=True, exist_ok=True)
    clear_dir(args.output_dir)
    for split in ("train", "val"):
        (args.output_dir / "images" / split).mkdir(parents=True, exist_ok=True)
        (args.output_dir / "labels" / split).mkdir(parents=True, exist_ok=True)

    patch_rows: list[dict[str, object]] = []
    labeled_rows = 0
    unlabeled_rows = 0
    total_instances = 0

    for row in rows:
        image_name = row_image_name(args.manifest, row)
        seg_path = resolve_label_path(label_dirs, image_name)
        if seg_path is None:
            unlabeled_rows += 1
            continue

        split = row_split(row)
        raw_path = row_source_image_path(args.manifest, row)
        require_exists(raw_path)
        image = np.asarray(Image.open(raw_path))
        image = ensure_rgb(image)
        masks = load_mask(
            seg_path,
            min_area=args.min_instance_area,
            max_area=args.max_instance_area,
            min_circularity=args.min_circularity,
            min_solidity=args.min_solidity,
            max_aspect_ratio=args.max_aspect_ratio,
        )
        if image.shape[:2] != masks.shape[:2]:
            raise ValueError(f"shape mismatch for {raw_path.name}: {image.shape[:2]} vs {masks.shape[:2]}")

        labeled_rows += 1
        y_starts = iter_starts(image.shape[0], args.patch_size, args.stride)
        x_starts = iter_starts(image.shape[1], args.patch_size, args.stride)
        for y0 in y_starts:
            for x0 in x_starts:
                crop_img = image[y0 : y0 + args.patch_size, x0 : x0 + args.patch_size]
                crop_mask = masks[y0 : y0 + args.patch_size, x0 : x0 + args.patch_size]
                crop_img = pad_to_patch(crop_img, args.patch_size, fill_value=0)
                crop_mask = pad_to_patch(crop_mask, args.patch_size, fill_value=0)

                lines = contour_lines(crop_mask, args.patch_size, args.min_instance_area)
                positive_instances = int(np.unique(crop_mask)[1:].size)
                total_instances += positive_instances

                stem = f"{Path(image_name).stem}__y{y0:05d}_x{x0:05d}"
                image_out = args.output_dir / "images" / split / f"{stem}.png"
                label_out = args.output_dir / "labels" / split / f"{stem}.txt"
                Image.fromarray(crop_img).save(image_out)
                label_out.write_text("\n".join(lines) + ("\n" if lines else ""))

                patch_rows.append(
                    {
                        "split": split,
                        "source_image": image_name,
                        "patch_image": str(image_out.relative_to(args.output_dir)),
                        "patch_label": str(label_out.relative_to(args.output_dir)),
                        "y0": y0,
                        "x0": x0,
                        "n_instances": positive_instances,
                        "n_segments": len(lines),
                    }
                )

    if labeled_rows == 0:
        raise SystemExit("no labeled rows found in label directory")

    manifest_out = args.output_dir / "patch_manifest.csv"
    with manifest_out.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(patch_rows[0].keys()))
        writer.writeheader()
        writer.writerows(patch_rows)

    summary = {
        "manifest": str(args.manifest.resolve()),
        "label_dirs": [str(path.resolve()) for path in label_dirs],
        "class_name": str(args.class_name),
        "n_labeled_source_tiles": labeled_rows,
        "n_unlabeled_source_tiles": unlabeled_rows,
        "n_train_patches": int(sum(1 for row in patch_rows if row["split"] == "train")),
        "n_val_patches": int(sum(1 for row in patch_rows if row["split"] == "val")),
        "n_positive_patches": int(sum(1 for row in patch_rows if int(row["n_instances"]) > 0)),
        "n_empty_patches": int(sum(1 for row in patch_rows if int(row["n_instances"]) == 0)),
        "n_total_instances_across_patches": int(total_instances),
        "patch_size": int(args.patch_size),
        "stride": int(args.stride),
        "min_instance_area": int(args.min_instance_area),
        "max_instance_area": int(args.max_instance_area),
        "min_circularity": float(args.min_circularity),
        "min_solidity": float(args.min_solidity),
        "max_aspect_ratio": float(args.max_aspect_ratio),
    }
    write_dataset_yaml(args.output_dir, args.class_name)
    (args.output_dir / "summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
