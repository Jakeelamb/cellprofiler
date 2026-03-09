#!/usr/bin/env python3
"""Run a trained Cellpose model on tile-bundle images and build an HTML review page."""

from __future__ import annotations

import argparse
import csv
import io
import json
import os
import subprocess
import shutil
import sys
from pathlib import Path

import numpy as np
from PIL import Image
from skimage.segmentation import find_boundaries

PROJECT = Path(__file__).resolve().parent.parent
DEFAULT_BUNDLE = PROJECT / "output" / "tile_annotation_bundle_v1"
DEFAULT_MODEL = PROJECT / "output" / "tile_training_round_v1" / "train" / "models" / "desmognathus_tile_round1"
DEFAULT_OUTPUT = PROJECT / "output" / "tile_bootstrap_review_v1"
PREVIEW_SIZE = 512
DEFAULT_CELLPOSE_PYTHON = os.environ.get("CELLPOSE_PYTHON", sys.executable)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tile-bundle", type=Path, default=DEFAULT_BUNDLE)
    parser.add_argument("--model-path", type=Path, default=DEFAULT_MODEL)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--cellpose-python", default=DEFAULT_CELLPOSE_PYTHON)
    parser.add_argument(
        "--label-dir",
        type=Path,
        action="append",
        default=[],
        help="Optional directory containing corrected labels that should count as already annotated",
    )
    parser.add_argument("--selection", choices=["unlabeled", "annotated", "all"], default="unlabeled")
    parser.add_argument("--use-gpu", action="store_true")
    parser.add_argument("--batch-size", type=int, default=4)
    parser.add_argument("--bsize", type=int, default=256)
    parser.add_argument("--flow-threshold", type=float, default=0.4)
    parser.add_argument("--cellprob-threshold", type=float, default=0.0)
    parser.add_argument("--min-size", type=int, default=15)
    parser.add_argument("--max-images", type=int, default=0, help="Optional limit for smoke tests")
    parser.add_argument("--stage-only", action="store_true", help="Only stage selected images into correction_bundle")
    parser.add_argument("--build-only", action="store_true", help="Only build HTML from an existing correction_bundle")
    return parser.parse_args()


def require_exists(path: Path) -> None:
    if not path.exists():
        raise FileNotFoundError(path)


def sanitize(text: str) -> str:
    keep = []
    for ch in str(text):
        keep.append(ch if ch.isalnum() or ch in {"_", "-", "."} else "_")
    return "".join(keep).strip("_") or "item"


def clear_dir(path: Path) -> None:
    if path.exists():
        shutil.rmtree(path)
    path.mkdir(parents=True, exist_ok=True)


def png_bytes(arr: np.ndarray) -> bytes:
    handle = io.BytesIO()
    Image.fromarray(arr).save(handle, format="PNG", compress_level=1)
    return handle.getvalue()


def jpeg_bytes(arr: np.ndarray, quality: int = 88) -> bytes:
    handle = io.BytesIO()
    Image.fromarray(arr).save(handle, format="JPEG", quality=quality, optimize=True)
    return handle.getvalue()


def resize_preview(arr: np.ndarray, max_side: int = PREVIEW_SIZE) -> np.ndarray:
    h, w = arr.shape[:2]
    scale = min(1.0, float(max_side) / float(max(h, w)))
    out_w = max(1, int(round(w * scale)))
    out_h = max(1, int(round(h * scale)))
    image = Image.fromarray(arr)
    if (out_w, out_h) != (w, h):
        image = image.resize((out_w, out_h), Image.Resampling.LANCZOS)
    return np.asarray(image)


def labeled_u16(mask: np.ndarray) -> np.ndarray:
    arr = np.asarray(mask)
    if arr.ndim > 2:
        arr = np.squeeze(arr)
    arr = np.rint(arr).astype(np.int64, copy=False)
    arr[arr < 0] = 0
    if arr.max(initial=0) > np.iinfo(np.uint16).max:
        raise ValueError("mask labels exceed uint16 range")
    return arr.astype(np.uint16, copy=False)


def has_any_label(image_path: Path, label_dirs: list[Path]) -> bool:
    candidates = [image_path]
    for label_dir in label_dirs:
        candidates.append(label_dir / image_path.name)
    for candidate in candidates:
        base = candidate.with_suffix("")
        if (base.parent / f"{base.name}_seg.npy").exists():
            return True
        if (base.parent / f"{base.name}_cp_masks.png").exists():
            return True
        if (base.parent / f"{base.name}_cp_masks.tif").exists():
            return True
    return False


def load_manifest_rows(
    bundle_dir: Path,
    selection: str,
    max_images: int,
    label_dirs: list[Path],
) -> list[dict[str, str]]:
    manifest_path = bundle_dir / "manifest.csv"
    require_exists(manifest_path)
    rows = list(csv.DictReader(manifest_path.open(newline="")))
    out = []
    for row in rows:
        image_path = bundle_dir / row["image_path"]
        is_annotated = has_any_label(image_path, label_dirs)
        if selection == "unlabeled" and is_annotated:
            continue
        if selection == "annotated" and not is_annotated:
            continue
        row = dict(row)
        row["image_path_abs"] = str(image_path.resolve())
        row["is_annotated"] = "1" if is_annotated else "0"
        out.append(row)
    if max_images > 0:
        out = out[:max_images]
    return out


def label_to_rgb(labels: np.ndarray) -> np.ndarray:
    labels = np.asarray(labels)
    rgb = np.zeros(labels.shape + (3,), dtype=np.uint8)
    if labels.max(initial=0) <= 0:
        return rgb
    palette = np.array(
        [
            [239, 83, 80],
            [255, 202, 40],
            [102, 187, 106],
            [66, 165, 245],
            [171, 71, 188],
            [255, 112, 67],
            [38, 198, 218],
            [124, 179, 66],
        ],
        dtype=np.uint8,
    )
    positive = labels > 0
    colors = palette[labels[positive] % len(palette)]
    rgb[positive] = colors
    return rgb


def overlay_labels(raw_u8: np.ndarray, labels: np.ndarray) -> np.ndarray:
    rgb = np.repeat(raw_u8[..., None], 3, axis=2)
    boundaries = find_boundaries(labels, mode="outer")
    rgb[boundaries, 0] = 255
    rgb[boundaries, 1] = 96
    rgb[boundaries, 2] = 48
    return rgb


def mask_metrics(labels: np.ndarray) -> tuple[int, float, int, int]:
    labels = np.asarray(labels)
    n_masks = int(labels.max(initial=0))
    if n_masks == 0:
        return 0, 0.0, 0, 0
    counts = np.bincount(labels.ravel())[1:]
    if len(counts) == 0:
        return 0, 0.0, 0, 0
    coverage = float((labels > 0).mean())
    return n_masks, coverage, int(np.median(counts)), int(np.max(counts))


def build_html(entries: list[dict[str, str | int | float]], output_path: Path) -> None:
    data_json = json.dumps(entries, separators=(",", ":"))
    html_text = f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>Tile Bootstrap Review</title>
  <style>
    * {{ box-sizing: border-box; }}
    body {{
      margin: 0;
      background: #111615;
      color: #e6eeea;
      font: 14px/1.45 Georgia, "Times New Roman", serif;
    }}
    .toolbar {{
      position: sticky;
      top: 0;
      z-index: 20;
      display: flex;
      gap: 10px;
      flex-wrap: wrap;
      align-items: center;
      padding: 12px 18px;
      background: rgba(15, 23, 21, 0.97);
      border-bottom: 1px solid #2e433e;
    }}
    .toolbar select {{
      border: 1px solid #45675d;
      background: #173129;
      color: #e6eeea;
      border-radius: 8px;
      padding: 8px 12px;
      font: inherit;
    }}
    .stats {{ color: #a7bbb4; }}
    .grid {{
      display: grid;
      grid-template-columns: repeat(auto-fit, minmax(380px, 1fr));
      gap: 16px;
      padding: 18px;
    }}
    .card {{
      background: #192623;
      border: 1px solid #355047;
      border-radius: 14px;
      overflow: hidden;
    }}
    .card-meta {{
      padding: 12px 14px 10px;
    }}
    .title {{
      font-size: 15px;
      margin-bottom: 6px;
      overflow-wrap: anywhere;
    }}
    .sub {{
      font-size: 12px;
      color: #a7bbb4;
      margin-bottom: 8px;
      overflow-wrap: anywhere;
    }}
    .metrics {{
      display: grid;
      grid-template-columns: repeat(2, minmax(0, 1fr));
      gap: 4px 10px;
      font-size: 12px;
      color: #d3dfda;
    }}
    .imgs {{
      display: grid;
      grid-template-columns: repeat(3, minmax(0, 1fr));
      gap: 6px;
      padding: 0 12px 12px;
    }}
    .imgs img {{
      width: 100%;
      display: block;
      background: #09100f;
      border-radius: 10px;
      cursor: zoom-in;
    }}
    .path {{
      padding: 0 14px 14px;
      font-size: 12px;
      color: #b6c8c2;
      overflow-wrap: anywhere;
    }}
    .lightbox {{
      position: fixed;
      inset: 0;
      display: none;
      align-items: center;
      justify-content: center;
      padding: 20px;
      background: rgba(0, 0, 0, 0.9);
      z-index: 40;
    }}
    .lightbox.active {{ display: flex; }}
    .lightbox img {{
      max-width: min(96vw, 1800px);
      max-height: 94vh;
      border-radius: 12px;
      background: #000;
    }}
  </style>
</head>
<body>
  <div class="toolbar">
    <select id="fileFilter" onchange="setFileFilter(this.value)"></select>
    <select id="speciesFilter" onchange="setSpeciesFilter(this.value)"></select>
    <span class="stats" id="stats"></span>
  </div>
  <main class="grid" id="content"></main>
  <div class="lightbox" id="lightbox" onclick="closeLightbox()"><img id="lightboxImg" alt="zoom"></div>
  <script>
    const ENTRIES = {data_json};
    let fileFilter = '__all__';
    let speciesFilter = '__all__';

    function setFileFilter(value) {{
      fileFilter = value;
      render();
    }}

    function setSpeciesFilter(value) {{
      speciesFilter = value;
      render();
    }}

    function filteredEntries() {{
      return ENTRIES.filter((entry) => {{
        if (fileFilter !== '__all__' && entry.filename !== fileFilter) return false;
        if (speciesFilter !== '__all__' && entry.species !== speciesFilter) return false;
        return true;
      }});
    }}

    function populateFilters() {{
      const files = ['__all__', ...Array.from(new Set(ENTRIES.map((entry) => entry.filename))).sort()];
      const species = ['__all__', ...Array.from(new Set(ENTRIES.map((entry) => entry.species))).sort()];
      document.getElementById('fileFilter').innerHTML = files.map((value) => {{
        const label = value === '__all__' ? 'All files' : value;
        const selected = value === fileFilter ? ' selected' : '';
        return `<option value="${{value}}"${{selected}}>${{label}}</option>`;
      }}).join('');
      document.getElementById('speciesFilter').innerHTML = species.map((value) => {{
        const label = value === '__all__' ? 'All species' : value;
        const selected = value === speciesFilter ? ' selected' : '';
        return `<option value="${{value}}"${{selected}}>${{label}}</option>`;
      }}).join('');
    }}

    function openLightbox(path) {{
      document.getElementById('lightboxImg').src = path;
      document.getElementById('lightbox').classList.add('active');
    }}

    function closeLightbox() {{
      document.getElementById('lightbox').classList.remove('active');
    }}

    function render() {{
      populateFilters();
      const rows = filteredEntries();
      document.getElementById('stats').textContent = `${{rows.length}} tiles | ${{rows.reduce((sum, row) => sum + row.n_masks, 0)}} predicted cells`;
      document.getElementById('content').innerHTML = rows.map((entry) => `
        <article class="card">
          <div class="card-meta">
            <div class="title">${{entry.filename}} :: ${{entry.tile_name}}</div>
            <div class="sub"><em>${{entry.species}}</em> | slide=${{entry.slide_id}} | decision=${{entry.decision}}</div>
            <div class="metrics">
              <div>pred cells: ${{entry.n_masks}}</div>
              <div>coverage: ${{(entry.coverage * 100).toFixed(2)}}%</div>
              <div>median area: ${{entry.median_area_px}}</div>
              <div>max area: ${{entry.max_area_px}}</div>
              <div>tile y,x: ${{entry.tile_y0}}, ${{entry.tile_x0}}</div>
              <div>source idx: ${{entry.bundle_index}}</div>
            </div>
          </div>
          <div class="imgs">
            <img src="${{entry.preview_raw_path}}" alt="raw" onclick="openLightbox('${{entry.full_raw_path}}')">
            <img src="${{entry.preview_overlay_path}}" alt="overlay" onclick="openLightbox('${{entry.full_overlay_path}}')">
            <img src="${{entry.preview_mask_path}}" alt="mask" onclick="openLightbox('${{entry.full_mask_path}}')">
          </div>
          <div class="path">Correction bundle: <code>${{entry.correction_image_path}}</code></div>
        </article>
      `).join('');
    }}

    document.addEventListener('keydown', (event) => {{
      if (event.key === 'Escape') closeLightbox();
    }});

    render();
  </script>
</body>
</html>
"""
    output_path.write_text(html_text)


def write_readme(output_dir: Path, summary: dict[str, int | float | str]) -> None:
    lines = [
        "# Tile Bootstrap Review",
        "",
        "This directory contains model predictions on tile-bundle images plus a correction-ready bundle.",
        "",
        "## Layout",
        "",
        "- `correction_bundle/`: raw tiles plus predicted `_masks.png` files for Cellpose GUI correction",
        "- `assets/`: browser preview images",
        "- `predictions_manifest.csv`: per-tile prediction metrics and file paths",
        "- `index.html`: static review page",
        "",
        "## Correct Predictions Faster",
        "",
        "1. Open Cellpose GUI with `python3 -m cellpose`.",
        f"2. Load images from `{output_dir / 'correction_bundle'}`.",
        "3. Enable autoload masks or load the matching `_masks.png` file, then correct instead of redrawing from scratch.",
        "",
        "## Summary",
        "",
    ]
    for key, value in sorted(summary.items()):
        lines.append(f"- `{key}`: `{value}`")
    (output_dir / "README.md").write_text("\n".join(lines) + "\n")


def cellpose_major_version(python_exe: str) -> int | None:
    cmd = [
        python_exe,
        "-c",
        "import importlib.metadata as m; print(m.version('cellpose').split('.')[0])",
    ]
    try:
        result = subprocess.run(cmd, check=True, capture_output=True, text=True)
    except Exception:
        return None
    text = result.stdout.strip()
    return int(text) if text.isdigit() else None


def run_cellpose_cli(correction_dir: Path, args: argparse.Namespace) -> None:
    major = cellpose_major_version(args.cellpose_python)
    cmd = [
        args.cellpose_python,
        "-m",
        "cellpose",
        "--dir",
        str(correction_dir),
        "--pretrained_model",
        str(args.model_path),
        "--save_png",
        "--output_name",
        "_masks",
        "--no_npy",
        "--batch_size",
        str(args.batch_size),
        "--flow_threshold",
        str(args.flow_threshold),
        "--cellprob_threshold",
        str(args.cellprob_threshold),
        "--min_size",
        str(args.min_size),
    ]
    if major is None or major >= 4:
        cmd.extend(["--bsize", str(args.bsize)])
    if args.use_gpu:
        cmd.append("--use_gpu")
    print("$", " ".join(cmd))
    subprocess.run(cmd, check=True)


def main() -> None:
    args = parse_args()
    if args.stage_only and args.build_only:
        raise SystemExit("--stage-only and --build-only are mutually exclusive")
    require_exists(args.tile_bundle)
    require_exists(args.model_path)
    label_dirs = [path.resolve() for path in args.label_dir]
    for label_dir in label_dirs:
        require_exists(label_dir)
    rows = load_manifest_rows(args.tile_bundle, args.selection, args.max_images, label_dirs)
    if not rows:
        raise SystemExit("No tiles matched the requested selection")

    correction_dir = args.output_dir / "correction_bundle"
    assets_dir = args.output_dir / "assets"
    if not args.build_only:
        clear_dir(args.output_dir)
        correction_dir.mkdir(parents=True, exist_ok=True)
        assets_dir.mkdir(parents=True, exist_ok=True)
        for row in rows:
            image_path = Path(str(row["image_path_abs"]))
            correction_image = correction_dir / image_path.name
            shutil.copy2(image_path, correction_image)
    else:
        require_exists(correction_dir)
        assets_dir.mkdir(parents=True, exist_ok=True)

    if args.stage_only:
        summary = {
            "tile_bundle": str(args.tile_bundle.resolve()),
            "selection": args.selection,
            "n_tiles": len(rows),
            "staged_dir": str(correction_dir.resolve()),
        }
        (args.output_dir / "summary.json").write_text(json.dumps(summary, indent=2))
        print(json.dumps(summary, indent=2))
        return

    if not args.build_only:
        run_cellpose_cli(correction_dir, args)

    entries: list[dict[str, str | int | float]] = []
    total_masks = 0
    missing_predictions = 0
    for idx, row in enumerate(rows, start=1):
        image_path = Path(str(row["image_path_abs"]))
        bundle_name = image_path.name
        correction_image = correction_dir / bundle_name
        correction_mask = correction_dir / f"{image_path.stem}_masks.png"
        if not correction_mask.exists():
            missing_predictions += 1
            print(f"SKIP missing mask: {correction_mask.name}")
            continue

        raw_u8 = np.asarray(Image.open(correction_image))
        masks_u16 = labeled_u16(np.asarray(Image.open(correction_mask)))
        overlay = overlay_labels(raw_u8, masks_u16)
        mask_rgb = label_to_rgb(masks_u16)

        n_masks, coverage, median_area_px, max_area_px = mask_metrics(masks_u16)
        total_masks += n_masks

        stem = sanitize(bundle_name)
        raw_preview_name = f"{stem}__raw.jpg"
        overlay_preview_name = f"{stem}__overlay.jpg"
        mask_preview_name = f"{stem}__mask.png"
        full_raw_name = f"{stem}__raw_full.png"
        full_overlay_name = f"{stem}__overlay_full.png"
        full_mask_name = f"{stem}__mask_full.png"

        (assets_dir / raw_preview_name).write_bytes(jpeg_bytes(resize_preview(raw_u8)))
        (assets_dir / overlay_preview_name).write_bytes(jpeg_bytes(resize_preview(overlay)))
        (assets_dir / mask_preview_name).write_bytes(png_bytes(resize_preview(mask_rgb)))
        (assets_dir / full_raw_name).write_bytes(png_bytes(raw_u8))
        (assets_dir / full_overlay_name).write_bytes(png_bytes(overlay))
        (assets_dir / full_mask_name).write_bytes(png_bytes(mask_rgb))

        entry = {
            **row,
            "n_masks": n_masks,
            "coverage": round(coverage, 6),
            "median_area_px": median_area_px,
            "max_area_px": max_area_px,
            "prediction_mask_path": str(correction_mask.relative_to(args.output_dir)),
            "correction_image_path": str(correction_image.relative_to(args.output_dir)),
            "preview_raw_path": f"assets/{raw_preview_name}",
            "preview_overlay_path": f"assets/{overlay_preview_name}",
            "preview_mask_path": f"assets/{mask_preview_name}",
            "full_raw_path": f"assets/{full_raw_name}",
            "full_overlay_path": f"assets/{full_overlay_name}",
            "full_mask_path": f"assets/{full_mask_name}",
        }
        entries.append(entry)
        print(
            f"[{idx}/{len(rows)}] {row['filename']} :: {row['tile_name']} -> "
            f"{n_masks} masks, coverage={coverage:.4f}"
        )

    manifest_path = args.output_dir / "predictions_manifest.csv"
    fieldnames = list(entries[0].keys())
    with manifest_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(entries)

    build_html(entries, args.output_dir / "index.html")
    summary = {
        "tile_bundle": str(args.tile_bundle.resolve()),
        "model_path": str(args.model_path.resolve()),
        "cellpose_python": args.cellpose_python,
        "selection": args.selection,
        "n_tiles": len(entries),
        "n_missing_predictions": missing_predictions,
        "n_total_predicted_masks": total_masks,
        "mean_masks_per_tile": round(total_masks / len(entries), 2) if entries else 0.0,
        "requested_gpu": bool(args.use_gpu),
        "inference_engine": "external_cellpose_cli" if args.build_only else "cellpose_cli",
    }
    (args.output_dir / "summary.json").write_text(json.dumps(summary, indent=2))
    write_readme(args.output_dir, summary)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
