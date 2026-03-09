# Mixed Cellpose + YOLO Workflow

## Current Decision

The current segmentation path forward is:

- **Cells:** custom Cellpose cell model
- **Nuclei:** custom YOLO nucleus model
- **Linkage:** direct nucleus-mask to cell-mask overlap, with centroid fallback only when overlap artifacts are missing

This replaces the previous idea of using YOLO for both objects. On the current smoke test, the mixed stack links substantially more usable cell-nucleus pairs.

## Active Models

- Cellpose cell model:
  - `output/tile_training_round_v1/train/models/desmognathus_tile_round1`
- YOLO nucleus model:
  - `runs/segment/output/yolo_nucleus_training_final_area500_shape/weights/best.pt`

## Why This Is The Chosen Stack

On the same 12 corrected cell training tiles, the mixed stack outperformed the YOLO-cell plus YOLO-nucleus smoke run for pair recovery:

- Mixed stack:
  - `554` matched nuclei
  - `388` strict-core linked pairs
- YOLO + YOLO smoke run:
  - `459` matched nuclei
  - `323` strict-core linked pairs

Mixed-stack linkage output:

- `output/cellpose_cell_yolo_nucleus_linkage_smoke/index.html`
- `output/cellpose_cell_yolo_nucleus_linkage_smoke/summary.json`

Important caveat:

- This is a **smoke validation**, not a held-out benchmark.
- The current comparison is based on `12` corrected tiles from `2` source images:
  - `Process_316_raw_green.ome.tiff`
  - `Process_366_raw_green.ome.tiff`

## Reproducible Smoke Workflow

### 1. Run Cellpose on the corrected cell tiles

```bash
python3 scripts/run_tile_prediction_review.py \
  --tile-bundle output/tile_annotation_bundle_v1 \
  --model-path output/tile_training_round_v1/train/models/desmognathus_tile_round1 \
  --output-dir output/cellpose_cell_eval_smoke \
  --cellpose-python python3 \
  --label-dir output/tile_bootstrap_review_v1/correction_bundle \
  --selection annotated \
  --use-gpu
```

This writes Cellpose `_masks.png` predictions into:

- `output/cellpose_cell_eval_smoke/correction_bundle`

### 2. Convert Cellpose cell masks into linkage-ready measurements

```bash
python3 scripts/measure_cellpose_tile_predictions.py \
  --manifest output/tile_training_round_v1/training_manifest.csv \
  --prediction-dir output/cellpose_cell_eval_smoke/correction_bundle \
  --output-dir output/cellpose_cell_measurements_smoke \
  --min-mask-area 3500 \
  --min-circularity 0.55 \
  --min-solidity 0.92 \
  --max-aspect-ratio 3.0
```

This produces:

- `output/cellpose_cell_measurements_smoke/all_measurements.csv`

The cell filters above are intentionally conservative and based on the corrected cell masks used for training.

### 3. Run the YOLO nucleus model on the same tile manifest

```bash
/home/jake/Projects/cellprofiler_test/.venv-yolo/bin/python scripts/run_yolo_tile_measurements.py \
  --manifest output/tile_training_round_v1/training_manifest.csv \
  --model runs/segment/output/yolo_nucleus_training_final_area500_shape/weights/best.pt \
  --output-dir output/yolo_nucleus_measurements_smoke \
  --object-kind nucleus \
  --background-cache output/nucleus_iod/nucleus_iod_measurements.csv \
  --min-mask-area 500 \
  --min-circularity 0.55 \
  --min-solidity 0.92 \
  --max-aspect-ratio 2.8 \
  --device 0
```

This produces:

- `output/yolo_nucleus_measurements_smoke/all_measurements.csv`

### 4. Link nuclei to cells by mask overlap

```bash
python3 scripts/build_cell_nucleus_linkage_report.py \
  --cell-csv output/cellpose_cell_measurements_smoke/all_measurements.csv \
  --nucleus-csv output/yolo_nucleus_measurements_smoke/all_measurements.csv \
  --output-dir output/cellpose_cell_yolo_nucleus_linkage_smoke
```

This produces:

- `output/cellpose_cell_yolo_nucleus_linkage_smoke/index.html`
- `output/cellpose_cell_yolo_nucleus_linkage_smoke/summary.json`
- `output/cellpose_cell_yolo_nucleus_linkage_smoke/linked_nucleus_pairs.csv.gz`
- `output/cellpose_cell_yolo_nucleus_linkage_smoke/cell_linkage_summary.csv.gz`

## Scripts Added Or Extended For This Path

- `scripts/measure_cellpose_tile_predictions.py`
  - Converts Cellpose `_masks.png` tile predictions into the same measurement contract used by the linker.
- `scripts/run_yolo_tile_measurements.py`
  - Runs the YOLO segmentation model on tile manifests and emits measurement CSVs plus label masks.
- `scripts/build_cell_nucleus_linkage_report.py`
  - Now prefers full nucleus-mask overlap before centroid fallback.

## Operational Notes

- The nucleus model is now shape-filtered and area-filtered before measurement export.
- The cell measurement export uses conservative morphology gates:
  - minimum area `3500 px`
  - minimum circularity `0.55`
  - minimum solidity `0.92`
  - maximum aspect ratio `3.0`
- The nucleus measurement export uses:
  - minimum area `500 px`
  - minimum circularity `0.55`
  - minimum solidity `0.92`
  - maximum aspect ratio `2.8`

## Next Steps

1. Run the mixed stack on a larger, non-training manifest so the linkage quality is measured off the bootstrap tiles.
2. Add a dedicated Cellpose tile measurement runner that predicts and measures in one step, instead of using a prediction bundle plus conversion step.
3. Promote the mixed stack into one wrapper command for routine runs.
4. Recompute image-level QC thresholds on the larger run, because the current "analysis ready" thresholds are based on only two source images.
5. Keep the YOLO cell model as a benchmark and fallback, but treat Cellpose as the primary cell backend unless a larger held-out run reverses the result.
