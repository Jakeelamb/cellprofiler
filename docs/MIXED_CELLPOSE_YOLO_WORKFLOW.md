# Mixed Cellpose + YOLO Workflow

## Current Stack

- Cells: custom Cellpose cell model
- Nuclei: custom YOLO nucleus model
- Nucleus intensity: brightfield nucleus IOD background cache
- Linkage: full nucleus-mask overlap to cell-mask overlap, with centroid fallback

## Active Models

- Cellpose cell model:
  - `output/tile_training_round_v1/train/models/desmognathus_tile_round1`
- YOLO nucleus model:
  - `runs/segment/output/yolo_nucleus_training_final_area500_shape/weights/best.pt`

## Retained Training Assets

- Cell training bundle:
  - `output/tile_annotation_bundle_v1`
  - `output/tile_bootstrap_review_v1`
  - `output/tile_training_round_v1`
- Nucleus training bundle:
  - `output/nucleus_label_manual_round1`
  - `output/yolo_nucleus_label_round2`
  - `output/yolo_nucleus_label_round3`
  - `output/yolo_nucleus_dataset_final_area500_shape`

## Full-Dataset Outputs

- Reused full Cellpose cell run:
  - `output/runs/full_dataset_v1`
- Final mixed run:
  - `output/runs/mixed_cellpose_yolo_full_dataset_v1`

Key final artifacts:

- `output/runs/mixed_cellpose_yolo_full_dataset_v1/linkage/index.html`
- `output/runs/mixed_cellpose_yolo_full_dataset_v1/linkage/summary.json`
- `output/runs/mixed_cellpose_yolo_full_dataset_v1/linkage/linked_nucleus_pairs.csv.gz`
- `output/runs/mixed_cellpose_yolo_full_dataset_v1/linkage/cell_linkage_summary.csv.gz`

## Reproducing The Current Mixed Run On New Brightfield Images

### 1. Run Cellpose cell segmentation

```bash
python3 cell_size_segmentation_pipeline/run_from_manifest.py \
  --manifest data/manifests/full_dataset_v1.csv \
  --output-dir output/runs/<run_tag>/cell_size_segmentation \
  --backend cellpose \
  --cellpose-model output/tile_training_round_v1/train/models/desmognathus_tile_round1 \
  --gpu \
  --image-type brightfield \
  --resume
```

### 2. Build a brightfield nucleus background cache

```bash
python3 nucleus_iod_estimate_pipeline/run_from_manifest.py \
  --manifest data/manifests/full_dataset_v1.csv \
  --output-dir output/runs/<run_tag>/nucleus_iod \
  --image-type brightfield \
  --backend imagej \
  --resume
```

This produces the background cache used for tile-level nucleus IOD measurement:

- `output/runs/<run_tag>/nucleus_iod/measurements/nucleus_iod_measurements.csv`

### 3. Extract the exact cell-linked tiles

```bash
python3 scripts/prepare_mixed_linkage_tiles.py \
  --cell-csv output/runs/<run_tag>/cell_size_segmentation/all_measurements.csv \
  --output-dir output/runs/<run_tag>/prepare
```

### 4. Run the YOLO nucleus model on those tiles

```bash
/home/jake/Projects/cellprofiler_test/.venv-yolo/bin/python scripts/run_yolo_tile_measurements.py \
  --manifest output/runs/<run_tag>/prepare/tile_manifest.csv \
  --model runs/segment/output/yolo_nucleus_training_final_area500_shape/weights/best.pt \
  --output-dir output/runs/<run_tag>/nucleus_measurements \
  --object-kind nucleus \
  --background-cache output/runs/<run_tag>/nucleus_iod/measurements/nucleus_iod_measurements.csv \
  --min-mask-area 500 \
  --min-circularity 0.55 \
  --min-solidity 0.92 \
  --max-aspect-ratio 2.8 \
  --device 0
```

### 5. Link nuclei to cells

```bash
python3 scripts/build_cell_nucleus_linkage_report.py \
  --cell-csv output/runs/<run_tag>/prepare/cell_measurements_backfilled.csv \
  --nucleus-csv output/runs/<run_tag>/nucleus_measurements/all_measurements.csv \
  --output-dir output/runs/<run_tag>/linkage
```

## Traceability Contract

The final linked table keeps:

- source image paths for both objects
- both run manifests
- both tile manifests
- both mask paths
- tile coordinates
- object IDs
- overlap/linkage statistics

## Archive Rule

Historical scripts, docs, experiments, and older outputs were moved to:

- `legacy_20260310/`
- `output/archive_20260310_active_prune/`
- `runs/segment/archive_20260310_active_prune/`

Those are preserved for reference only and are not part of the active workflow.
