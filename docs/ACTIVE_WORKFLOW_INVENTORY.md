# Active Workflow Inventory

This is the retained working set for the salamander brightfield pipeline that produced the current full-dataset mixed results.

## Active Data Inputs

- `data/brightfield/`
- `data/manifests/full_dataset_v1.csv`
- `data/metadata/`
- `master_image_metadata.csv`
- `microscope_to_imageid_desmognathus.csv`
- `slide_glass_mapping.csv`
- `slide_species_mapping.csv`

## Active Code

### Cell segmentation

- `cell_size_segmentation_pipeline/run_from_manifest.py`
- `Cellsize_segmentation_cellpose_pipeline/scripts/segment_cells.py`
- `scripts/prepare_cellpose_training_round.py`
- `scripts/prepare_tile_training_round.py`
- `scripts/run_tile_prediction_review.py`
- `scripts/train_cellpose_with_bsize.py`

### Nucleus training and inference

- `scripts/prepare_nucleus_tile_training_round.py`
- `scripts/build_yolo_nucleus_dataset.py`
- `scripts/stage_yolo_nucleus_predictions.py`
- `scripts/run_yolo_nucleus_training.py`
- `scripts/run_cellpose3_gui.sh`

### Nucleus IOD

- `nucleus_iod_estimate_pipeline/run_from_manifest.py`
- `scripts/imagej_nucleus_iod.py`
- `scripts/nucleus_iod_python.py`
- `scripts/nucleus_iod.ijm`

### Mixed linkage

- `scripts/prepare_mixed_linkage_tiles.py`
- `scripts/run_yolo_tile_measurements.py`
- `scripts/build_cell_nucleus_linkage_report.py`

### Shared support

- `src/cellprofiler_tools/pipeline_runs.py`
- `src/cellprofiler_tools/convergence.py`

## Active Model and Training Artifacts

### Cell model

- `output/tile_annotation_bundle_v1`
- `output/tile_bootstrap_review_v1`
- `output/tile_training_round_v1`

### Nucleus model

- `output/nucleus_label_manual_round1`
- `output/yolo_nucleus_label_round2`
- `output/yolo_nucleus_label_round3`
- `output/yolo_nucleus_dataset_final_area500_shape`
- `output/yolo_models/yolo26n-seg.pt`
- `runs/segment/output/yolo_nucleus_training_final_area500_shape`

## Active Full-Run Outputs

- `output/runs/full_dataset_v1`
- `output/runs/mixed_cellpose_yolo_full_dataset_v1`

## Cleanup Rule

Anything moved into `legacy_20260310/`, `output/archive_20260310_active_prune/`, or `runs/segment/archive_20260310_active_prune/` is not part of the current working pipeline and should be treated as historical only.
