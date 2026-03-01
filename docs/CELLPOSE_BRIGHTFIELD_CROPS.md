# Cellpose Brightfield Crops (Validated Workflow)

## Scope

This workflow is for **brightfield crop images only**.

- Do use: `data/tuning/brightfield/*.tif*`, `data/tuning/subset/*brightfield*.tif*`
- Do not use: full-resolution whole-slide images under `data/brightfield/`

## Objective

Segment both nuclei and cell boundaries with a published pretrained Cellpose model, then keep only isolated, non-clumped cells.

## Runtime Environment

Validated with:

- Python: `/home/jake/bin/miniconda3/bin/python`
- Cellpose: v4 (`cpsam` pretrained model)
- CPU mode (`gpu=False`)

The model file is expected at `~/.cellpose/models/cpsam`.

## Command

```bash
/home/jake/bin/miniconda3/bin/python scripts/run_cellpose_brightfield_crops.py \
  --max-images 3 \
  --run-name validated
```

## Segmentation + Filtering Logic

1. Nucleus-scale segmentation (Cellpose `cpsam`, diameter 12, invert brightfield)
2. Cell-scale segmentation (Cellpose `cpsam`, diameter 62, invert brightfield)
3. Keep only cells that pass all of:
   - area in `[700, 20000]` px
   - solidity `>= 0.85`
   - eccentricity `<= 0.92`
   - extent `>= 0.40`
   - exactly one overlapping nucleus
   - nucleus overlap fraction in `[0.03, 0.55]`
4. Reject edge-touching objects and likely overlaps/clumps

## Output Structure

Each run writes:

```text
results/cellpose_brightfield_crops_<timestamp>_<run_name>/
  config.json
  summary_metrics.csv
  object_metrics.csv
  run_report.md
  images/<image_name>/
    raw_crop.png
    nucleus_mask_raw.tiff
    nucleus_mask_kept.tiff
    cell_mask_raw.tiff
    cell_mask_kept.tiff
    overlay_raw_masks.png
    overlay_kept_vs_rejected.png
```

## QC Metrics

`summary_metrics.csv` includes per-image:

- raw object counts (`raw_cell_count`, `raw_nucleus_count`)
- kept/rejected counts (`kept_count`, `rejected_count`)
- `rejected_overlap_count`
- coverage fractions (`raw_cell_mask_coverage`, `kept_cell_mask_coverage`)
- kept area range (`kept_area_min_px`, `kept_area_max_px`, `kept_area_median_px`)

`object_metrics.csv` includes per-cell status and reject reason for traceable QC.

