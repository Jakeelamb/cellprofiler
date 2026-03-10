# Cellsize Segmentation Pipeline (Cellpose)

Automated cell-size measurement pipeline for salamander red blood cells. Uses Cellpose for instance segmentation, then filters for isolated cells and measures area + IOD (integrated optical density).

## Directory Structure

```
Cellsize_segmentation_cellpose_pipeline/
  README.md                  # This file
  scripts/
    segment_cells.py         # Core: Cellpose segmentation + shape/isolation filtering + measurement
    run_batch.py             # Orchestrator: convergence-aware GPU batch processing
    tile_selector.py         # Optional: curses TUI for manual tile pre-filtering
  output/
    brightfield/             # Per-image measurement CSVs + masks + tile manifests
    pmount/                  # Per-image measurement CSVs + masks + tile manifests
    all_measurements.csv     # Combined results across all images
    convergence_summary.csv  # Per-species n, mean, SD, SEM, SEM%
    run_manifest.json        # Full provenance: params, git hash, timestamp
    progress.json            # Resume checkpoint
```

## Data Flow

```
Source image (grayscale green channel, OME-TIFF)
  |
  v
Tiling (4096x4096, auto-skip empty tiles)
  |
  v
Pre-scan: score tiles by sparseness (moderate std = isolated cells)
  |
  v
Process tiles in best-first order (sparsest first):
  1. Cellpose instance segmentation
  2. Shape filter: solidity > 0.7, circularity > 0.4
  3. Isolation filter: nearest neighbor centroid > 50px apart
  4. Edge filter: discard objects within 200px of tile boundary
  |
  v
Measure per cell: area (px + um^2), IOD, mean OD, centroid, i_bg
  |
  v
Per-image CSV + label masks + tile manifest + combined all_measurements.csv
  |
  v
Convergence tracking: stop species when SEM% < threshold
```

## Parameters

### Segmentation (segment_cells.py)

| Parameter | Value | Rationale |
|---|---|---|
| `PIXEL_SIZE_UM` | 0.12 | Olympus APX100-HCU, 40x, confirmed from OME XML |
| `TILE_SIZE` | 4096 | Fits in GPU VRAM (~200MB per tile) |
| `EDGE_MARGIN` | 200 px | = MAX_NUCLEUS_DIAMETER_PX; discard partial objects at tile edges |
| `MIN_AREA_PX` | ~353 | pi * (15px)^2 * 0.5; reject tiny debris |
| `MIN_SOLIDITY` | 0.7 | Reject irregular/fragmented objects |
| `MIN_CIRCULARITY` | 0.4 | 4*pi*area/perimeter^2; reject elongated objects |
| `NEIGHBOR_DISTANCE_PX` | 50 | ~6um; reject touching/clustered cells |
| `PMOUNT_DIAMETER_PX` | 13 | Tuned via d5-d15 sweep for nucleus detection |

Brightfield mode uses auto-diameter (Cellpose decides). Pmount mode uses fixed diameter=13 for nuclei.

### Batch Processing (run_batch.py)

| Flag | Default | Description |
|---|---|---|
| `--gpu` | off | Use CUDA for Cellpose inference |
| `--resume` | off | Skip images already in progress.json |
| `--image-type` | all | `brightfield`, `pmount`, or `all` |
| `--tile-filter` | auto | `green`=high-content, `auto`=skip empties, `yellow`=borderline, `all`=none |
| `--no-crops` | off | Skip saving cell crop images |
| `--max-cells-image` | 500 | Cap cells extracted per image |
| `--min-cells` | 30 | Min cells per species before convergence check |
| `--max-cells-species` | 500 | Hard cap on cells per species |
| `--sem-threshold` | 3.0 | SEM% below which a species is converged |
| `--no-convergence` | off | Disable convergence-based early stopping |
| `--rerun-unconverged` | off | Resume + exempt capped-but-unconverged species from the cap |
| `--output-dir` | `<pipeline>/output` | Output directory |

## Traceability

Every cell measurement row traces back to its source:

1. **Cell -> Tile**: `tile_y0`, `tile_x0`, `tile_h`, `tile_w` define the exact tile region in the source image
2. **Tile -> Image**: `filename` identifies the source OME-TIFF
3. **Image -> Specimen**: `slide_id`, `specimen_id`, `species` from master_image_metadata.csv
4. **Cell location**: `centroid_y`, `centroid_x` in full-image coordinates
5. **Mask provenance**: `mask_path`, `raw_mask_path`, `mask_label_id`, and `tile_manifest_path` point back to the saved label image and tile manifest
6. **Tile priority**: `tile_score` (0-1) indicates sparseness ranking

The `run_manifest.json` records git hash, all parameters, timing, and error log for full reproducibility.

Canonical manifest-driven wrapper:

```bash
uv run python cell_size_segmentation_pipeline/run_from_manifest.py \
    --manifest data/manifests/fail_fast_panel_v1.csv \
    --output-dir output/runs/fail_fast_panel_v1_canonical/cell_size_segmentation
```

## Convergence Methodology

The batch orchestrator tracks per-species measurement stability:

- **Metric**: SEM% = 100 * (SD / sqrt(n)) / mean of cell area
- **Threshold**: Default 3.0% — additional cells won't meaningfully change the species mean
- **Minimum**: At least 30 cells before convergence is checked
- **Maximum**: Hard cap at 500 cells per species (unless exempted)
- **Interleaving**: Jobs are round-robined across species so all get sampled before any is deeply processed
- **Dynamic budgeting**: Per-image cell cap adjusts based on how many cells a species still needs

With typical RBC CV ~15%, convergence at SEM% < 3% requires ~25 cells per species.

### Handling unconverged species

Species that hit the 500-cell cap without reaching the SEM% threshold are flagged as `capped_not_converged` in `convergence_summary.csv`. To continue collecting for these species:

```bash
python Cellsize_segmentation_cellpose_pipeline/scripts/run_batch.py \
    --gpu --no-crops --rerun-unconverged
```

This reloads prior measurements into the tracker, identifies species that were capped but never truly converged, and exempts them from the cap. Processing resumes only for images not yet completed, and exempt species collect cells until SEM% is met (or all images are exhausted).

On resume (`--resume` or `--rerun-unconverged`), the tracker is seeded from the existing `all_measurements.csv` so convergence state is continuous across runs.

## Usage

```bash
# Single image
python Cellsize_segmentation_cellpose_pipeline/scripts/segment_cells.py \
    data/brightfield/Process_312_raw_green.ome.tiff --gpu

# Full batch run (GPU, no crops for speed)
python Cellsize_segmentation_cellpose_pipeline/scripts/run_batch.py \
    --gpu --no-crops --image-type brightfield

# Resume interrupted run
python Cellsize_segmentation_cellpose_pipeline/scripts/run_batch.py \
    --gpu --no-crops --resume

# Continue collecting for species that hit the cap without converging
python Cellsize_segmentation_cellpose_pipeline/scripts/run_batch.py \
    --gpu --no-crops --rerun-unconverged

# Pre-filter tiles interactively (optional)
python Cellsize_segmentation_cellpose_pipeline/scripts/tile_selector.py

# Pmount nuclei (fixed diameter=13)
python Cellsize_segmentation_cellpose_pipeline/scripts/run_batch.py \
    --gpu --no-crops --image-type pmount
```

## Output CSV Columns

| Column | Description |
|---|---|
| `label` | Object ID within the image |
| `area_px` | Cell area in pixels |
| `area_um2` | Cell area in um^2 (area_px * 0.0144) |
| `solidity` | Convex hull fill ratio (0-1) |
| `circularity` | 4*pi*area/perimeter^2 (0-1) |
| `iod` | Integrated optical density (sum of log10(I_bg/I_pixel)) |
| `mean_od` | Mean optical density per pixel |
| `centroid_y` | Y coordinate in full image (pixels) |
| `centroid_x` | X coordinate in full image (pixels) |
| `i_bg` | Background intensity (95th percentile of non-object pixels) |
| `tile_y0`, `tile_x0` | Tile origin in full image |
| `tile_h`, `tile_w` | Tile dimensions |
| `tile_score` | Sparseness score (0-1) |
| `mask_path`, `raw_mask_path` | Saved filtered/raw label mask TIFFs |
| `mask_label_id` | Object label within the saved mask |
| `tile_manifest_path` | Tile-manifest CSV for the image |
| `overlay_path` | Optional debug overlay image |
| `filename` | Source image filename |
| `slide_id` | Slide identifier |
| `specimen_id` | Specimen identifier |
| `species` | Species name |
| `image_type` | `brightfield` or `pmount` |

## Known Limitations

- Deliberately excludes clustered cells — only isolated cells pass the neighbor distance filter
- Tile edge objects are discarded (200px margin), so some valid cells near tile boundaries are lost
- Cellpose auto-diameter can be inconsistent across tiles with very different cell densities
- IOD comparisons across preservation types (fixed vs dried blood) are NOT valid due to differential Feulgen staining
- 2 brightfield images missing from data: Process_342, Process_344 (skipped automatically)
