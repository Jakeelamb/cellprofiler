# Cell Size Segmentation Pipeline

## Scope

This folder is the canonical documentation anchor for the Cellpose-based cell-size segmentation pipeline.

Repo-level shared architecture:

- `/home/jake/Projects/cellprofiler_test/three_pipelines/README.md`

It does not move or duplicate the working code or outputs. It defines:

- the canonical data sources
- the scripts that own the segmentation pass
- the expected output locations
- the traceability target
- the current gaps that still need to be closed

## Current Status

- The canonical wrapper now supports two brightfield backends:
  - `rules`: fast nucleus-first watershed segmentation with geometry filters
  - `cellpose`: the older model-based pass
- The fast `rules` backend is now the default for manifest-driven brightfield runs.
- Another agent is actively cleaning and normalizing the completed segmentation outputs for better traceability.
- Because that cleanup is in progress, this README treats the code and source data as canonical and treats the current `output/segmentation/` tree as an in-progress frozen results area.

Current output holding area:

- `/home/jake/Projects/cellprofiler_test/output/segmentation`

Currently visible top-level files in that holding area:

- `/home/jake/Projects/cellprofiler_test/output/segmentation/batch_run.log`
- `/home/jake/Projects/cellprofiler_test/output/segmentation/brightfield_batch.log`
- `/home/jake/Projects/cellprofiler_test/output/segmentation/cell_contact_sheet.png`
- `/home/jake/Projects/cellprofiler_test/output/segmentation/measurement_histograms.png`
- `/home/jake/Projects/cellprofiler_test/output/segmentation/pmount_blur_comparison.png`
- `/home/jake/Projects/cellprofiler_test/output/segmentation/pmount_mask_overlay.png`
- `/home/jake/Projects/cellprofiler_test/output/segmentation/visual_summary.png`

## Canonical Inputs

Metadata and image sources:

- `/home/jake/Projects/cellprofiler_test/data/metadata/master_image_metadata.csv`
- `/home/jake/Projects/cellprofiler_test/data/brightfield/`
- `/home/jake/Projects/cellprofiler_test/data/pmount/`

Related metadata tables that may be needed for downstream joins:

- `/home/jake/Projects/cellprofiler_test/slide_species_mapping.csv`
- `/home/jake/Projects/cellprofiler_test/slide_glass_mapping.csv`
- `/home/jake/Projects/cellprofiler_test/microscope_to_imageid_desmognathus.csv`

## Canonical Scripts

Canonical wrapper for manifest-driven runs:

- `/home/jake/Projects/cellprofiler_test/cell_size_segmentation_pipeline/run_from_manifest.py`

Rule-based brightfield backend:

- `/home/jake/Projects/cellprofiler_test/scripts/rule_based_cell_size.py`

Existing Cellpose core implementation:

- `/home/jake/Projects/cellprofiler_test/Cellsize_segmentation_cellpose_pipeline/scripts/segment_cells.py`
- `/home/jake/Projects/cellprofiler_test/Cellsize_segmentation_cellpose_pipeline/scripts/run_batch.py`
- `/home/jake/Projects/cellprofiler_test/Cellsize_segmentation_cellpose_pipeline/scripts/tile_selector.py`

Related supporting scripts that exist in the repo but should not be treated as a second canonical production path unless explicitly promoted:

- `/home/jake/Projects/cellprofiler_test/scripts/cellpose_cell_size.py`
- `/home/jake/Projects/cellprofiler_test/scripts/correlate_cells_nuclei.py`
- `/home/jake/Projects/cellprofiler_test/scripts/analyze_results.py`
- `/home/jake/Projects/cellprofiler_test/scripts/visualize_results.py`

## Current Measurement Logic

Default brightfield path (`--backend rules`):

1. Reads a brightfield image tile-by-tile.
2. Reuses the tuned nucleus thresholding logic from the Python IOD pipeline.
3. Grows cell masks from nucleus seeds by watershed inside a thresholded cell foreground mask.
4. Filters by realistic cell geometry and nucleus-within-cell constraints.
5. Measures cell area and cell-level IOD-style values on the original green image.
6. Writes per-image measurement CSVs with tile and mask provenance fields.
7. Writes per-image or per-tile label masks.
8. Writes tile manifests for tiled images.
9. Optionally writes crops and debug overlays.

Legacy brightfield/pmount path (`--backend cellpose`):

1. Reads a brightfield or pmount image.
2. Runs Cellpose instance segmentation.
3. Filters objects by shape.
4. Filters objects by nearest-neighbor isolation.
5. Measures area and IOD-style per-object values.
6. Writes the same measurement/mask/tile artifacts.

## Current Output Contract

Today, the segmentation pipeline is expected to produce some or all of:

- per-image measurement CSVs
- per-image or per-tile label masks
- per-image tile manifest CSVs
- optional crop directories
- optional debug overlays
- batch logs
- summary figures

The intended canonical future location for normalized segmentation outputs is:

- `/home/jake/Projects/cellprofiler_test/output/runs/<run_tag>/cell_size_segmentation/`

Recommended substructure for the cleaned canonical segmentation run:

- `measurements/`
- `masks/`
- `overlays/`
- `crops/`
- `logs/`
- `manifests/`

## Traceability Plan

End-state traceability for this pipeline should be:

1. source image
2. segmentation run manifest row
3. per-image measurement CSV
4. saved label mask or object map
5. optional overlay / crop for visual QC
6. object row used downstream by the IOD pipeline or by size summaries

The minimum reproducible segmentation artifact set should be:

- source filename
- image type
- slide ID
- specimen ID
- species
- segmentation settings used
- per-object measurements CSV
- label mask for the image or tile
- tile manifest with tile coordinates
- overlay image for visual checking

## Relationship To The IOD Pipeline

This pipeline provides the segmentation / object-definition side of the project.

The nucleus IOD estimate pipeline is documented separately in:

- `/home/jake/Projects/cellprofiler_test/nucleus_iod_estimate_pipeline/README.md`

The two pipelines should converge on a shared traceability contract:

- source image path
- image metadata
- object or tile provenance
- mask path
- final analysis row path

## Reproducibility Rules

Use these rules going forward:

- do not create copied subset image folders under `data/`
- do not create new ad hoc output roots for each experiment
- normalize completed runs under `output/runs/<run_tag>/cell_size_segmentation/`
- preserve one manifest per run
- preserve one settings snapshot per run
- invoke the wrapper via `uv run python cell_size_segmentation_pipeline/run_from_manifest.py ...`
- for brightfield rule tuning, preserve any per-image nucleus threshold CSV passed via `--threshold-csv`

Shared data/output docs:

- `/home/jake/Projects/cellprofiler_test/data/README.md`
- `/home/jake/Projects/cellprofiler_test/output/README.md`

## Known Gaps

- the currently cleaned output tree is still in transition
- the completed segmentation run is not yet documented as a single canonical run bundle
- Cellpose mask persistence is now implemented in the canonical wrapper/core, but the historical completed run in `output/segmentation/` has not yet been re-normalized into the shared `output/runs/<run_tag>/` layout

## Immediate Next Step

After the cleanup agent finishes, update this pipeline so it points to:

- the finalized canonical segmentation run root
- the real per-image measurement directory
- the real saved mask directory
- the real tile manifest directory
- the real run manifest and settings snapshot
