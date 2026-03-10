# Nucleus IOD Estimate Pipeline

## Scope

This folder is the canonical documentation anchor for the nucleus-area and nucleus-IOD analysis pipeline.

Repo-level shared architecture:

- `/home/jake/Projects/cellprofiler_test/three_pipelines/README.md`

It documents:

- the current working analysis baseline
- the data sources used to build it
- the scripts that own the QC and estimate logic
- the current traceability coverage
- the remaining gap to full mask-level reproducibility

## Current Status

The current best automated baseline is:

- `/home/jake/Projects/cellprofiler_test/output/qc_report_blockbalanced`

The compact working estimates bundle is:

- `/home/jake/Projects/cellprofiler_test/output/qc_report_blockbalanced/estimates`

The current primary analysis tier is:

- `strict_core`

The current best interpretation is:

- block-balanced dense-region exclusion is active
- dense / clumped regions are automatically downweighted or excluded at the spatial-block level
- species-level outputs are already traceable to image, tile, and nucleus-row level
- mask-level traceability is not yet complete for the legacy run

## Canonical Inputs

Primary source data:

- `/home/jake/Projects/cellprofiler_test/data/metadata/master_image_metadata.csv`
- `/home/jake/Projects/cellprofiler_test/data/brightfield/`
- `/home/jake/Projects/cellprofiler_test/data/pmount/`

Current measurement table used by the legacy QC baseline:

- `/home/jake/Projects/cellprofiler_test/output/nucleus_iod/nucleus_iod_measurements.csv`

Fail-fast validation manifest for the next canonical artifact-producing run:

- `/home/jake/Projects/cellprofiler_test/data/manifests/fail_fast_panel_v1.csv`

## Canonical Scripts

Current analysis / reporting scripts:

- `/home/jake/Projects/cellprofiler_test/nucleus_iod_estimate_pipeline/run_from_manifest.py`
- `/home/jake/Projects/cellprofiler_test/scripts/qc_analysis.py`
- `/home/jake/Projects/cellprofiler_test/scripts/build_final_species_results.py`
- `/home/jake/Projects/cellprofiler_test/scripts/build_traceability_bundle.py`
- `/home/jake/Projects/cellprofiler_test/scripts/build_traceability_viewer.py`
- `/home/jake/Projects/cellprofiler_test/scripts/build_interim_estimate_bundle.py`

Prototype or validation-side measurement scripts that exist in the repo but are not yet the final canonical production engine:

- `/home/jake/Projects/cellprofiler_test/scripts/imagej_nucleus_iod.py`
- `/home/jake/Projects/cellprofiler_test/scripts/nucleus_iod.ijm`
- `/home/jake/Projects/cellprofiler_test/scripts/measure_pmount_iod.py`
- `/home/jake/Projects/cellprofiler_test/scripts/run_reproducible_imagej_bundle.sh`

## Current Working Outputs

Primary analysis outputs:

- `/home/jake/Projects/cellprofiler_test/output/qc_report_blockbalanced/qc_report.html`
- `/home/jake/Projects/cellprofiler_test/output/qc_report_blockbalanced/final_species_results.csv`
- `/home/jake/Projects/cellprofiler_test/output/qc_report_blockbalanced/species_flags.csv`
- `/home/jake/Projects/cellprofiler_test/output/qc_report_blockbalanced/preparation_sensitivity.csv`
- `/home/jake/Projects/cellprofiler_test/output/qc_report_blockbalanced/species_tier_stability.csv`

Working estimates bundle:

- `/home/jake/Projects/cellprofiler_test/output/qc_report_blockbalanced/estimates/species_primary_estimates.csv`
- `/home/jake/Projects/cellprofiler_test/output/qc_report_blockbalanced/estimates/species_state_iod_summary.csv`
- `/home/jake/Projects/cellprofiler_test/output/qc_report_blockbalanced/estimates/image_qc_review.csv`
- `/home/jake/Projects/cellprofiler_test/output/qc_report_blockbalanced/estimates/same_slide_iod_validity.csv`
- `/home/jake/Projects/cellprofiler_test/output/qc_report_blockbalanced/estimates/estimate_outlier_summary.csv`

Visual verification figures:

- `/home/jake/Projects/cellprofiler_test/output/qc_report_blockbalanced/genome_sizes.png`
- `/home/jake/Projects/cellprofiler_test/output/qc_report_blockbalanced/area_boxplots.png`
- `/home/jake/Projects/cellprofiler_test/output/qc_report_blockbalanced/iod_boxplots_by_type.png`
- `/home/jake/Projects/cellprofiler_test/output/qc_report_blockbalanced/per_image_qc.png`
- `/home/jake/Projects/cellprofiler_test/output/qc_report_blockbalanced/area_vs_genome.png`
- `/home/jake/Projects/cellprofiler_test/output/qc_report_blockbalanced/species_area_histograms.png`

## Current Traceability Coverage

Current traceability bundle:

- `/home/jake/Projects/cellprofiler_test/output/qc_report_blockbalanced/traceability`

Key files:

- `/home/jake/Projects/cellprofiler_test/output/qc_report_blockbalanced/traceability/index.html`
- `/home/jake/Projects/cellprofiler_test/output/qc_report_blockbalanced/traceability/species_index.csv`
- `/home/jake/Projects/cellprofiler_test/output/qc_report_blockbalanced/traceability/image_traceability.csv`
- `/home/jake/Projects/cellprofiler_test/output/qc_report_blockbalanced/traceability/tile_traceability.csv`
- `/home/jake/Projects/cellprofiler_test/output/qc_report_blockbalanced/traceability/nucleus_trace_index.csv.gz`
- `/home/jake/Projects/cellprofiler_test/output/qc_report_blockbalanced/traceability/finding_traceability.csv`

Current legacy coverage already supports:

1. species result -> contributing images
2. image -> contributing tile or spatial block
3. tile -> nucleus rows
4. nucleus row -> source image path

Current legacy gap:

- `mask_path = 0`
- `roi_zip_path = 0`
- `tile_manifest_path = 0`

That gap is documented in:

- `/home/jake/Projects/cellprofiler_test/output/qc_report_blockbalanced/traceability/traceability_coverage.json`

## Scientific Guardrails

Important interpretation rules already established in the current analysis:

- `D. fuscus = 16.36 pg` is the reference genome size
- reference anchoring is used within preservation groups
- `brightfield` vs `pmount` is an experimental-state difference, not just a viewer setting difference
- mixed-state same-slide IOD comparisons are not valid direct gold-standard calibration
- dense spatial blocks are excluded automatically before image-level aggregation

## End-State Traceability Plan

The final canonical IOD run must support:

1. species result row
2. image row contributing to that result
3. tile or block contributing to that image summary
4. tile input image path
5. saved mask TIFF path
6. saved ROI ZIP path
7. nucleus rows used or excluded by each QC tier

In short:

- result -> image -> tile -> mask -> nucleus row -> source image

## Canonical Future Run Layout

Future canonical runs should live under:

- `/home/jake/Projects/cellprofiler_test/output/runs/<run_tag>/`

Expected subdirectories:

- `iod/`
- `qc/`
- `traceability/`
- `logs/`
- `manifests/`

The next fail-fast canonical run should be:

- `output/runs/fail_fast_panel_v1/`

and should use:

- `/home/jake/Projects/cellprofiler_test/data/manifests/fail_fast_panel_v1.csv`

## Known Gaps

- the current best baseline is still a legacy run, not the final artifact-complete canonical run
- same-slide IOD rows in the current dataset are mixed image-type comparisons and should not be treated as direct calibration
- mask-level traceability still requires the artifact-producing rerun

## Immediate Next Step

Use the completed segmentation baseline plus the fail-fast panel to produce the next canonical IOD run that writes:

- tile manifests
- tile inputs
- mask TIFFs
- ROI ZIPs
- per-tile and per-image measurement tables
- QC and traceability outputs from the same run root

Shared data/output docs:

- `/home/jake/Projects/cellprofiler_test/data/README.md`
- `/home/jake/Projects/cellprofiler_test/output/README.md`

Recommended invocation:

```bash
uv run python nucleus_iod_estimate_pipeline/run_from_manifest.py \
  --manifest data/manifests/fail_fast_panel_v1.csv \
  --output-dir output/runs/fail_fast_panel_v1_canonical/nucleus_iod \
  --qc-dir output/runs/fail_fast_panel_v1_canonical/qc \
  --trace-dir output/runs/fail_fast_panel_v1_canonical/traceability \
  --run-qc
```
