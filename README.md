# CellProfiler Tools

This repo is now trimmed to the active salamander brightfield workflow:

- Cell segmentation with the custom Cellpose cell model
- Nucleus segmentation with the custom YOLO nucleus model
- Nucleus IOD measurement
- Cell+nucleus linkage with full traceability

The main workflow document is:

- [docs/MIXED_CELLPOSE_YOLO_WORKFLOW.md](/home/jake/Projects/cellprofiler_test/docs/MIXED_CELLPOSE_YOLO_WORKFLOW.md)

The active code and artifact inventory is:

- [docs/ACTIVE_WORKFLOW_INVENTORY.md](/home/jake/Projects/cellprofiler_test/docs/ACTIVE_WORKFLOW_INVENTORY.md)

Key active entrypoints:

- [cell_size_segmentation_pipeline/run_from_manifest.py](/home/jake/Projects/cellprofiler_test/cell_size_segmentation_pipeline/run_from_manifest.py)
- [nucleus_iod_estimate_pipeline/run_from_manifest.py](/home/jake/Projects/cellprofiler_test/nucleus_iod_estimate_pipeline/run_from_manifest.py)
- [scripts/run_yolo_tile_measurements.py](/home/jake/Projects/cellprofiler_test/scripts/run_yolo_tile_measurements.py)
- [scripts/build_cell_nucleus_linkage_report.py](/home/jake/Projects/cellprofiler_test/scripts/build_cell_nucleus_linkage_report.py)
- [scripts/build_species_linked_stats_report.py](/home/jake/Projects/cellprofiler_test/scripts/build_species_linked_stats_report.py)

Current full-dataset mixed outputs:

- [output/runs/full_dataset_v1](/home/jake/Projects/cellprofiler_test/output/runs/full_dataset_v1)
- [output/runs/mixed_cellpose_yolo_full_dataset_v1](/home/jake/Projects/cellprofiler_test/output/runs/mixed_cellpose_yolo_full_dataset_v1)
