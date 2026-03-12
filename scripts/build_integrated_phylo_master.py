#!/usr/bin/env python3
"""Build an integrated comparative master table inside this repo.

This merges the current bg-clean linked cell/nucleus/genome summaries with the
older Desmognathus_TE multi-trait master table so downstream phylogenetic
analyses can run locally without mutating the sibling workspace.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OVERVIEW = (
    PROJECT_ROOT
    / "output"
    / "runs"
    / "mixed_cellpose_yolo_full_dataset_v1_bgclean"
    / "linked_species_stats"
    / "species_overview.csv"
)
DEFAULT_LEGACY_MASTER = (
    PROJECT_ROOT.parent
    / "Desmognathus_TE"
    / "path_analysis"
    / "data"
    / "derived"
    / "path_input_master.csv"
)
DEFAULT_OUTPUT_DIR = (
    PROJECT_ROOT
    / "output"
    / "runs"
    / "mixed_cellpose_yolo_full_dataset_v1_bgclean"
    / "integrated_multinode_phylo"
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--overview", type=Path, default=DEFAULT_OVERVIEW)
    parser.add_argument("--legacy-master", type=Path, default=DEFAULT_LEGACY_MASTER)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    return parser.parse_args()


def standardize_species(label: str) -> str:
    text = str(label).strip()
    if text.startswith("D. "):
        text = text[3:]
    if text.startswith("Desmognathus "):
        text = text.split(" ", 1)[1]
    return text.strip().lower()


def require_exists(path: Path) -> None:
    if not path.exists():
        raise FileNotFoundError(path)


def load_overview(path: Path) -> pd.DataFrame:
    overview = pd.read_csv(path)
    overview["species_key"] = overview["species"].map(standardize_species)
    overview = overview.rename(
        columns={
            "median_estimated_genome_pg_strict": "genome_size_pg_current",
            "median_cell_area_um2_strict": "morph_cell_area_um2_current",
            "median_nuc_area_um2_strict": "morph_nucleus_area_um2_current",
            "median_nc_ratio_strict": "morph_nc_ratio_current",
            "n_pairs_strict_core": "cf_n_pairs_strict_core_current",
            "n_pairs_one_to_one": "cf_n_pairs_one_to_one_current",
            "n_images": "cf_n_images_current",
            "n_specimens": "cf_n_specimens_current",
            "strict_retention_fraction": "cf_strict_retention_fraction_current",
        }
    )
    return overview


def load_legacy_master(path: Path) -> pd.DataFrame:
    master = pd.read_csv(path, low_memory=False)
    master["species_key"] = master["species"].map(standardize_species)
    return master


def merge_master(master: pd.DataFrame, overview: pd.DataFrame) -> pd.DataFrame:
    merged = master.merge(
        overview[
            [
                "species_key",
                "genome_size_pg_current",
                "morph_cell_area_um2_current",
                "morph_nucleus_area_um2_current",
                "morph_nc_ratio_current",
                "cf_n_pairs_strict_core_current",
                "cf_n_pairs_one_to_one_current",
                "cf_n_images_current",
                "cf_n_specimens_current",
                "cf_strict_retention_fraction_current",
            ]
        ],
        on="species_key",
        how="left",
    )

    current_updates = {
        "genome_size_pg": "genome_size_pg_current",
        "morph_cell_area_um2": "morph_cell_area_um2_current",
        "morph_nucleus_area_um2": "morph_nucleus_area_um2_current",
        "morph_nc_ratio": "morph_nc_ratio_current",
    }
    for target_col, current_col in current_updates.items():
        legacy_col = f"{target_col}_legacy"
        merged[legacy_col] = merged[target_col]
        merged[target_col] = merged[current_col].combine_first(merged[target_col])

    merged["has_genome"] = merged["genome_size_pg"].notna()
    merged["has_morphology"] = (
        merged["morph_cell_area_um2"].notna()
        & merged["morph_nucleus_area_um2"].notna()
        & merged["morph_nc_ratio"].notna()
    )
    merged["current_cellprofiler_bridge"] = merged["genome_size_pg_current"].notna()
    merged["current_cellprofiler_source"] = (
        "mixed_cellpose_yolo_full_dataset_v1_bgclean"
    )
    merged.loc[~merged["current_cellprofiler_bridge"], "current_cellprofiler_source"] = pd.NA
    return merged


def build_overlap_summary(merged: pd.DataFrame) -> pd.DataFrame:
    panel_defs = {
        "observed_genome_morphology": [
            "genome_size_pg",
            "morph_cell_area_um2",
            "morph_nucleus_area_um2",
        ],
        "observed_te_genome": [
            "genome_size_pg",
            "order_pielou",
            "ltr_line_logratio",
        ],
        "observed_te_genome_ectopic": [
            "genome_size_pg",
            "order_pielou",
            "ltr_line_logratio",
            "ectopic_mean_ratio",
        ],
        "observed_te_genome_ltr_history": [
            "genome_size_pg",
            "order_pielou",
            "ltr_line_logratio",
            "ltr_history_age_central_mya",
        ],
        "observed_te_genome_morphology": [
            "genome_size_pg",
            "order_pielou",
            "ltr_line_logratio",
            "morph_cell_area_um2",
            "morph_nucleus_area_um2",
        ],
        "observed_te_genome_organismal": [
            "genome_size_pg",
            "order_pielou",
            "ltr_line_logratio",
            "body_size_proxy_mm",
            "aquaticity_index",
        ],
    }

    rows: list[dict[str, object]] = []
    for panel, cols in panel_defs.items():
        subset = merged.loc[
            merged["has_tree_tip"].fillna(False) & merged[cols].notna().all(axis=1),
            "species",
        ].sort_values()
        rows.append(
            {
                "panel_name": panel,
                "n_species": int(len(subset)),
                "species_list": ";".join(subset.astype(str)),
            }
        )
    return pd.DataFrame(rows)


def main() -> None:
    args = parse_args()
    require_exists(args.overview)
    require_exists(args.legacy_master)
    args.output_dir.mkdir(parents=True, exist_ok=True)

    overview = load_overview(args.overview)
    master = load_legacy_master(args.legacy_master)
    merged = merge_master(master, overview)
    overlap = build_overlap_summary(merged)

    merged.to_csv(args.output_dir / "integrated_master_observed.csv", index=False)
    overlap.to_csv(args.output_dir / "integrated_overlap_observed.csv", index=False)

    current_bridge = merged.loc[merged["current_cellprofiler_bridge"], [
        "species",
        "genome_size_pg",
        "genome_size_pg_legacy",
        "morph_cell_area_um2",
        "morph_cell_area_um2_legacy",
        "morph_nucleus_area_um2",
        "morph_nucleus_area_um2_legacy",
        "morph_nc_ratio",
        "morph_nc_ratio_legacy",
        "cf_n_pairs_strict_core_current",
        "cf_n_images_current",
    ]].sort_values("species")
    current_bridge.to_csv(args.output_dir / "cellprofiler_current_vs_legacy_bridge.csv", index=False)

    print(f"Wrote master: {args.output_dir / 'integrated_master_observed.csv'}")
    print(f"Wrote overlap summary: {args.output_dir / 'integrated_overlap_observed.csv'}")


if __name__ == "__main__":
    main()
