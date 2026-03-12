#!/usr/bin/env python3
"""Build a species-level HTML report from linked cell-nucleus pairs."""

from __future__ import annotations

import argparse
import html
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats as sp_stats

PROJECT = Path(__file__).resolve().parent.parent
DEFAULT_LINKED_PAIRS_CSV = (
    PROJECT
    / "output"
    / "runs"
    / "mixed_cellpose_yolo_full_dataset_v1"
    / "linkage"
    / "linked_nucleus_pairs.csv.gz"
)
DEFAULT_OUTPUT_DIR = (
    PROJECT
    / "output"
    / "runs"
    / "mixed_cellpose_yolo_full_dataset_v1"
    / "linked_species_stats"
)
DEFAULT_BACKGROUND_CACHE_CSV = (
    PROJECT
    / "output"
    / "runs"
    / "full_dataset_v1"
    / "nucleus_iod"
    / "brightfield"
    / "measurements"
    / "nucleus_iod_measurements.csv"
)
DEFAULT_REFERENCE_SPECIES = "D. fuscus"
DEFAULT_REFERENCE_GENOME_PG = 16.36

METRICS = [
    ("cell_area_um2", "Cell area (um^2)"),
    ("nuc_area_um2", "Nucleus area (um^2)"),
    ("nuc_iod", "Nucleus IOD (image-background normalized)"),
    ("estimated_genome_pg", "Estimated genome size (pg)"),
]
PRIMARY_METRICS = [
    ("cell_area_um2", "Cell area (um^2)"),
    ("nuc_area_um2", "Nucleus area (um^2)"),
    ("estimated_genome_pg", "Estimated genome size (pg)"),
]
NUMERIC_COLUMNS = [
    "cell_area_um2",
    "nuc_area_um2",
    "nuc_area_px",
    "nuc_iod",
    "nuc_i_bg",
    "nc_area_ratio",
    "cytoplasm_area_um2",
    "specimen_id",
]
BOOL_COLUMNS = [
    "has_cell_match",
    "physical_pair_ok",
    "one_to_one_cell",
    "keep_strict_core",
]
REQUIRED_COLUMNS = {
    "species",
    "filename",
    "cell_area_um2",
    "nuc_area_um2",
    "nuc_area_px",
    "nuc_iod",
    "nuc_i_bg",
    "has_cell_match",
    "physical_pair_ok",
    "one_to_one_cell",
    "keep_strict_core",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--linked-pairs-csv", type=Path, default=DEFAULT_LINKED_PAIRS_CSV)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--background-cache-csv", type=Path, default=DEFAULT_BACKGROUND_CACHE_CSV)
    parser.add_argument("--reference-species", default=DEFAULT_REFERENCE_SPECIES)
    parser.add_argument("--reference-genome-pg", type=float, default=DEFAULT_REFERENCE_GENOME_PG)
    return parser.parse_args()


def require_exists(path: Path) -> None:
    if not path.exists():
        raise FileNotFoundError(path)


def parse_bool(series: pd.Series) -> pd.Series:
    if pd.api.types.is_bool_dtype(series):
        return series.fillna(False)
    lowered = series.astype(str).str.strip().str.lower()
    return lowered.isin({"1", "true", "t", "yes", "y"})


def load_pairs(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, keep_default_na=False)
    missing = sorted(REQUIRED_COLUMNS.difference(df.columns))
    if missing:
        raise ValueError(f"Missing required columns in {path}: {missing}")
    for col in NUMERIC_COLUMNS:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")
    for col in BOOL_COLUMNS:
        if col in df.columns:
            df[col] = parse_bool(df[col])
    df["species"] = df["species"].astype(str)
    df["filename"] = df["filename"].astype(str)
    if "specimen_id" not in df.columns:
        df["specimen_id"] = np.nan
    return df


def load_background_reference(path: Path | None) -> dict[str, float]:
    if path is None or not path.exists():
        return {}
    df = pd.read_csv(path, usecols=["filename", "i_bg"], keep_default_na=False)
    df["filename"] = df["filename"].astype(str)
    df["i_bg"] = pd.to_numeric(df["i_bg"], errors="coerce")
    df = df.dropna(subset=["filename", "i_bg"])
    grouped = df.groupby("filename", sort=False)["i_bg"].median()
    return {str(name): float(val) for name, val in grouped.items() if np.isfinite(val) and val > 0}


def normalize_nucleus_iod(
    df: pd.DataFrame,
    background_reference: dict[str, float],
) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, float | int | str]]:
    out = df.copy()
    out["nuc_iod_raw"] = out["nuc_iod"]
    out["nuc_i_bg_raw"] = out["nuc_i_bg"]
    image_max_bg = out.groupby("filename", sort=False)["nuc_i_bg"].max()
    out["nuc_i_bg_reference"] = out["filename"].map(background_reference)
    out["nuc_i_bg_reference_source"] = np.where(
        out["filename"].isin(background_reference),
        "background_cache",
        "linked_pairs_image_max",
    )
    missing_ref = out["nuc_i_bg_reference"].isna()
    out.loc[missing_ref, "nuc_i_bg_reference"] = out.loc[missing_ref, "filename"].map(image_max_bg)
    out["nuc_i_bg_reference"] = pd.to_numeric(out["nuc_i_bg_reference"], errors="coerce")
    valid_ref = out["nuc_i_bg_reference"].gt(0) & out["nuc_i_bg_raw"].gt(0) & out["nuc_area_px"].gt(0)
    out["nuc_iod_background_shift"] = 0.0
    out.loc[valid_ref, "nuc_iod_background_shift"] = (
        out.loc[valid_ref, "nuc_area_px"]
        * np.log10(out.loc[valid_ref, "nuc_i_bg_reference"] / out.loc[valid_ref, "nuc_i_bg_raw"])
    )
    out["nuc_iod"] = out["nuc_iod_raw"] + out["nuc_iod_background_shift"]
    out["nuc_mean_od_raw"] = out.get("nuc_mean_od", np.nan)
    out["nuc_mean_od"] = np.where(
        out["nuc_area_px"].gt(0),
        out["nuc_iod"] / out["nuc_area_px"],
        np.nan,
    )
    out["nuc_i_bg_was_corrected"] = np.abs(out["nuc_iod_background_shift"]) > 1e-9

    audit = (
        out.groupby("filename", sort=True)
        .agg(
            species=("species", "first"),
            n_pairs=("filename", "size"),
            n_unique_raw_i_bg=("nuc_i_bg_raw", "nunique"),
            raw_i_bg_min=("nuc_i_bg_raw", "min"),
            raw_i_bg_median=("nuc_i_bg_raw", "median"),
            raw_i_bg_max=("nuc_i_bg_raw", "max"),
            reference_i_bg=("nuc_i_bg_reference", "first"),
            reference_source=("nuc_i_bg_reference_source", "first"),
            n_corrected_pairs=("nuc_i_bg_was_corrected", "sum"),
            median_iod_shift=("nuc_iod_background_shift", "median"),
            max_abs_iod_shift=("nuc_iod_background_shift", lambda s: float(np.max(np.abs(s))) if len(s) else np.nan),
        )
        .reset_index()
        .sort_values(["n_unique_raw_i_bg", "n_corrected_pairs", "filename"], ascending=[False, False, True])
        .reset_index(drop=True)
    )

    summary = {
        "n_images_in_report": int(out["filename"].nunique()),
        "n_reference_background_images": int(len(background_reference)),
        "n_images_with_mixed_raw_i_bg": int(audit["n_unique_raw_i_bg"].gt(1).sum()),
        "n_pairs_with_background_shift": int(out["nuc_i_bg_was_corrected"].sum()),
        "n_images_without_cache_reference": int(audit["reference_source"].eq("linked_pairs_image_max").sum()),
        "max_abs_iod_shift": float(audit["max_abs_iod_shift"].max()) if len(audit) else np.nan,
    }
    return out, audit, summary


def mode_bin_center(values: np.ndarray) -> float:
    vals = np.asarray(values, dtype=float)
    vals = vals[np.isfinite(vals)]
    if vals.size == 0:
        return np.nan
    if vals.size == 1:
        return float(vals[0])
    q1, q3 = np.percentile(vals, [25, 75])
    iqr = q3 - q1
    data_range = float(vals.max() - vals.min())
    if not np.isfinite(data_range) or data_range <= 0:
        return float(np.median(vals))
    if not np.isfinite(iqr) or iqr <= 0:
        bins = min(24, max(6, int(np.sqrt(vals.size))))
    else:
        bin_width = 2 * iqr * np.power(vals.size, -1 / 3)
        if not np.isfinite(bin_width) or bin_width <= 0:
            bins = min(24, max(6, int(np.sqrt(vals.size))))
        else:
            bins = int(np.ceil(data_range / bin_width))
            bins = max(6, min(80, bins))
    hist, edges = np.histogram(vals, bins=bins)
    idx = int(hist.argmax())
    return float((edges[idx] + edges[idx + 1]) / 2.0)


def summarize_series(series: pd.Series) -> dict[str, float]:
    vals = pd.to_numeric(series, errors="coerce").dropna().to_numpy(dtype=float)
    if vals.size == 0:
        return {
            "n": 0,
            "mean": np.nan,
            "median": np.nan,
            "mode_fd_bin_center": np.nan,
            "std": np.nan,
            "sem": np.nan,
            "mad": np.nan,
            "min": np.nan,
            "q1": np.nan,
            "q3": np.nan,
            "iqr": np.nan,
            "max": np.nan,
            "cv": np.nan,
            "skewness": np.nan,
        }
    q1, q3 = np.percentile(vals, [25, 75])
    mean = float(np.mean(vals))
    std = float(np.std(vals, ddof=1)) if vals.size > 1 else 0.0
    skewness = float(sp_stats.skew(vals, bias=False)) if vals.size >= 3 and std > 0 else np.nan
    return {
        "n": int(vals.size),
        "mean": mean,
        "median": float(np.median(vals)),
        "mode_fd_bin_center": mode_bin_center(vals),
        "std": std,
        "sem": float(std / np.sqrt(vals.size)) if vals.size > 0 else np.nan,
        "mad": float(np.median(np.abs(vals - np.median(vals)))),
        "min": float(np.min(vals)),
        "q1": float(q1),
        "q3": float(q3),
        "iqr": float(q3 - q1),
        "max": float(np.max(vals)),
        "cv": float(std / mean) if mean else np.nan,
        "skewness": skewness,
    }


def corr_pair(df: pd.DataFrame, left: str, right: str, method: str) -> float:
    sub = df[[left, right]].replace([np.inf, -np.inf], np.nan).dropna()
    if len(sub) < 3:
        return np.nan
    if method == "pearson":
        return float(sub[left].corr(sub[right], method="pearson"))
    if method == "spearman":
        return float(sub[left].corr(sub[right], method="spearman"))
    raise ValueError(f"Unsupported method: {method}")


def prepare_analysis_sets(
    df: pd.DataFrame,
    reference_species: str,
    reference_genome_pg: float,
) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, float | str]]:
    positive = (
        df["has_cell_match"]
        & df["physical_pair_ok"]
        & df["one_to_one_cell"]
        & df["cell_area_um2"].gt(0)
        & df["nuc_area_um2"].gt(0)
        & df["nuc_iod"].gt(0)
        & df["species"].ne("")
    )
    one_to_one = df.loc[positive].copy()
    strict = one_to_one.loc[one_to_one["keep_strict_core"]].copy()

    reference_rows = strict.loc[strict["species"] == reference_species].copy()
    reference_set = "strict_core"
    if reference_rows.empty:
        reference_rows = one_to_one.loc[one_to_one["species"] == reference_species].copy()
        reference_set = "one_to_one_physical"
    if reference_rows.empty:
        raise ValueError(
            f"No usable linked rows found for reference species {reference_species!r}."
        )

    reference_median_iod = float(reference_rows["nuc_iod"].median())
    if not np.isfinite(reference_median_iod) or reference_median_iod <= 0:
        raise ValueError(f"Invalid reference median IOD for species {reference_species!r}.")
    scale_pg_per_iod = float(reference_genome_pg / reference_median_iod)

    for frame in (one_to_one, strict):
        frame["estimated_genome_pg"] = frame["nuc_iod"] * scale_pg_per_iod

    calibration = {
        "reference_species": reference_species,
        "reference_genome_pg": float(reference_genome_pg),
        "reference_set": reference_set,
        "reference_rows": int(len(reference_rows)),
        "reference_median_iod": reference_median_iod,
        "scale_pg_per_iod": scale_pg_per_iod,
    }
    return one_to_one, strict, calibration


def summarize_metric_table(df: pd.DataFrame, dataset_label: str) -> pd.DataFrame:
    rows: list[dict[str, float | int | str]] = []
    for species, grp in df.groupby("species", sort=True):
        for metric_key, metric_label in METRICS:
            stats = summarize_series(grp[metric_key])
            rows.append(
                {
                    "dataset": dataset_label,
                    "species": species,
                    "metric_key": metric_key,
                    "metric_label": metric_label,
                    **stats,
                }
            )
    return pd.DataFrame(rows)


def summarize_overview(one_to_one: pd.DataFrame, strict: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, float | int | str]] = []
    species_list = sorted(set(one_to_one["species"]).union(strict["species"]))
    for species in species_list:
        all_grp = one_to_one.loc[one_to_one["species"] == species].copy()
        strict_grp = strict.loc[strict["species"] == species].copy()
        n_all = int(len(all_grp))
        n_strict = int(len(strict_grp))
        rows.append(
            {
                "species": species,
                "n_pairs_one_to_one": n_all,
                "n_pairs_strict_core": n_strict,
                "strict_retention_fraction": float(n_strict / n_all) if n_all else np.nan,
                "n_images": int(all_grp["filename"].nunique()) if n_all else 0,
                "n_specimens": int(all_grp["specimen_id"].dropna().nunique()) if n_all else 0,
                "median_cell_area_um2_strict": float(strict_grp["cell_area_um2"].median()) if n_strict else np.nan,
                "median_nuc_area_um2_strict": float(strict_grp["nuc_area_um2"].median()) if n_strict else np.nan,
                "median_estimated_genome_pg_strict": float(strict_grp["estimated_genome_pg"].median()) if n_strict else np.nan,
                "median_nc_ratio_strict": float(strict_grp["nc_area_ratio"].median()) if n_strict else np.nan,
                "pearson_cell_vs_nucleus_area_strict": corr_pair(
                    strict_grp, "cell_area_um2", "nuc_area_um2", "pearson"
                ),
                "spearman_cell_vs_nucleus_area_strict": corr_pair(
                    strict_grp, "cell_area_um2", "nuc_area_um2", "spearman"
                ),
                "pearson_cell_vs_genome_pg_strict": corr_pair(
                    strict_grp, "cell_area_um2", "estimated_genome_pg", "pearson"
                ),
                "spearman_cell_vs_genome_pg_strict": corr_pair(
                    strict_grp, "cell_area_um2", "estimated_genome_pg", "spearman"
                ),
                "pearson_nucleus_area_vs_genome_pg_strict": corr_pair(
                    strict_grp, "nuc_area_um2", "estimated_genome_pg", "pearson"
                ),
                "spearman_nucleus_area_vs_genome_pg_strict": corr_pair(
                    strict_grp, "nuc_area_um2", "estimated_genome_pg", "spearman"
                ),
            }
        )
    return pd.DataFrame(rows).sort_values("n_pairs_strict_core", ascending=False).reset_index(drop=True)


def summarize_global(
    one_to_one: pd.DataFrame,
    strict: pd.DataFrame,
    calibration: dict[str, float | str],
    background_audit_summary: dict[str, float | int | str],
) -> dict[str, float | int | str]:
    summary: dict[str, float | int | str] = {
        "n_species": int(one_to_one["species"].nunique()),
        "n_pairs_one_to_one": int(len(one_to_one)),
        "n_pairs_strict_core": int(len(strict)),
        "strict_retention_fraction": float(len(strict) / len(one_to_one)) if len(one_to_one) else np.nan,
        "median_cell_area_um2_strict": float(strict["cell_area_um2"].median()) if len(strict) else np.nan,
        "median_nuc_area_um2_strict": float(strict["nuc_area_um2"].median()) if len(strict) else np.nan,
        "median_estimated_genome_pg_strict": float(strict["estimated_genome_pg"].median()) if len(strict) else np.nan,
        "median_nc_ratio_strict": float(strict["nc_area_ratio"].median()) if len(strict) else np.nan,
        "pearson_cell_vs_nucleus_area_strict": corr_pair(strict, "cell_area_um2", "nuc_area_um2", "pearson"),
        "spearman_cell_vs_nucleus_area_strict": corr_pair(strict, "cell_area_um2", "nuc_area_um2", "spearman"),
        "pearson_cell_vs_genome_pg_strict": corr_pair(strict, "cell_area_um2", "estimated_genome_pg", "pearson"),
        "spearman_cell_vs_genome_pg_strict": corr_pair(strict, "cell_area_um2", "estimated_genome_pg", "spearman"),
        "pearson_nucleus_area_vs_genome_pg_strict": corr_pair(strict, "nuc_area_um2", "estimated_genome_pg", "pearson"),
        "spearman_nucleus_area_vs_genome_pg_strict": corr_pair(strict, "nuc_area_um2", "estimated_genome_pg", "spearman"),
    }
    summary.update(calibration)
    summary.update(background_audit_summary)
    return summary


def fmt_num(value: float | int | str, digits: int = 2) -> str:
    if isinstance(value, str):
        return value
    if value is None or (isinstance(value, float) and not np.isfinite(value)):
        return "NA"
    if isinstance(value, (int, np.integer)):
        return f"{int(value):,}"
    return f"{float(value):,.{digits}f}"


def metric_table_for_html(metric_stats: pd.DataFrame, metric_key: str) -> str:
    sub = (
        metric_stats.loc[metric_stats["metric_key"] == metric_key]
        .sort_values("median", ascending=False)
        .reset_index(drop=True)
    )
    rows = []
    for _, row in sub.iterrows():
        rows.append(
            "<tr>"
            f"<td>{html.escape(str(row['species']))}</td>"
            f"<td>{fmt_num(row['n'], 0)}</td>"
            f"<td>{fmt_num(row['mean'])}</td>"
            f"<td>{fmt_num(row['median'])}</td>"
            f"<td>{fmt_num(row['mode_fd_bin_center'])}</td>"
            f"<td>{fmt_num(row['q1'])}</td>"
            f"<td>{fmt_num(row['q3'])}</td>"
            f"<td>{fmt_num(row['iqr'])}</td>"
            f"<td>{fmt_num(row['std'])}</td>"
            "</tr>"
        )
    return "".join(rows)


def overview_rows_for_html(overview: pd.DataFrame) -> str:
    rows = []
    for _, row in overview.iterrows():
        rows.append(
            "<tr>"
            f"<td>{html.escape(str(row['species']))}</td>"
            f"<td>{fmt_num(row['n_pairs_one_to_one'], 0)}</td>"
            f"<td>{fmt_num(row['n_pairs_strict_core'], 0)}</td>"
            f"<td>{fmt_num(row['strict_retention_fraction'], 3)}</td>"
            f"<td>{fmt_num(row['median_cell_area_um2_strict'])}</td>"
            f"<td>{fmt_num(row['median_nuc_area_um2_strict'])}</td>"
            f"<td>{fmt_num(row['median_estimated_genome_pg_strict'], 3)}</td>"
            f"<td>{fmt_num(row['spearman_cell_vs_nucleus_area_strict'], 3)}</td>"
            f"<td>{fmt_num(row['spearman_cell_vs_genome_pg_strict'], 3)}</td>"
            "</tr>"
        )
    return "".join(rows)


def audit_rows_for_html(audit: pd.DataFrame) -> str:
    rows = []
    for _, row in audit.iterrows():
        rows.append(
            "<tr>"
            f"<td>{html.escape(str(row['filename']))}</td>"
            f"<td>{html.escape(str(row['species']))}</td>"
            f"<td>{fmt_num(row['n_pairs'], 0)}</td>"
            f"<td>{fmt_num(row['n_unique_raw_i_bg'], 0)}</td>"
            f"<td>{fmt_num(row['raw_i_bg_min'])}</td>"
            f"<td>{fmt_num(row['raw_i_bg_median'])}</td>"
            f"<td>{fmt_num(row['raw_i_bg_max'])}</td>"
            f"<td>{fmt_num(row['reference_i_bg'])}</td>"
            f"<td>{html.escape(str(row['reference_source']))}</td>"
            f"<td>{fmt_num(row['n_corrected_pairs'], 0)}</td>"
            f"<td>{fmt_num(row['median_iod_shift'], 3)}</td>"
            "</tr>"
        )
    return "".join(rows)


def plot_species_counts(overview: pd.DataFrame, output_path: Path) -> None:
    sub = overview.sort_values("n_pairs_strict_core", ascending=True)
    fig_h = max(7.0, 0.35 * len(sub) + 1.5)
    fig, ax = plt.subplots(figsize=(11, fig_h))
    y = np.arange(len(sub))
    ax.barh(y, sub["n_pairs_one_to_one"], color="#d8c1aa", label="One-to-one physical")
    ax.barh(y, sub["n_pairs_strict_core"], color="#8c5632", label="Strict core")
    ax.set_yticks(y)
    ax.set_yticklabels(sub["species"])
    ax.set_xlabel("Linked pairs")
    ax.set_title("Per-species linked pair counts")
    ax.grid(axis="x", alpha=0.25)
    ax.legend(frameon=False, loc="lower right")
    fig.tight_layout()
    fig.savefig(output_path, dpi=180, bbox_inches="tight")
    plt.close(fig)


def plot_metric_intervals(
    metric_stats: pd.DataFrame,
    metric_key: str,
    title: str,
    xlabel: str,
    output_path: Path,
) -> None:
    sub = (
        metric_stats.loc[metric_stats["metric_key"] == metric_key]
        .sort_values("median", ascending=True)
        .reset_index(drop=True)
    )
    fig_h = max(7.0, 0.35 * len(sub) + 1.5)
    fig, ax = plt.subplots(figsize=(12, fig_h))
    y = np.arange(len(sub))
    ax.hlines(y, sub["q1"], sub["q3"], color="#c58c63", linewidth=5, alpha=0.95, label="IQR")
    ax.scatter(sub["mean"], y, s=44, marker="D", color="#d8a94f", label="Mean", zorder=3)
    ax.scatter(sub["mode_fd_bin_center"], y, s=52, marker="^", color="#476b6b", label="Mode (FD bin center)", zorder=3)
    ax.scatter(sub["median"], y, s=48, marker="o", color="#8c5632", label="Median", zorder=4)
    ax.set_yticks(y)
    ax.set_yticklabels(sub["species"])
    ax.set_xlabel(xlabel)
    ax.set_title(title)
    ax.grid(axis="x", alpha=0.25)
    ax.legend(frameon=False, loc="lower right")
    fig.tight_layout()
    fig.savefig(output_path, dpi=180, bbox_inches="tight")
    plt.close(fig)


def plot_species_median_scatter(
    overview: pd.DataFrame,
    x_col: str,
    y_col: str,
    title: str,
    xlabel: str,
    ylabel: str,
    output_path: Path,
) -> None:
    sub = overview[[x_col, y_col, "species", "n_pairs_strict_core"]].replace([np.inf, -np.inf], np.nan).dropna()
    fig, ax = plt.subplots(figsize=(10.5, 8.0))
    sizes = 30 + 0.18 * sub["n_pairs_strict_core"].to_numpy(dtype=float)
    ax.scatter(
        sub[x_col],
        sub[y_col],
        s=sizes,
        color="#8c5632",
        alpha=0.78,
        edgecolors="#f3e8dc",
        linewidths=0.8,
    )
    for _, row in sub.iterrows():
        ax.annotate(
            row["species"],
            (row[x_col], row[y_col]),
            xytext=(5, 4),
            textcoords="offset points",
            fontsize=8,
            color="#3e2f22",
        )
    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig(output_path, dpi=180, bbox_inches="tight")
    plt.close(fig)


def build_html(
    summary: dict[str, float | int | str],
    overview: pd.DataFrame,
    strict_metrics: pd.DataFrame,
    background_audit: pd.DataFrame,
) -> str:
    figures = [
        ("Species linked pair counts", "species_pair_counts.png"),
        ("Cell area mean / median / mode / IQR", "species_cell_area_stats.png"),
        ("Nucleus area mean / median / mode / IQR", "species_nucleus_area_stats.png"),
        ("Estimated genome size mean / median / mode / IQR", "species_genome_size_stats.png"),
        ("Species median cell size vs median genome size", "species_cell_vs_genome_scatter.png"),
        ("Species median cell size vs median nucleus size", "species_cell_vs_nucleus_scatter.png"),
    ]

    figure_cards = "".join(
        "<div class='card figure'>"
        f"<a href='{filename}'><img src='{filename}' alt='{html.escape(label)}'></a>"
        f"<p><a href='{filename}'>{html.escape(label)}</a></p>"
        "</div>"
        for label, filename in figures
    )

    metric_sections = []
    for metric_key, metric_label in METRICS:
        metric_sections.append(
            f"""
            <section>
              <h2>{html.escape(metric_label)} by species (strict core)</h2>
              <table>
                <tr>
                  <th>Species</th>
                  <th>n</th>
                  <th>Mean</th>
                  <th>Median</th>
                  <th>Mode</th>
                  <th>Q1</th>
                  <th>Q3</th>
                  <th>IQR</th>
                  <th>Std</th>
                </tr>
                {metric_table_for_html(strict_metrics, metric_key)}
              </table>
            </section>
            """
        )

    return f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>Linked Species Statistics Report</title>
  <style>
    :root {{
      --bg: #f3ecdf;
      --panel: #fffcf8;
      --ink: #2e261d;
      --muted: #70604d;
      --line: #d7c6b1;
      --accent: #8c5632;
      --accent-soft: #c58c63;
    }}
    body {{
      margin: 0;
      color: var(--ink);
      background:
        radial-gradient(circle at top left, rgba(140, 86, 50, 0.12), transparent 28%),
        linear-gradient(180deg, #f8f3eb 0%, var(--bg) 100%);
      font: 15px/1.55 Georgia, serif;
    }}
    .wrap {{
      max-width: 1520px;
      margin: 0 auto;
      padding: 24px;
    }}
    .hero, .card, table {{
      background: var(--panel);
      border: 1px solid var(--line);
      border-radius: 18px;
      box-shadow: 0 12px 28px rgba(43, 31, 19, 0.06);
    }}
    .hero, .card {{
      padding: 18px;
    }}
    .grid {{
      display: grid;
      grid-template-columns: repeat(auto-fit, minmax(220px, 1fr));
      gap: 14px;
      margin: 18px 0 24px;
    }}
    .metric-label {{
      color: var(--muted);
      text-transform: uppercase;
      letter-spacing: 0.08em;
      font-size: 12px;
    }}
    .metric-value {{
      margin-top: 10px;
      color: var(--accent);
      font-size: 28px;
    }}
    .links a {{
      display: inline-block;
      margin: 0 10px 10px 0;
      padding: 8px 12px;
      border-radius: 999px;
      background: #efe1d2;
      border: 1px solid #dfc9b5;
      color: var(--ink);
      text-decoration: none;
    }}
    .fig-grid {{
      display: grid;
      grid-template-columns: repeat(auto-fit, minmax(320px, 1fr));
      gap: 14px;
      margin-bottom: 24px;
    }}
    img {{
      width: 100%;
      height: auto;
      display: block;
      border-radius: 12px;
      border: 1px solid #e4d7ca;
      background: #faf5ef;
    }}
    table {{
      width: 100%;
      border-collapse: collapse;
      overflow: hidden;
      margin-bottom: 24px;
    }}
    th, td {{
      padding: 8px 10px;
      border-bottom: 1px solid var(--line);
      text-align: left;
      font-size: 13px;
    }}
    th {{
      background: #f1e3d5;
      position: sticky;
      top: 0;
    }}
    .note {{
      color: var(--muted);
    }}
  </style>
</head>
<body>
  <div class="wrap">
    <section class="hero">
      <h1>Linked Species Statistics Report</h1>
      <p>Species-level descriptive statistics derived from linked one-to-one cell and nucleus pairs. The primary tables and figures use the <strong>strict-core</strong> subset so cell size, nucleus size, and nucleus-content estimates stay tied to the same QC-filtered pair rows.</p>
      <p class="note">Before any species summary is computed, linked nucleus IOD values are renormalized to a single per-image background reference so the optical-density calculation follows the same incident-light baseline within each source image. Genome-size estimates are then inferred from the normalized nucleus IOD and scaled so that <strong>{html.escape(str(summary['reference_species']))} median {html.escape(str(summary['reference_set']))} normalized IOD maps to {fmt_num(summary['reference_genome_pg'], 2)} pg</strong>. Mode is reported as the peak histogram bin center using Freedman-Diaconis binning because these measurements are continuous.</p>
      <div class="links">
        <a href="linked_pairs_with_estimated_genome.csv.gz">linked pairs with genome estimate</a>
        <a href="iod_background_audit.csv">IOD background audit</a>
        <a href="species_overview.csv">species overview</a>
        <a href="species_metric_stats_strict_core.csv">strict-core metric stats</a>
        <a href="species_metric_stats_one_to_one.csv">one-to-one metric stats</a>
        <a href="overall_metric_stats_strict_core.csv">overall strict-core stats</a>
        <a href="overall_metric_stats_one_to_one.csv">overall one-to-one stats</a>
        <a href="summary.json">summary json</a>
      </div>
    </section>

    <section class="grid">
      <div class="card"><div class="metric-label">Species</div><div class="metric-value">{fmt_num(summary['n_species'], 0)}</div></div>
      <div class="card"><div class="metric-label">One-to-one pairs</div><div class="metric-value">{fmt_num(summary['n_pairs_one_to_one'], 0)}</div></div>
      <div class="card"><div class="metric-label">Strict-core pairs</div><div class="metric-value">{fmt_num(summary['n_pairs_strict_core'], 0)}</div></div>
      <div class="card"><div class="metric-label">Strict retention</div><div class="metric-value">{fmt_num(summary['strict_retention_fraction'], 3)}</div></div>
      <div class="card"><div class="metric-label">Median strict cell area</div><div class="metric-value">{fmt_num(summary['median_cell_area_um2_strict'])}</div></div>
      <div class="card"><div class="metric-label">Median strict nucleus area</div><div class="metric-value">{fmt_num(summary['median_nuc_area_um2_strict'])}</div></div>
      <div class="card"><div class="metric-label">Median strict genome size</div><div class="metric-value">{fmt_num(summary['median_estimated_genome_pg_strict'], 3)}</div></div>
      <div class="card"><div class="metric-label">Spearman cell vs genome</div><div class="metric-value">{fmt_num(summary['spearman_cell_vs_genome_pg_strict'], 3)}</div></div>
      <div class="card"><div class="metric-label">Images with mixed raw I_bg</div><div class="metric-value">{fmt_num(summary['n_images_with_mixed_raw_i_bg'], 0)}</div></div>
      <div class="card"><div class="metric-label">Pairs background-shifted</div><div class="metric-value">{fmt_num(summary['n_pairs_with_background_shift'], 0)}</div></div>
      <div class="card"><div class="metric-label">Max abs IOD shift</div><div class="metric-value">{fmt_num(summary['max_abs_iod_shift'], 3)}</div></div>
    </section>

    <section class="fig-grid">
      {figure_cards}
    </section>

    <section>
      <h2>Species overview</h2>
      <table>
        <tr>
          <th>Species</th>
          <th>One-to-one pairs</th>
          <th>Strict pairs</th>
          <th>Strict retention</th>
          <th>Median cell area</th>
          <th>Median nucleus area</th>
          <th>Median genome size</th>
          <th>Spearman cell vs nucleus</th>
          <th>Spearman cell vs genome</th>
        </tr>
        {overview_rows_for_html(overview)}
      </table>
    </section>

    <section>
      <h2>Background normalization audit</h2>
      <p class="note">Rows below show the source-image background spread observed in the linked pairs file and the per-image reference I_bg used to renormalize nucleus IOD before genome-size estimation. This is the main computational check against the densitometry method.</p>
      <table>
        <tr>
          <th>Filename</th>
          <th>Species</th>
          <th>Pairs</th>
          <th>Unique raw I_bg</th>
          <th>Raw min</th>
          <th>Raw median</th>
          <th>Raw max</th>
          <th>Reference I_bg</th>
          <th>Reference source</th>
          <th>Corrected pairs</th>
          <th>Median IOD shift</th>
        </tr>
        {audit_rows_for_html(background_audit.head(24))}
      </table>
    </section>

    {''.join(metric_sections)}
  </div>
</body>
</html>
"""


def main() -> None:
    args = parse_args()
    require_exists(args.linked_pairs_csv)
    if args.background_cache_csv is not None and args.background_cache_csv.exists():
        require_exists(args.background_cache_csv)
    args.output_dir.mkdir(parents=True, exist_ok=True)

    pairs = load_pairs(args.linked_pairs_csv)
    background_reference = load_background_reference(args.background_cache_csv)
    pairs, background_audit, background_audit_summary = normalize_nucleus_iod(
        pairs,
        background_reference=background_reference,
    )
    one_to_one, strict, calibration = prepare_analysis_sets(
        pairs,
        reference_species=args.reference_species,
        reference_genome_pg=args.reference_genome_pg,
    )

    one_to_one.to_csv(args.output_dir / "linked_pairs_with_estimated_genome.csv.gz", index=False)

    species_metric_stats_one = summarize_metric_table(one_to_one, "one_to_one_physical")
    species_metric_stats_strict = summarize_metric_table(strict, "strict_core")
    overall_metric_stats_one = pd.DataFrame(
        [{"dataset": "one_to_one_physical", "species": "ALL", "metric_key": key, "metric_label": label, **summarize_series(one_to_one[key])} for key, label in METRICS]
    )
    overall_metric_stats_strict = pd.DataFrame(
        [{"dataset": "strict_core", "species": "ALL", "metric_key": key, "metric_label": label, **summarize_series(strict[key])} for key, label in METRICS]
    )
    overview = summarize_overview(one_to_one, strict)
    summary = summarize_global(one_to_one, strict, calibration, background_audit_summary)

    species_metric_stats_one.to_csv(
        args.output_dir / "species_metric_stats_one_to_one.csv",
        index=False,
        float_format="%.6f",
    )
    species_metric_stats_strict.to_csv(
        args.output_dir / "species_metric_stats_strict_core.csv",
        index=False,
        float_format="%.6f",
    )
    overall_metric_stats_one.to_csv(
        args.output_dir / "overall_metric_stats_one_to_one.csv",
        index=False,
        float_format="%.6f",
    )
    overall_metric_stats_strict.to_csv(
        args.output_dir / "overall_metric_stats_strict_core.csv",
        index=False,
        float_format="%.6f",
    )
    background_audit.to_csv(args.output_dir / "iod_background_audit.csv", index=False, float_format="%.6f")
    overview.to_csv(args.output_dir / "species_overview.csv", index=False, float_format="%.6f")
    (args.output_dir / "summary.json").write_text(json.dumps(summary, indent=2))

    plot_species_counts(overview, args.output_dir / "species_pair_counts.png")
    for metric_key, metric_label in PRIMARY_METRICS:
        name = {
            "cell_area_um2": "species_cell_area_stats.png",
            "nuc_area_um2": "species_nucleus_area_stats.png",
            "estimated_genome_pg": "species_genome_size_stats.png",
        }[metric_key]
        plot_metric_intervals(
            species_metric_stats_strict,
            metric_key=metric_key,
            title=f"{metric_label} by species (strict core)",
            xlabel=metric_label,
            output_path=args.output_dir / name,
        )

    plot_species_median_scatter(
        overview,
        x_col="median_cell_area_um2_strict",
        y_col="median_estimated_genome_pg_strict",
        title="Species median cell size vs estimated genome size",
        xlabel="Median strict-core cell area (um^2)",
        ylabel="Median strict-core estimated genome size (pg)",
        output_path=args.output_dir / "species_cell_vs_genome_scatter.png",
    )
    plot_species_median_scatter(
        overview,
        x_col="median_cell_area_um2_strict",
        y_col="median_nuc_area_um2_strict",
        title="Species median cell size vs nucleus size",
        xlabel="Median strict-core cell area (um^2)",
        ylabel="Median strict-core nucleus area (um^2)",
        output_path=args.output_dir / "species_cell_vs_nucleus_scatter.png",
    )

    html_text = build_html(summary, overview, species_metric_stats_strict, background_audit)
    (args.output_dir / "index.html").write_text(html_text)

    print(f"One-to-one physical pairs: {len(one_to_one)}")
    print(f"Strict-core pairs: {len(strict)}")
    print(f"Reference species: {calibration['reference_species']}")
    print(f"Reference median IOD: {calibration['reference_median_iod']:.6f}")
    print(f"Images with mixed raw I_bg: {background_audit_summary['n_images_with_mixed_raw_i_bg']}")
    print(f"Pairs background-shifted: {background_audit_summary['n_pairs_with_background_shift']}")
    print(f"Wrote report: {args.output_dir / 'index.html'}")


if __name__ == "__main__":
    main()
