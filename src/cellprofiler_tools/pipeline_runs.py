from __future__ import annotations

import csv
import json
import subprocess
import time
from collections import Counter
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

PROJECT = Path(__file__).resolve().parents[2]
DATA_ROOT = PROJECT / "data"
OUTPUT_ROOT = PROJECT / "output"
RUNS_ROOT = OUTPUT_ROOT / "runs"

REQUIRED_MANIFEST_COLUMNS = {
    "filename",
    "image_type",
    "species",
    "slide_id",
    "specimen_id",
}


@dataclass(frozen=True)
class CanonicalRunPaths:
    run_tag: str
    root: Path
    cell_size_segmentation: Path
    nucleus_iod: Path
    qc: Path
    traceability: Path
    logs: Path
    manifests: Path

    def as_dict(self) -> dict[str, str]:
        return {key: str(value) for key, value in asdict(self).items()}


def get_run_paths(run_tag: str) -> CanonicalRunPaths:
    root = RUNS_ROOT / run_tag
    return CanonicalRunPaths(
        run_tag=run_tag,
        root=root,
        cell_size_segmentation=root / "cell_size_segmentation",
        nucleus_iod=root / "nucleus_iod",
        qc=root / "qc",
        traceability=root / "traceability",
        logs=root / "logs",
        manifests=root / "manifests",
    )


def ensure_run_directories(paths: CanonicalRunPaths) -> list[Path]:
    dirs = [
        paths.root,
        paths.cell_size_segmentation,
        paths.nucleus_iod,
        paths.qc,
        paths.traceability,
        paths.logs,
        paths.manifests,
    ]
    for directory in dirs:
        directory.mkdir(parents=True, exist_ok=True)
    return dirs


def load_metadata_lookup() -> dict[tuple[str, str], dict[str, str]]:
    metadata_path = DATA_ROOT / "metadata" / "master_image_metadata.csv"
    with open(metadata_path, newline="") as handle:
        rows = list(csv.DictReader(handle))
    return {
        (row["image_type"], row["filename"]): row
        for row in rows
        if row.get("filename") and row.get("image_type")
    }


def load_manifest_rows(manifest_path: Path) -> list[dict[str, str]]:
    manifest_path = Path(manifest_path)
    with open(manifest_path, newline="") as handle:
        rows = list(csv.DictReader(handle))
    if not rows:
        raise ValueError(f"Manifest is empty: {manifest_path}")

    lookup = load_metadata_lookup()
    enriched = []
    for row in rows:
        key = (row.get("image_type", ""), row.get("filename", ""))
        meta = lookup.get(key, {})
        merged = dict(meta)
        merged.update({k: v for k, v in row.items() if v not in (None, "")})
        enriched.append(merged)

    missing = set()
    for required in REQUIRED_MANIFEST_COLUMNS:
        if any(not row.get(required) for row in enriched):
            missing.add(required)
    if missing:
        cols = ", ".join(sorted(missing))
        raise ValueError(f"Manifest rows still missing required fields after metadata join: {cols}")
    return enriched


def resolve_image_path(row: dict[str, str]) -> Path:
    return DATA_ROOT / row["image_type"] / row["filename"]


def manifest_summary(rows: list[dict[str, str]]) -> dict[str, Any]:
    image_types = Counter(row.get("image_type", "") for row in rows)
    preservation = Counter(row.get("preservation", "") for row in rows)
    species = Counter(row.get("species", "") for row in rows)
    mounting = Counter(row.get("mounting", "") for row in rows)
    return {
        "image_count": len(rows),
        "image_types": dict(sorted(image_types.items())),
        "preservation": dict(sorted(preservation.items())),
        "mounting": dict(sorted(mounting.items())),
        "species": dict(sorted(species.items())),
    }


def rows_missing_images(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    missing = []
    for row in rows:
        image_path = resolve_image_path(row)
        if not image_path.exists():
            missing.append(
                {
                    "filename": row["filename"],
                    "image_type": row["image_type"],
                    "species": row.get("species", ""),
                    "expected_path": str(image_path),
                }
            )
    return missing


def current_git_hash() -> str:
    try:
        return subprocess.check_output(
            ["git", "rev-parse", "--short", "HEAD"],
            cwd=PROJECT,
            text=True,
            timeout=5,
        ).strip()
    except Exception:
        return ""


def copy_manifest_into_run(manifest_path: Path, paths: CanonicalRunPaths) -> Path:
    dest = paths.manifests / Path(manifest_path).name
    dest.write_text(Path(manifest_path).read_text())
    return dest


def build_run_config(
    run_tag: str,
    manifest_path: Path,
    rows: list[dict[str, str]],
    paths: CanonicalRunPaths,
    extra: dict[str, Any] | None = None,
) -> dict[str, Any]:
    config = {
        "run_tag": run_tag,
        "created_at": time.strftime("%Y-%m-%dT%H:%M:%S"),
        "git_hash": current_git_hash(),
        "repo_root": str(PROJECT),
        "manifest_path": str(Path(manifest_path).resolve()),
        "run_paths": paths.as_dict(),
        "manifest_summary": manifest_summary(rows),
        "shared_inputs": {
            "data_root": str(DATA_ROOT),
            "metadata_root": str(DATA_ROOT / "metadata"),
            "manifests_root": str(DATA_ROOT / "manifests"),
        },
        "pipeline_entrypoints": {
            "prepare": str(PROJECT / "scripts" / "prepare_canonical_run.py"),
            "preflight": str(PROJECT / "scripts" / "preflight_three_pipelines.py"),
            "run_three_pipelines": str(PROJECT / "scripts" / "run_three_pipelines.py"),
            "cell_size": str(PROJECT / "cell_size_segmentation_pipeline" / "run_from_manifest.py"),
            "nucleus_iod": str(PROJECT / "nucleus_iod_estimate_pipeline" / "run_from_manifest.py"),
        },
        "traceability_contract": [
            "result_row",
            "object_row",
            "mask_with_coordinates",
            "tile_or_block_coordinates",
            "source_image_file",
        ],
    }
    if extra:
        config.update(extra)
    return config


def write_json(path: Path, data: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=2, sort_keys=False))


def write_run_readme(paths: CanonicalRunPaths, config: dict[str, Any]) -> Path:
    summary = config["manifest_summary"]
    lines = [
        f"# Canonical Run: {config['run_tag']}",
        "",
        "This run root is the shared execution/output surface for the active production pipelines.",
        "",
        "## Manifest",
        "",
        f"- Source manifest: `{config['manifest_path']}`",
        f"- Copied manifest: `{paths.manifests / Path(config['manifest_path']).name}`",
        "",
        "## Summary",
        "",
        f"- Images: {summary['image_count']}",
        f"- Image types: {summary['image_types']}",
        f"- Preservation: {summary['preservation']}",
        "",
        "## Pipeline Roots",
        "",
        f"- Cell size: `{paths.cell_size_segmentation}`",
        f"- Nucleus IOD: `{paths.nucleus_iod}`",
        f"- QC: `{paths.qc}`",
        f"- Traceability: `{paths.traceability}`",
        f"- Logs: `{paths.logs}`",
    ]
    readme = paths.root / "README.md"
    readme.write_text("\n".join(lines) + "\n")
    return readme
