#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
import subprocess
import sys
import time
import traceback
from concurrent.futures import ProcessPoolExecutor, as_completed, wait, FIRST_COMPLETED
from pathlib import Path

import numpy as np

PROJECT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT / "src"))
sys.path.insert(0, str(PROJECT / "scripts"))

from cellprofiler_tools.convergence import (  # noqa: E402
    ConvergenceTracker,
    interleave_by_species,
)
from cellprofiler_tools.pipeline_runs import (  # noqa: E402
    current_git_hash,
    load_manifest_rows,
    resolve_image_path,
)
def _load_backend(name: str):
    if name == "python":
        from nucleus_iod_python import CSV_COLUMNS, load_background_cache, process_image, save_csv
    elif name == "imagej":
        from imagej_nucleus_iod import CSV_COLUMNS, load_background_cache, process_image, save_csv
    else:
        raise ValueError(f"Unknown backend: {name}")
    return CSV_COLUMNS, load_background_cache, process_image, save_csv


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--artifact-dir", type=Path)
    parser.add_argument("--qc-dir", type=Path)
    parser.add_argument("--trace-dir", type=Path)
    parser.add_argument("--background-cache", type=Path)
    parser.add_argument("--tile-filter", choices=["green", "auto", "all"], default="auto")
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--run-qc", action="store_true")
    parser.add_argument("--write-nucleus-qc", action="store_true")
    parser.add_argument("--primary-tier", default="strict_core")
    parser.add_argument("--segmentation-dir", type=Path, default=PROJECT / "output" / "segmentation")
    # Convergence control
    parser.add_argument("--min-cells", type=int, default=30,
                        help="Min nuclei per species before checking convergence")
    parser.add_argument("--max-cells-species", type=int, default=500,
                        help="Hard cap on nuclei per species")
    parser.add_argument("--sem-threshold", type=float, default=3.0,
                        help="SEM%% threshold for convergence")
    parser.add_argument("--convergence-key", default="iod",
                        help="Measurement key to track for convergence (default: iod)")
    parser.add_argument("--no-convergence", action="store_true",
                        help="Disable convergence-based early stopping")
    parser.add_argument("--image-type", choices=["brightfield", "pmount", "all"],
                        default="all",
                        help="Process only this image type (default: all)")
    parser.add_argument("--workers", type=int, default=1,
                        help="Number of parallel image workers (default: 1)")
    parser.add_argument("--backend", choices=["python", "imagej"], default="python",
                        help="Segmentation backend (default: python)")
    parser.add_argument("--threshold-csv", type=Path,
                        help="Per-image threshold CSV (columns: filename, threshold, excluded)")
    return parser.parse_args()


def _process_one(image_path, tile_filter, cached_i_bg, artifact_dir, backend="python",
                 manual_threshold=None):
    """Worker function for ProcessPoolExecutor."""
    _, _, process_image_fn, _ = _load_backend(backend)
    return process_image_fn(
        image_path,
        tile_filter=tile_filter,
        cached_i_bg=cached_i_bg,
        artifact_dir=artifact_dir,
        manual_threshold=manual_threshold,
    )


def run_command(cmd: list[str]) -> None:
    result = subprocess.run(cmd, cwd=PROJECT, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"Command failed ({result.returncode}): {' '.join(cmd)}")


def main() -> None:
    args = parse_args()
    CSV_COLUMNS, load_background_cache, process_image, save_csv = _load_backend(args.backend)
    print(f"Nucleus IOD backend: {args.backend}")

    # Load per-image thresholds and exclusions
    threshold_map = {}  # filename -> int threshold
    excluded_set = set()
    if args.threshold_csv and args.threshold_csv.exists():
        with open(args.threshold_csv) as f:
            for tr in csv.DictReader(f):
                fn = tr.get("filename", "").strip()
                if not fn:
                    continue
                if tr.get("excluded", "").strip() == "true":
                    excluded_set.add(fn)
                else:
                    try:
                        threshold_map[fn] = int(tr["threshold"])
                    except (KeyError, ValueError):
                        pass
        print(f"Thresholds: {len(threshold_map)} per-image, {len(excluded_set)} excluded")

    rows = load_manifest_rows(args.manifest)
    if args.image_type != "all":
        rows = [r for r in rows if r["image_type"] == args.image_type]
        print(f"Nucleus IOD: filtered to {len(rows)} {args.image_type} images")
    if excluded_set:
        before = len(rows)
        rows = [r for r in rows if r["filename"] not in excluded_set]
        print(f"Nucleus IOD: excluded {before - len(rows)} images, {len(rows)} remaining")
    args.output_dir.mkdir(parents=True, exist_ok=True)
    measurement_dir = args.output_dir / "measurements"
    measurement_dir.mkdir(parents=True, exist_ok=True)
    measurements_csv = measurement_dir / "nucleus_iod_measurements.csv"
    artifact_dir = args.artifact_dir or (args.output_dir / "artifacts")
    qc_dir = args.qc_dir
    trace_dir = args.trace_dir
    bg_cache = load_background_cache(args.background_cache) if args.background_cache else {}
    index_rows = []

    # Convergence tracker
    tracker = None
    if not args.no_convergence and not args.dry_run:
        tracker = ConvergenceTracker(
            min_cells=args.min_cells,
            max_cells=args.max_cells_species,
            sem_pct=args.sem_threshold,
        )
        print(f"IOD convergence: min={args.min_cells}, max={args.max_cells_species}/species, "
              f"SEM%<{args.sem_threshold}%, key={args.convergence_key}")

    # Interleave rows by species for even sampling
    rows = interleave_by_species(rows)
    species_set = {r["species"] for r in rows if r.get("species")}
    print(f"Nucleus IOD: {len(rows)} images, {len(species_set)} species")

    completed = set()
    if args.resume and measurements_csv.exists():
        with open(measurements_csv, newline="") as handle:
            prior = list(csv.DictReader(handle))
            completed = {row["filename"] for row in prior}
        if tracker:
            for row in prior:
                sp = row.get("species", "")
                val = row.get(args.convergence_key)
                if sp and val not in (None, ""):
                    try:
                        tracker.species_values[sp].append(float(val))
                    except ValueError:
                        pass
            print(f"Seeded IOD tracker: {sum(len(v) for v in tracker.species_values.values())} "
                  f"nuclei across {len(tracker.species_values)} species")

    t_start = time.time()
    total_nuclei = 0
    skipped_converged = 0
    errors = []
    all_converged = False

    def _make_index_row(row, status, image_path, image_artifact_dir):
        return {
            "filename": row["filename"], "image_type": row["image_type"],
            "species": row.get("species", ""), "slide_id": row["slide_id"],
            "specimen_id": row["specimen_id"], "status": status,
            "source_image_path": str(image_path),
            "artifact_dir": str(image_artifact_dir),
            "tile_manifest_path": str(image_artifact_dir / "tile_manifest.csv"),
            "raw_imagej_results_path": str(image_artifact_dir / "raw_imagej_results.csv"),
        }

    def _handle_result(nucleus_rows, row, species, t0):
        nonlocal total_nuclei
        save_csv(
            nucleus_rows,
            measurements_csv,
            metadata={
                "slide_id": row["slide_id"],
                "specimen_id": row["specimen_id"],
                "species": species,
                "image_type": row["image_type"],
            },
        )
        total_nuclei += len(nucleus_rows)
        dt = time.time() - t0
        print(f"  {row['filename']}: {len(nucleus_rows)} nuclei [{dt:.0f}s] (total={total_nuclei})")
        if tracker and species and nucleus_rows:
            values = []
            for nr in nucleus_rows:
                val = nr.get(args.convergence_key)
                if val not in (None, ""):
                    try:
                        values.append(float(val))
                    except ValueError:
                        pass
            if values:
                tracker.add(species, values)
                print(f"  {species}: {tracker.status(species)}")

    def _should_skip(row, i):
        nonlocal skipped_converged
        species = row.get("species", "")
        image_path = resolve_image_path(row)
        image_artifact_dir = artifact_dir / Path(row["filename"]).stem

        if args.resume and row["filename"] in completed:
            index_rows.append(_make_index_row(row, "skipped_existing", image_path, image_artifact_dir))
            return True
        if tracker and species and tracker.is_done(species):
            print(f"  SKIP [{i+1}/{len(rows)}] {row['filename']}: "
                  f"{species} converged ({tracker.status(species)})")
            skipped_converged += 1
            index_rows.append(_make_index_row(row, "skipped_converged", image_path, image_artifact_dir))
            return True
        return False

    if args.workers > 1 and not args.dry_run:
        print(f"Parallel mode: {args.workers} workers")
        pool = ProcessPoolExecutor(max_workers=args.workers)
        # future -> (row, species, submit_time)
        pending: dict = {}

        def _drain_completed():
            """Drain all completed futures, update tracker and CSV."""
            nonlocal all_converged
            if not pending:
                return
            done_futures = [f for f in pending if f.done()]
            for future in done_futures:
                row, species, t0 = pending.pop(future)
                image_path = resolve_image_path(row)
                image_artifact_dir = artifact_dir / Path(row["filename"]).stem
                try:
                    nucleus_rows = future.result()
                    _handle_result(nucleus_rows, row, species, t0)
                    index_rows.append(_make_index_row(row, "completed", image_path, image_artifact_dir))
                except Exception as e:
                    dt = time.time() - t0
                    errors.append((row["filename"], str(e)))
                    print(f"  ERROR {row['filename']} [{dt:.0f}s]: {e}")
                    index_rows.append(_make_index_row(row, "error", image_path, image_artifact_dir))
            if tracker and species_set and all(tracker.is_done(sp) for sp in species_set):
                all_converged = True

        for i, row in enumerate(rows):
            if all_converged:
                species = row.get("species", "")
                image_path = resolve_image_path(row)
                image_artifact_dir = artifact_dir / Path(row["filename"]).stem
                index_rows.append(_make_index_row(row, "skipped_converged", image_path, image_artifact_dir))
                skipped_converged += 1
                continue

            if _should_skip(row, i):
                continue

            species = row.get("species", "")
            image_path = resolve_image_path(row)
            image_artifact_dir = artifact_dir / Path(row["filename"]).stem
            print(f"[{i+1}/{len(rows)}] {row['filename']} ({species}, {row['image_type']})")

            # If pool is full, wait for one to finish
            while len(pending) >= args.workers:
                done, _ = wait(pending.keys(), return_when=FIRST_COMPLETED)
                _drain_completed()

            future = pool.submit(
                _process_one,
                str(image_path),
                args.tile_filter,
                bg_cache.get(row["filename"]),
                str(image_artifact_dir),
                args.backend,
                threshold_map.get(row["filename"]),
            )
            pending[future] = (row, species, time.time())

        # Drain remaining
        if pending:
            wait(pending.keys())
            _drain_completed()
        pool.shutdown(wait=False)

        if all_converged:
            print(f"\nAll {len(species_set)} species converged! Stopping early.")
            if tracker:
                print(tracker.summary())

    else:
        for i, row in enumerate(rows):
            if _should_skip(row, i):
                continue

            species = row.get("species", "")
            image_path = resolve_image_path(row)
            image_artifact_dir = artifact_dir / Path(row["filename"]).stem

            if not args.dry_run:
                print(f"[{i+1}/{len(rows)}] {row['filename']} ({species}, {row['image_type']})")
                t0 = time.time()
                try:
                    nucleus_rows = process_image(
                        image_path,
                        tile_filter=args.tile_filter,
                        cached_i_bg=bg_cache.get(row["filename"]),
                        artifact_dir=image_artifact_dir,
                        manual_threshold=threshold_map.get(row["filename"]),
                    )
                    _handle_result(nucleus_rows, row, species, t0)
                except Exception as e:
                    dt = time.time() - t0
                    errors.append((row["filename"], str(e)))
                    print(f"  ERROR [{dt:.0f}s]: {e}")
                    traceback.print_exc()

            index_rows.append(_make_index_row(
                row,
                "planned" if args.dry_run else "completed",
                image_path, image_artifact_dir,
            ))

            # Check if ALL species converged
            if tracker and species_set:
                if all(tracker.is_done(sp) for sp in species_set):
                    print(f"\nAll {len(species_set)} species converged! Stopping early.")
                    print(tracker.summary())
                    break

    elapsed = time.time() - t_start

    index_path = args.output_dir / "image_index.csv"
    with open(index_path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=index_rows[0].keys() if index_rows else [
            "filename", "image_type", "species", "slide_id", "specimen_id", "status",
            "source_image_path", "artifact_dir", "tile_manifest_path", "raw_imagej_results_path",
        ])
        writer.writeheader()
        writer.writerows(index_rows)

    manifest = {
        "timestamp": time.strftime("%Y-%m-%dT%H:%M:%S"),
        "git_hash": current_git_hash(),
        "manifest_path": str(args.manifest.resolve()),
        "output_dir": str(args.output_dir),
        "measurements_csv": str(measurements_csv),
        "artifact_dir": str(artifact_dir),
        "image_count": len(rows),
        "tile_filter": args.tile_filter,
        "background_cache": str(args.background_cache) if args.background_cache else "",
        "elapsed_seconds": round(elapsed, 1),
        "total_nuclei": total_nuclei,
        "images_skipped_converged": skipped_converged,
        "convergence_enabled": not args.no_convergence,
        "convergence_key": args.convergence_key,
        "errors": [{"filename": fn, "error": err} for fn, err in errors],
        "dry_run": args.dry_run,
    }
    (args.output_dir / "run_manifest.json").write_text(json.dumps(manifest, indent=2))
    print(f"Nucleus IOD image index: {index_path}")
    print(f"Nucleus IOD manifest: {args.output_dir / 'run_manifest.json'}")

    if tracker and tracker.species_values:
        print(f"\nConvergence summary ({args.convergence_key}):")
        print(tracker.summary())

        # Save convergence summary CSV
        conv_rows = []
        for sp in sorted(tracker.species_values):
            vals = np.array(tracker.species_values[sp])
            n = len(vals)
            sem = vals.std(ddof=1) / np.sqrt(n) if n > 1 else float("inf")
            sem_pct = 100 * sem / vals.mean() if vals.mean() > 0 else float("inf")
            conv_rows.append({
                "species": sp,
                "n_nuclei": n,
                f"mean_{args.convergence_key}": round(float(vals.mean()), 3),
                f"sd_{args.convergence_key}": round(float(vals.std(ddof=1)), 3) if n > 1 else 0,
                "sem_pct": round(float(sem_pct), 2) if n > 1 else 0,
                "converged": tracker.truly_converged(sp),
                "capped_not_converged": tracker.is_capped(sp),
            })
        conv_path = args.output_dir / "convergence_summary.csv"
        with open(conv_path, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=conv_rows[0].keys())
            w.writeheader()
            w.writerows(conv_rows)
        print(f"Convergence CSV: {conv_path}")

    if skipped_converged:
        print(f"Skipped {skipped_converged} images (species already converged)")
    if errors:
        print(f"\n{len(errors)} errors:")
        for fn, err in errors:
            print(f"  {fn}: {err}")

    if not args.run_qc or args.dry_run:
        return

    if qc_dir is None:
        raise ValueError("--run-qc requires --qc-dir")
    if trace_dir is None:
        raise ValueError("--run-qc requires --trace-dir")

    qc_dir.mkdir(parents=True, exist_ok=True)
    trace_dir.mkdir(parents=True, exist_ok=True)
    qc_trace_link = qc_dir / "traceability"
    if qc_trace_link.exists() or qc_trace_link.is_symlink():
        if qc_trace_link.resolve() != trace_dir.resolve():
            qc_trace_link.unlink()
            qc_trace_link.symlink_to(trace_dir)
    else:
        qc_trace_link.symlink_to(trace_dir)

    qc_cmd = [
        sys.executable,
        str(PROJECT / "scripts" / "qc_analysis.py"),
        "--nucleus-iod",
        str(measurements_csv),
        "--output-dir",
        str(qc_dir),
    ]
    if args.write_nucleus_qc:
        qc_cmd.append("--write-nucleus-qc")
    run_command(qc_cmd)
    run_command(
        [
            sys.executable,
            str(PROJECT / "scripts" / "build_final_species_results.py"),
            "--qc-dir",
            str(qc_dir),
            "--primary-tier",
            args.primary_tier,
        ]
    )
    run_command(
        [
            sys.executable,
            str(PROJECT / "scripts" / "build_traceability_bundle.py"),
            "--qc-dir",
            str(qc_dir),
            "--nucleus-iod",
            str(measurements_csv),
            "--artifact-root",
            str(artifact_dir),
            "--segmentation-dir",
            str(args.segmentation_dir),
            "--output-dir",
            str(trace_dir),
        ]
    )
    run_command(
        [
            sys.executable,
            str(PROJECT / "scripts" / "build_traceability_viewer.py"),
            "--trace-dir",
            str(qc_dir / "traceability"),
        ]
    )
    run_command(
        [
            sys.executable,
            str(PROJECT / "scripts" / "build_interim_estimate_bundle.py"),
            "--qc-dir",
            str(qc_dir),
        ]
    )
    run_command(
        [
            sys.executable,
            str(PROJECT / "scripts" / "build_analysis_landing_page.py"),
            "--qc-dir",
            str(qc_dir),
        ]
    )


if __name__ == "__main__":
    main()
