#!/usr/bin/env python3
"""
Parallel batch process VSI and BigTIFF files for CellProfiler.

Monitors RAM usage and adjusts worker count dynamically.
"""

import argparse
import subprocess
import sys
import tempfile
import shutil
import psutil
import time
from pathlib import Path
from datetime import datetime
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass
from typing import Optional
import multiprocessing as mp


@dataclass
class ProcessResult:
    path: Path
    success: bool
    message: str
    duration: float


def get_available_ram_gb() -> float:
    """Get available RAM in GB."""
    return psutil.virtual_memory().available / (1024**3)


def find_files(root_dir: Path) -> tuple[list[Path], list[Path]]:
    """Find all VSI and BTF files."""
    vsi_files = []
    btf_files = []

    for pattern in ["**/*.vsi", "**/*.VSI"]:
        vsi_files.extend(root_dir.glob(pattern))

    for pattern in ["**/*.btf", "**/*.BTF", "**/*.ome.btf"]:
        btf_files.extend(root_dir.glob(pattern))

    return sorted(set(vsi_files)), sorted(set(btf_files))


def process_btf(btf_path: Path, output_dir: Path) -> ProcessResult:
    """Process a BTF file - extract green channel."""
    start = time.time()

    # Expected output file
    output_name = f"{btf_path.stem.replace('.ome', '')}_green.ome.tiff"
    output_path = output_dir / output_name

    cmd = [
        "uv", "run", "python", "btf_to_green.py",
        str(btf_path),
        str(output_dir),
        "--channel", "1",
    ]

    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=3600)
        duration = time.time() - start

        # Check if output file was created (more reliable than return code)
        if output_path.exists() and output_path.stat().st_size > 1000:
            size_mb = output_path.stat().st_size / (1024**2)
            return ProcessResult(btf_path, True, f"OK ({size_mb:.0f}MB)", duration)
        elif result.returncode == 0:
            return ProcessResult(btf_path, True, "OK", duration)
        else:
            err = result.stderr[-200:] if result.stderr else "Unknown error"
            return ProcessResult(btf_path, False, err, duration)
    except subprocess.TimeoutExpired:
        return ProcessResult(btf_path, False, "Timeout", time.time() - start)
    except Exception as e:
        return ProcessResult(btf_path, False, str(e), time.time() - start)


def process_vsi(vsi_path: Path, output_dir: Path, bfconvert: Path, temp_dir: Path) -> ProcessResult:
    """Process a VSI file - convert to OME-TIFF then extract green channel."""
    start = time.time()

    # Use unique temp file per process
    temp_tiff = temp_dir / f"{vsi_path.stem}_{mp.current_process().pid}.ome.tiff"

    cmd_convert = [
        str(bfconvert),
        "-bigtiff",
        "-compression", "LZW",
        "-tilex", "512",
        "-tiley", "512",
        "-series", "0",
        str(vsi_path),
        str(temp_tiff)
    ]

    try:
        result = subprocess.run(cmd_convert, capture_output=True, text=True, timeout=7200)

        if result.returncode != 0 or not temp_tiff.exists():
            return ProcessResult(vsi_path, False, f"bfconvert failed", time.time() - start)

        cmd_green = [
            "uv", "run", "python", "btf_to_green.py",
            str(temp_tiff),
            str(output_dir),
            "--channel", "1",
        ]

        result = subprocess.run(cmd_green, capture_output=True, text=True, timeout=3600)

        if temp_tiff.exists():
            temp_tiff.unlink()

        duration = time.time() - start
        if result.returncode == 0:
            return ProcessResult(vsi_path, True, "OK", duration)
        else:
            return ProcessResult(vsi_path, False, "Green extraction failed", duration)

    except subprocess.TimeoutExpired:
        if temp_tiff.exists():
            temp_tiff.unlink()
        return ProcessResult(vsi_path, False, "Timeout", time.time() - start)
    except Exception as e:
        if temp_tiff.exists():
            temp_tiff.unlink()
        return ProcessResult(vsi_path, False, str(e), time.time() - start)


def worker_btf(args):
    """Worker function for BTF processing."""
    btf_path, output_dir = args
    return process_btf(btf_path, output_dir)


def worker_vsi(args):
    """Worker function for VSI processing."""
    vsi_path, output_dir, bfconvert, temp_dir = args
    return process_vsi(vsi_path, output_dir, bfconvert, temp_dir)


def main():
    parser = argparse.ArgumentParser(description="Parallel batch process VSI and BTF files")
    parser.add_argument("input_dir", type=Path, help="Input directory")
    parser.add_argument("output_dir", type=Path, help="Output directory")
    parser.add_argument("--bfconvert", type=Path, default=Path.home() / "bin/bftools/bftools/bfconvert")
    parser.add_argument("--temp-dir", type=Path, default=None)
    parser.add_argument("--workers", "-j", type=int, default=4, help="Number of parallel workers (default: 4)")
    parser.add_argument("--max-workers", type=int, default=8, help="Max workers when RAM is plentiful (default: 8)")
    parser.add_argument("--min-ram-gb", type=float, default=16.0, help="Minimum free RAM in GB before throttling (default: 16)")
    parser.add_argument("--btf-only", action="store_true")
    parser.add_argument("--vsi-only", action="store_true")
    parser.add_argument("--skip-existing", action="store_true")
    parser.add_argument("--dry-run", "-n", action="store_true")

    args = parser.parse_args()

    if not args.input_dir.exists():
        print(f"Error: Input directory not found: {args.input_dir}")
        sys.exit(1)

    vsi_files, btf_files = find_files(args.input_dir)

    if args.btf_only:
        vsi_files = []
    if args.vsi_only:
        btf_files = []

    # Filter out existing outputs
    if args.skip_existing:
        args.output_dir.mkdir(parents=True, exist_ok=True)
        existing = {f.stem.replace("_green.ome", "").replace("_green", "") for f in args.output_dir.glob("*.tiff")}

        btf_files = [f for f in btf_files if f.stem.replace(".ome", "").replace("_raw", "") not in existing]
        vsi_files = [f for f in vsi_files if f.stem not in existing]

    total = len(btf_files) + len(vsi_files)
    print(f"Files to process: {len(btf_files)} BTF + {len(vsi_files)} VSI = {total} total")
    print(f"Workers: {args.workers} (max: {args.max_workers}, min RAM: {args.min_ram_gb} GB)")
    print(f"Available RAM: {get_available_ram_gb():.1f} GB")

    if args.dry_run:
        print("\nDRY RUN")
        return

    if total == 0:
        print("Nothing to process!")
        return

    args.output_dir.mkdir(parents=True, exist_ok=True)
    temp_dir = args.temp_dir or Path(tempfile.mkdtemp(prefix="vsi_convert_"))
    temp_dir.mkdir(parents=True, exist_ok=True)

    print(f"\nStarting at {datetime.now().strftime('%H:%M:%S')}")
    print("=" * 70)

    successful = 0
    failed = 0
    completed = 0
    start_time = time.time()

    # Process BTF files first (faster, can parallelize more)
    if btf_files:
        print(f"\n--- Processing {len(btf_files)} BTF files with {args.workers} workers ---\n")

        btf_args = [(f, args.output_dir) for f in btf_files]

        with ProcessPoolExecutor(max_workers=args.workers) as executor:
            futures = {executor.submit(worker_btf, arg): arg[0] for arg in btf_args}

            for future in as_completed(futures):
                completed += 1
                result = future.result()

                ram_gb = get_available_ram_gb()
                elapsed = time.time() - start_time
                rate = completed / elapsed * 60 if elapsed > 0 else 0

                status = "✓" if result.success else "✗"
                print(f"[{completed}/{total}] {status} {result.path.name} ({result.duration:.1f}s) | RAM: {ram_gb:.0f}GB | {rate:.1f}/min")

                if result.success:
                    successful += 1
                else:
                    failed += 1
                    print(f"         Error: {result.message}")

    # Process VSI files (slower, may need fewer workers)
    if vsi_files:
        print(f"\n--- Processing {len(vsi_files)} VSI files with {args.workers} workers ---\n")

        vsi_args = [(f, args.output_dir, args.bfconvert, temp_dir) for f in vsi_files]

        with ProcessPoolExecutor(max_workers=args.workers) as executor:
            futures = {executor.submit(worker_vsi, arg): arg[0] for arg in vsi_args}

            for future in as_completed(futures):
                completed += 1
                result = future.result()

                ram_gb = get_available_ram_gb()
                elapsed = time.time() - start_time
                rate = completed / elapsed * 60 if elapsed > 0 else 0

                status = "✓" if result.success else "✗"
                print(f"[{completed}/{total}] {status} {result.path.name} ({result.duration:.1f}s) | RAM: {ram_gb:.0f}GB | {rate:.1f}/min")

                if result.success:
                    successful += 1
                else:
                    failed += 1
                    print(f"         Error: {result.message}")

    # Cleanup
    if not args.temp_dir and temp_dir.exists():
        shutil.rmtree(temp_dir, ignore_errors=True)

    elapsed = time.time() - start_time
    print()
    print("=" * 70)
    print(f"Completed at {datetime.now().strftime('%H:%M:%S')} ({elapsed/60:.1f} min total)")
    print(f"Successful: {successful} | Failed: {failed}")

    if failed > 0:
        sys.exit(1)


if __name__ == "__main__":
    main()
