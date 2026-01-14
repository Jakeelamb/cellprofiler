#!/usr/bin/env python3
"""
Batch process VSI and BigTIFF files for CellProfiler.

Workflow:
1. BTF files -> Extract green channel -> Output green OME-TIFF
2. VSI files -> Convert to OME-TIFF (via bfconvert) -> Extract green channel -> Output green OME-TIFF

Usage:
    uv run python batch_process.py /path/to/easystore /output/dir
"""

import argparse
import subprocess
import sys
import tempfile
import shutil
from pathlib import Path
from datetime import datetime


def find_files(root_dir: Path) -> tuple[list[Path], list[Path]]:
    """Find all VSI and BTF files."""
    vsi_files = []
    btf_files = []

    for pattern in ["**/*.vsi", "**/*.VSI"]:
        vsi_files.extend(root_dir.glob(pattern))

    for pattern in ["**/*.btf", "**/*.BTF", "**/*.ome.btf"]:
        btf_files.extend(root_dir.glob(pattern))

    return sorted(set(vsi_files)), sorted(set(btf_files))


def process_btf(btf_path: Path, output_dir: Path, verbose: bool = False) -> tuple[bool, str]:
    """Process a BTF file - extract green channel."""
    cmd = [
        "uv", "run", "python", "btf_to_green.py",
        str(btf_path),
        str(output_dir),
        "--channel", "1",  # Green
    ]
    if verbose:
        cmd.append("-v")

    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=3600)
        if result.returncode == 0:
            return True, "Green channel extracted"
        else:
            return False, result.stderr[-500:] if result.stderr else "Unknown error"
    except subprocess.TimeoutExpired:
        return False, "Timeout after 1 hour"
    except Exception as e:
        return False, str(e)


def process_vsi(vsi_path: Path, output_dir: Path, bfconvert: Path, temp_dir: Path, verbose: bool = False) -> tuple[bool, str]:
    """Process a VSI file - convert to OME-TIFF then extract green channel."""

    # Step 1: Convert VSI to OME-TIFF using bfconvert
    temp_tiff = temp_dir / (vsi_path.stem + ".ome.tiff")

    cmd_convert = [
        str(bfconvert),
        "-bigtiff",
        "-compression", "LZW",
        "-tilex", "512",
        "-tiley", "512",
        "-series", "0",  # Full resolution only
        str(vsi_path),
        str(temp_tiff)
    ]

    try:
        if verbose:
            print(f"    Converting VSI to OME-TIFF...")
        result = subprocess.run(cmd_convert, capture_output=True, text=True, timeout=7200)

        if result.returncode != 0 or not temp_tiff.exists():
            return False, f"bfconvert failed: {result.stderr[-500:] if result.stderr else 'Unknown error'}"

        # Step 2: Extract green channel from the temp TIFF
        if verbose:
            print(f"    Extracting green channel...")

        cmd_green = [
            "uv", "run", "python", "btf_to_green.py",
            str(temp_tiff),
            str(output_dir),
            "--channel", "1",
        ]
        if verbose:
            cmd_green.append("-v")

        result = subprocess.run(cmd_green, capture_output=True, text=True, timeout=3600)

        # Clean up temp file
        if temp_tiff.exists():
            temp_tiff.unlink()

        if result.returncode == 0:
            return True, "Converted and green channel extracted"
        else:
            return False, result.stderr[-500:] if result.stderr else "Green extraction failed"

    except subprocess.TimeoutExpired:
        if temp_tiff.exists():
            temp_tiff.unlink()
        return False, "Timeout"
    except Exception as e:
        if temp_tiff.exists():
            temp_tiff.unlink()
        return False, str(e)


def main():
    parser = argparse.ArgumentParser(description="Batch process VSI and BTF files for CellProfiler")
    parser.add_argument("input_dir", type=Path, help="Input directory containing VSI/BTF files")
    parser.add_argument("output_dir", type=Path, help="Output directory for green channel TIFFs")
    parser.add_argument("--bfconvert", type=Path, default=Path.home() / "bin/bftools/bftools/bfconvert",
                        help="Path to bfconvert")
    parser.add_argument("--temp-dir", type=Path, default=None, help="Temp directory for intermediate files")
    parser.add_argument("--btf-only", action="store_true", help="Only process BTF files")
    parser.add_argument("--vsi-only", action="store_true", help="Only process VSI files")
    parser.add_argument("--dry-run", "-n", action="store_true", help="Show what would be done")
    parser.add_argument("--verbose", "-v", action="store_true", help="Verbose output")
    parser.add_argument("--skip-existing", action="store_true", help="Skip files that already have output")

    args = parser.parse_args()

    # Validate inputs
    if not args.input_dir.exists():
        print(f"Error: Input directory not found: {args.input_dir}")
        sys.exit(1)

    if not args.bfconvert.exists() and not args.btf_only:
        print(f"Error: bfconvert not found at {args.bfconvert}")
        print("Install Bio-Formats tools or use --btf-only")
        sys.exit(1)

    # Find files
    vsi_files, btf_files = find_files(args.input_dir)

    if args.btf_only:
        vsi_files = []
    if args.vsi_only:
        btf_files = []

    total_files = len(btf_files) + len(vsi_files)
    print(f"Found {len(btf_files)} BTF files and {len(vsi_files)} VSI files")
    print(f"Total: {total_files} files to process")
    print(f"Output directory: {args.output_dir}")

    if args.dry_run:
        print("\nDRY RUN - no files will be processed")
        print("\nBTF files:")
        for f in btf_files[:10]:
            print(f"  {f}")
        if len(btf_files) > 10:
            print(f"  ... and {len(btf_files) - 10} more")
        print("\nVSI files:")
        for f in vsi_files[:10]:
            print(f"  {f}")
        if len(vsi_files) > 10:
            print(f"  ... and {len(vsi_files) - 10} more")
        sys.exit(0)

    # Create output directory
    args.output_dir.mkdir(parents=True, exist_ok=True)

    # Create temp directory for VSI conversion
    temp_dir = args.temp_dir or Path(tempfile.mkdtemp(prefix="vsi_convert_"))
    temp_dir.mkdir(parents=True, exist_ok=True)

    print(f"\nStarting processing at {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("=" * 60)

    successful = 0
    failed = 0
    skipped = 0
    current = 0

    # Process BTF files first (faster)
    if btf_files:
        print(f"\n--- Processing {len(btf_files)} BTF files ---\n")

        for btf_file in btf_files:
            current += 1

            # Check if output exists
            output_name = f"{btf_file.stem.replace('.ome', '')}_green.ome.tiff"
            output_path = args.output_dir / output_name

            if args.skip_existing and output_path.exists():
                skipped += 1
                print(f"[{current}/{total_files}] SKIP {btf_file.name} (output exists)")
                continue

            print(f"[{current}/{total_files}] {btf_file.name}")

            success, msg = process_btf(btf_file, args.output_dir, args.verbose)

            if success:
                successful += 1
                print(f"  ✓ {msg}")
            else:
                failed += 1
                print(f"  ✗ {msg}")

    # Process VSI files
    if vsi_files:
        print(f"\n--- Processing {len(vsi_files)} VSI files ---\n")

        for vsi_file in vsi_files:
            current += 1

            # Check if output exists
            output_name = f"{vsi_file.stem}_green.ome.tiff"
            output_path = args.output_dir / output_name

            if args.skip_existing and output_path.exists():
                skipped += 1
                print(f"[{current}/{total_files}] SKIP {vsi_file.name} (output exists)")
                continue

            print(f"[{current}/{total_files}] {vsi_file.name}")

            success, msg = process_vsi(vsi_file, args.output_dir, args.bfconvert, temp_dir, args.verbose)

            if success:
                successful += 1
                print(f"  ✓ {msg}")
            else:
                failed += 1
                print(f"  ✗ {msg}")

    # Cleanup temp directory
    if not args.temp_dir and temp_dir.exists():
        shutil.rmtree(temp_dir, ignore_errors=True)

    print()
    print("=" * 60)
    print(f"Completed at {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"Successful: {successful}")
    print(f"Failed: {failed}")
    print(f"Skipped: {skipped}")

    if failed > 0:
        sys.exit(1)


if __name__ == "__main__":
    main()
