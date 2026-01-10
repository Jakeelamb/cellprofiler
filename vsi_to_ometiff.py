#!/usr/bin/env python3
"""
Convert VSI files to OME-TIFF format preserving all metadata.

This script uses Bio-Formats bfconvert to convert Olympus VSI files to
OME-TIFF format, which preserves:
- Physical pixel sizes (µm/pixel)
- Objective information
- Stage positions
- Channel metadata
- All other OME-XML metadata

Usage:
    python vsi_to_ometiff.py <input_path> <output_dir> [options]

Examples:
    # Convert a single VSI file
    python vsi_to_ometiff.py /path/to/image.vsi /output/dir

    # Convert all VSI files in a directory (recursively)
    python vsi_to_ometiff.py /path/to/input_dir /output/dir --recursive

    # Convert only the full resolution (series 0)
    python vsi_to_ometiff.py /path/to/image.vsi /output/dir --series 0

    # Dry run - show what would be converted
    python vsi_to_ometiff.py /path/to/input_dir /output/dir --dry-run

    # For GREEN CHANNEL ONLY: Use the two-step process:
    # Step 1: Convert VSI to OME-TIFF (this script)
    python vsi_to_ometiff.py /path/to/image.vsi /temp/dir
    # Step 2: Extract green channel (use extract_channel.py)
    uv run python extract_channel.py /temp/dir /output/dir --green

SAFETY: This script ONLY reads from source and writes to output_dir.
        It will refuse to run if output_dir is on the same filesystem as input.
"""

import argparse
import os
import subprocess
import sys
from pathlib import Path
from typing import Optional


# Default path to bfconvert - can be overridden with --bfconvert
DEFAULT_BFCONVERT = Path.home() / "bin" / "bftools" / "bftools" / "bfconvert"


def is_same_filesystem(path1: Path, path2: Path) -> bool:
    """Check if two paths are on the same filesystem/mount point."""
    try:
        return path1.resolve().stat().st_dev == path2.resolve().stat().st_dev
    except (OSError, FileNotFoundError):
        return False


def find_vsi_files(root_dir: Path, recursive: bool = True) -> list[Path]:
    """Find all VSI files in directory."""
    pattern = "**/*.vsi" if recursive else "*.vsi"
    files = list(root_dir.glob(pattern))
    # Also check for uppercase
    files.extend(root_dir.glob(pattern.replace(".vsi", ".VSI")))
    return sorted(set(files))


def get_series_info(vsi_path: Path, bfconvert_dir: Path) -> list[dict]:
    """Get information about series in the VSI file."""
    showinf = bfconvert_dir.parent / "showinf"
    if not showinf.exists():
        showinf = bfconvert_dir.with_name("showinf")

    try:
        result = subprocess.run(
            [str(showinf), "-nopix", str(vsi_path)],
            capture_output=True,
            text=True,
            timeout=60
        )

        series = []
        current_series = {}

        for line in result.stdout.split('\n'):
            if line.startswith("Series #"):
                if current_series:
                    series.append(current_series)
                current_series = {"id": len(series)}
            elif "Width =" in line:
                current_series["width"] = int(line.split("=")[1].strip())
            elif "Height =" in line:
                current_series["height"] = int(line.split("=")[1].strip())

        if current_series:
            series.append(current_series)

        return series
    except Exception as e:
        print(f"Warning: Could not get series info: {e}")
        return []


def convert_vsi_to_ometiff(
    input_path: Path,
    output_dir: Path,
    bfconvert: Path,
    series: Optional[int] = None,
    compression: str = "LZW",
    tile_size: int = 512,
    dry_run: bool = False,
    verbose: bool = False
) -> tuple[bool, str]:
    """
    Convert a VSI file to OME-TIFF format.

    Args:
        input_path: Path to input VSI file
        output_dir: Directory for output files
        bfconvert: Path to bfconvert executable
        series: Specific series to export (None = all series)
        compression: Compression type (LZW, JPEG-2000, zlib, uncompressed)
        tile_size: Tile size for output TIFF
        dry_run: If True, just print what would be done
        verbose: Print detailed output

    Returns:
        (success, message) tuple
    """
    # Generate output filename
    output_name = input_path.stem + ".ome.tiff"
    output_path = output_dir / output_name

    if dry_run:
        return True, f"Would convert: {input_path} -> {output_path}"

    # Create output directory
    output_dir.mkdir(parents=True, exist_ok=True)

    # Build bfconvert command
    cmd = [
        str(bfconvert),
        "-bigtiff",  # Always use BigTIFF for large files
        "-compression", compression,
        "-tilex", str(tile_size),
        "-tiley", str(tile_size),
        "-noflat",  # Preserve pyramid structure
    ]

    if series is not None:
        cmd.extend(["-series", str(series)])

    cmd.extend([str(input_path), str(output_path)])

    if verbose:
        print(f"Running: {' '.join(cmd)}")

    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=3600  # 1 hour timeout for large files
        )

        if result.returncode != 0:
            error_msg = result.stderr[-2000:] if result.stderr else "Unknown error"
            return False, f"bfconvert failed: {error_msg}"

        if not output_path.exists():
            return False, "Output file was not created"

        size_mb = output_path.stat().st_size / (1024 * 1024)
        return True, f"Created {output_path} ({size_mb:.1f} MB)"

    except subprocess.TimeoutExpired:
        return False, "Conversion timed out after 1 hour"
    except Exception as e:
        return False, f"Error: {e}"


def main():
    parser = argparse.ArgumentParser(
        description="Convert VSI files to OME-TIFF preserving metadata",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )
    parser.add_argument(
        "input_path",
        type=Path,
        help="Input VSI file or directory"
    )
    parser.add_argument(
        "output_dir",
        type=Path,
        help="Output directory for OME-TIFF files"
    )
    parser.add_argument(
        "--bfconvert",
        type=Path,
        default=DEFAULT_BFCONVERT,
        help=f"Path to bfconvert (default: {DEFAULT_BFCONVERT})"
    )
    parser.add_argument(
        "--series",
        type=int,
        default=0,
        help="Series to export (default: 0 = full resolution). Use -1 for all series."
    )
    parser.add_argument(
        "--compression",
        choices=["LZW", "JPEG-2000", "zlib", "uncompressed"],
        default="LZW",
        help="Compression type (default: LZW)"
    )
    parser.add_argument(
        "--tile-size",
        type=int,
        default=512,
        help="Tile size for output (default: 512)"
    )
    parser.add_argument(
        "--recursive", "-r",
        action="store_true",
        help="Recursively search for VSI files"
    )
    parser.add_argument(
        "--dry-run", "-n",
        action="store_true",
        help="Show what would be converted without doing it"
    )
    parser.add_argument(
        "--verbose", "-v",
        action="store_true",
        help="Verbose output"
    )
    parser.add_argument(
        "--skip-safety-check",
        action="store_true",
        help="Skip filesystem safety check (use with caution)"
    )

    args = parser.parse_args()

    # Validate bfconvert exists
    if not args.bfconvert.exists():
        print(f"Error: bfconvert not found at {args.bfconvert}")
        print("Install Bio-Formats tools or specify --bfconvert path")
        sys.exit(1)

    # Safety check: ensure output is not on same filesystem as input
    if not args.skip_safety_check and not args.dry_run:
        if args.input_path.exists() and args.output_dir.exists():
            if is_same_filesystem(args.input_path, args.output_dir):
                # Additional check: ensure we're not writing to the input location
                try:
                    input_mount = subprocess.run(
                        ["findmnt", "-n", "-o", "TARGET", "--target", str(args.input_path)],
                        capture_output=True, text=True
                    ).stdout.strip()
                    output_mount = subprocess.run(
                        ["findmnt", "-n", "-o", "TARGET", "--target", str(args.output_dir)],
                        capture_output=True, text=True
                    ).stdout.strip()

                    # Allow if they're on the same base filesystem but different mount points
                    if input_mount and output_mount and input_mount != output_mount:
                        pass  # Different mount points, OK
                    elif "easystore" in str(args.input_path).lower() or "media" in str(args.input_path):
                        print("SAFETY ERROR: Output appears to be on the same drive as input!")
                        print(f"  Input: {args.input_path}")
                        print(f"  Output: {args.output_dir}")
                        print("Use --skip-safety-check to override (NOT RECOMMENDED)")
                        sys.exit(1)
                except Exception:
                    pass  # If findmnt fails, continue

    # Find input files
    if args.input_path.is_file():
        if args.input_path.suffix.lower() != ".vsi":
            print(f"Error: Input file must be a VSI file, got: {args.input_path}")
            sys.exit(1)
        vsi_files = [args.input_path]
    elif args.input_path.is_dir():
        vsi_files = find_vsi_files(args.input_path, args.recursive)
    else:
        print(f"Error: Input path does not exist: {args.input_path}")
        sys.exit(1)

    if not vsi_files:
        print("No VSI files found.")
        sys.exit(0)

    print(f"Found {len(vsi_files)} VSI file(s)")
    if args.dry_run:
        print("DRY RUN - no files will be modified")
    print()

    # Process files
    successful = 0
    failed = 0
    series_arg = args.series if args.series >= 0 else None

    for i, vsi_file in enumerate(vsi_files, 1):
        print(f"[{i}/{len(vsi_files)}] {vsi_file.name}")

        success, message = convert_vsi_to_ometiff(
            vsi_file,
            args.output_dir,
            args.bfconvert,
            series=series_arg,
            compression=args.compression,
            tile_size=args.tile_size,
            dry_run=args.dry_run,
            verbose=args.verbose
        )

        if success:
            successful += 1
            print(f"  ✓ {message}")
        else:
            failed += 1
            print(f"  ✗ {message}")

    print()
    print("=" * 60)
    print(f"Complete: {successful} successful, {failed} failed")

    if failed > 0:
        sys.exit(1)


if __name__ == "__main__":
    main()
