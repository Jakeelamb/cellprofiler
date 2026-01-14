#!/usr/bin/env python3
"""
Extract green channel from large BigTIFF files with compression.

This script processes large OME-BigTIFF files by:
1. Reading the image (can handle multi-GB files with sufficient RAM)
2. Extracting the green channel
3. Writing with LZW/deflate compression
4. Preserving OME-XML metadata

Typical size reduction: 8GB uncompressed RGB -> 500MB-1GB compressed green

Usage:
    # Single file
    uv run python btf_to_green.py /path/to/file.ome.btf /output/dir

    # Directory (all .btf files)
    uv run python btf_to_green.py /path/to/input_dir /output/dir --recursive

    # Dry run
    uv run python btf_to_green.py /path/to/input_dir /output/dir --dry-run
"""

import argparse
import gc
import os
import sys
import xml.etree.ElementTree as ET
from pathlib import Path

import numpy as np
import tifffile
import zarr
from tqdm import tqdm


def extract_pixel_metadata(ome_xml: str) -> dict:
    """Extract pixel size and other metadata from OME-XML."""
    metadata = {}

    try:
        xml_start = ome_xml.find("<?xml")
        if xml_start < 0:
            xml_start = ome_xml.find("<OME")
        if xml_start >= 0:
            ome_xml = ome_xml[xml_start:]

        root = ET.fromstring(ome_xml)

        # Try different namespace patterns
        for ns_prefix in [
            {"ome": "http://www.openmicroscopy.org/Schemas/OME/2015-01"},
            {"ome": "http://www.openmicroscopy.org/Schemas/OME/2016-06"},
            {},
        ]:
            pixels = root.find(".//ome:Pixels", ns_prefix) if ns_prefix else root.find(".//Pixels")
            if pixels is not None:
                break

        if pixels is not None:
            phys_x = pixels.get("PhysicalSizeX")
            phys_y = pixels.get("PhysicalSizeY")
            phys_x_unit = pixels.get("PhysicalSizeXUnit", "µm")
            phys_y_unit = pixels.get("PhysicalSizeYUnit", "µm")

            if phys_x:
                metadata["PhysicalSizeX"] = float(phys_x)
                metadata["PhysicalSizeXUnit"] = phys_x_unit
            if phys_y:
                metadata["PhysicalSizeY"] = float(phys_y)
                metadata["PhysicalSizeYUnit"] = phys_y_unit

            # Get image type
            metadata["Type"] = pixels.get("Type", "uint8")
            metadata["SizeX"] = pixels.get("SizeX")
            metadata["SizeY"] = pixels.get("SizeY")

        # Try to get objective info
        objective = root.find(".//*[@NominalMagnification]")
        if objective is not None:
            metadata["Magnification"] = objective.get("NominalMagnification")
            metadata["ObjectiveModel"] = objective.get("Model")

    except Exception as e:
        print(f"  Warning: Could not parse OME-XML: {e}")

    return metadata


def create_ome_xml(metadata: dict, width: int, height: int) -> str:
    """Create minimal OME-XML for single-channel grayscale output."""
    phys_x = metadata.get("PhysicalSizeX", "")
    phys_y = metadata.get("PhysicalSizeY", "")
    unit = metadata.get("PhysicalSizeXUnit", "um")
    # Use ASCII-safe unit
    if unit in ["µm", "um", "micron"]:
        unit = "um"

    dtype = metadata.get("Type", "uint8")

    ome_xml = f'''<?xml version="1.0" encoding="UTF-8"?>
<OME xmlns="http://www.openmicroscopy.org/Schemas/OME/2016-06">
  <Image ID="Image:0" Name="green_channel">
    <Pixels DimensionOrder="XYC" ID="Pixels:0"
            PhysicalSizeX="{phys_x}" PhysicalSizeXUnit="{unit}"
            PhysicalSizeY="{phys_y}" PhysicalSizeYUnit="{unit}"
            SizeX="{width}" SizeY="{height}" SizeC="1" SizeZ="1" SizeT="1"
            Type="{dtype}">
      <Channel ID="Channel:0:0" Name="Green" SamplesPerPixel="1"/>
    </Pixels>
  </Image>
</OME>'''
    return ome_xml


def process_btf_file_chunked(
    input_path: Path,
    output_dir: Path,
    compression: str = "deflate",
    compression_level: int = 6,
    tile_size: int = 512,
    channel: int = 1,
    dry_run: bool = False,
    verbose: bool = False,
    chunk_rows: int = 2048,
) -> tuple[bool, str]:
    """
    Extract green channel from a BigTIFF file with compression.

    Uses memory-mapped temporary file to process the image in chunks
    without loading the entire file into RAM at once.

    Returns (success, message) tuple.
    """
    channel_names = {0: "red", 1: "green", 2: "blue"}
    channel_name = channel_names.get(channel, f"ch{channel}")

    output_name = f"{input_path.stem.replace('.ome', '')}_{channel_name}.ome.tiff"
    output_path = output_dir / output_name

    if dry_run:
        return True, f"Would process: {input_path.name} -> {output_name}"

    temp_path = None
    try:
        with tifffile.TiffFile(input_path) as tif:
            page = tif.pages[0]
            shape = page.shape

            if len(shape) < 3 or shape[-1] not in [3, 4]:
                return False, f"Not an RGB image (shape: {shape})"

            height, width = shape[0], shape[1]
            dtype = page.dtype
            input_size_gb = input_path.stat().st_size / (1024**3)

            if verbose:
                print(f"  Input: {width}x{height}, {shape[-1]} channels, {input_size_gb:.1f} GB")

            # Extract metadata
            metadata = {}
            if tif.ome_metadata:
                metadata = extract_pixel_metadata(tif.ome_metadata)
                if verbose and metadata.get("PhysicalSizeX"):
                    print(f"  Pixel size: {metadata['PhysicalSizeX']} {metadata.get('PhysicalSizeXUnit', 'µm')}")

            # Create output directory
            output_dir.mkdir(parents=True, exist_ok=True)

            # Create OME-XML
            ome_xml = create_ome_xml(metadata, width, height)

            # Align chunk_rows to tile_size for efficient writing
            chunk_rows = ((chunk_rows + tile_size - 1) // tile_size) * tile_size

            # Calculate number of chunks
            num_chunks = (height + chunk_rows - 1) // chunk_rows

            if verbose:
                chunk_mem_mb = (chunk_rows * width * np.dtype(dtype).itemsize) / (1024**2)
                print(f"  Processing in {num_chunks} chunks (~{chunk_mem_mb:.0f} MB each)")

            # Use zarr store for memory-efficient access (reads tiles on demand)
            store = tif.aszarr()
            z = zarr.open(store, mode='r')

            # Handle zarr array structure - may be nested for OME-TIFF
            if hasattr(z, 'shape'):
                zarray = z
            elif '0' in z:
                zarray = z['0']
            else:
                keys = list(z.keys())
                if keys:
                    zarray = z[keys[0]]
                else:
                    return False, "Could not access image data via zarr"

            # Create memory-mapped temp file for the single-channel output
            # Use output_dir for temp files to avoid filling /tmp (which may be RAM-based)
            temp_path = output_dir / f".tmp_{input_path.stem}_{os.getpid()}.dat"
            temp_file = None  # Will be created by memmap

            # Create memory-mapped array
            memmap_shape = (height, width)
            memmap_array = np.memmap(temp_path, dtype=dtype, mode='w+', shape=memmap_shape)

            if verbose:
                print(f"  Extracting {channel_name} channel...")

            # Fill the memmap chunk by chunk
            for chunk_idx in tqdm(range(num_chunks), desc="  Reading", disable=not verbose):
                row_start = chunk_idx * chunk_rows
                row_end = min(row_start + chunk_rows, height)

                # Read chunk via zarr (only loads necessary tiles from disk)
                memmap_array[row_start:row_end, :] = zarray[row_start:row_end, :, channel]

            # Flush to disk
            memmap_array.flush()

            store.close()

            if verbose:
                print(f"  Writing compressed output...")

            # Write the memmap to compressed TIFF
            tifffile.imwrite(
                output_path,
                memmap_array,
                bigtiff=True,
                tile=(tile_size, tile_size),
                compression=compression,
                compressionargs={"level": compression_level},
                photometric="minisblack",
                description=ome_xml,
            )

            # Clean up memmap
            del memmap_array
            gc.collect()

            output_size_mb = output_path.stat().st_size / (1024**2)
            reduction = input_size_gb * 1024 / output_size_mb

            return True, f"Created {output_path.name} ({output_size_mb:.0f} MB, {reduction:.1f}x reduction)"

    except MemoryError:
        return False, "Out of memory - try reducing --chunk-rows"
    except Exception as e:
        import traceback
        return False, f"Error: {e}\n{traceback.format_exc()}"
    finally:
        # Clean up temp file
        if temp_path is not None:
            try:
                os.unlink(temp_path)
            except Exception:
                pass


def process_btf_file(
    input_path: Path,
    output_dir: Path,
    compression: str = "deflate",
    compression_level: int = 6,
    tile_size: int = 512,
    channel: int = 1,
    dry_run: bool = False,
    verbose: bool = False,
) -> tuple[bool, str]:
    """
    Extract green channel from a BigTIFF file with compression.
    Legacy function - loads entire image into memory.

    Returns (success, message) tuple.
    """
    channel_names = {0: "red", 1: "green", 2: "blue"}
    channel_name = channel_names.get(channel, f"ch{channel}")

    output_name = f"{input_path.stem.replace('.ome', '')}_{channel_name}.ome.tiff"
    output_path = output_dir / output_name

    if dry_run:
        return True, f"Would process: {input_path.name} -> {output_name}"

    try:
        with tifffile.TiffFile(input_path) as tif:
            page = tif.pages[0]
            shape = page.shape

            if len(shape) < 3 or shape[-1] not in [3, 4]:
                return False, f"Not an RGB image (shape: {shape})"

            height, width = shape[0], shape[1]
            input_size_gb = input_path.stat().st_size / (1024**3)

            if verbose:
                print(f"  Input: {width}x{height}, {shape[-1]} channels, {input_size_gb:.1f} GB")

            # Extract metadata
            metadata = {}
            if tif.ome_metadata:
                metadata = extract_pixel_metadata(tif.ome_metadata)
                if verbose and metadata.get("PhysicalSizeX"):
                    print(f"  Pixel size: {metadata['PhysicalSizeX']} {metadata.get('PhysicalSizeXUnit', 'µm')}")

            # Read full image - this requires sufficient RAM
            if verbose:
                print(f"  Reading image data...")

            data = page.asarray()

            # Extract channel
            if verbose:
                print(f"  Extracting {channel_name} channel...")

            green_data = data[:, :, channel]

            # Free RGB data
            del data
            gc.collect()

            # Create output directory
            output_dir.mkdir(parents=True, exist_ok=True)

            # Create OME-XML
            ome_xml = create_ome_xml(metadata, width, height)

            # Write with compression
            if verbose:
                print(f"  Writing compressed output...")

            tifffile.imwrite(
                output_path,
                green_data,
                bigtiff=True,
                tile=(tile_size, tile_size),
                compression=compression,
                compressionargs={"level": compression_level},
                photometric="minisblack",
                description=ome_xml,
            )

            del green_data
            gc.collect()

            output_size_mb = output_path.stat().st_size / (1024**2)
            reduction = input_size_gb * 1024 / output_size_mb

            return True, f"Created {output_path.name} ({output_size_mb:.0f} MB, {reduction:.1f}x reduction)"

    except MemoryError:
        return False, "Out of memory - file too large for available RAM"
    except Exception as e:
        return False, f"Error: {e}"


def find_btf_files(root_dir: Path, recursive: bool = True) -> list[Path]:
    """Find all BigTIFF files in directory."""
    patterns = ["*.btf", "*.ome.btf", "*.BTF"]
    files = []

    for pattern in patterns:
        if recursive:
            files.extend(root_dir.glob(f"**/{pattern}"))
        else:
            files.extend(root_dir.glob(pattern))

    return sorted(set(files))


def main():
    parser = argparse.ArgumentParser(
        description="Extract green channel from BigTIFF files with compression",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )
    parser.add_argument(
        "input_path",
        type=Path,
        help="Input BigTIFF file or directory"
    )
    parser.add_argument(
        "output_dir",
        type=Path,
        help="Output directory"
    )
    parser.add_argument(
        "--channel",
        type=int,
        default=1,
        help="Channel to extract (0=Red, 1=Green, 2=Blue). Default: 1 (green)"
    )
    parser.add_argument(
        "--compression",
        choices=["deflate", "lzw", "zlib", "zstd"],
        default="deflate",
        help="Compression type (default: deflate)"
    )
    parser.add_argument(
        "--compression-level",
        type=int,
        default=6,
        help="Compression level 1-9 (default: 6)"
    )
    parser.add_argument(
        "--tile-size",
        type=int,
        default=512,
        help="Tile size (default: 512)"
    )
    parser.add_argument(
        "--recursive", "-r",
        action="store_true",
        help="Recursively search for BTF files"
    )
    parser.add_argument(
        "--dry-run", "-n",
        action="store_true",
        help="Show what would be done"
    )
    parser.add_argument(
        "--verbose", "-v",
        action="store_true",
        help="Verbose output"
    )
    parser.add_argument(
        "--chunk-rows",
        type=int,
        default=2048,
        help="Number of rows to process at once (default: 2048). Lower = less RAM"
    )
    parser.add_argument(
        "--no-chunked",
        action="store_true",
        help="Disable chunked processing (loads entire image into RAM)"
    )

    args = parser.parse_args()

    # Find input files
    if args.input_path.is_file():
        files = [args.input_path]
    elif args.input_path.is_dir():
        files = find_btf_files(args.input_path, args.recursive)
    else:
        print(f"Error: {args.input_path} not found")
        sys.exit(1)

    if not files:
        print("No BigTIFF files found")
        sys.exit(0)

    channel_names = {0: "Red", 1: "Green", 2: "Blue"}
    print(f"Extracting {channel_names.get(args.channel, f'channel {args.channel}')} channel")
    print(f"Found {len(files)} file(s)")
    if args.dry_run:
        print("DRY RUN - no files will be created")
    print()

    successful = 0
    failed = 0

    for i, f in enumerate(files, 1):
        print(f"[{i}/{len(files)}] {f.name}")

        if args.no_chunked:
            success, msg = process_btf_file(
                f,
                args.output_dir,
                compression=args.compression,
                compression_level=args.compression_level,
                tile_size=args.tile_size,
                channel=args.channel,
                dry_run=args.dry_run,
                verbose=args.verbose,
            )
        else:
            success, msg = process_btf_file_chunked(
                f,
                args.output_dir,
                compression=args.compression,
                compression_level=args.compression_level,
                tile_size=args.tile_size,
                channel=args.channel,
                dry_run=args.dry_run,
                verbose=args.verbose,
                chunk_rows=args.chunk_rows,
            )

        if success:
            successful += 1
            print(f"  ✓ {msg}")
        else:
            failed += 1
            print(f"  ✗ {msg}")

    print()
    print("=" * 60)
    print(f"Complete: {successful} successful, {failed} failed")

    if failed > 0:
        sys.exit(1)


if __name__ == "__main__":
    main()
