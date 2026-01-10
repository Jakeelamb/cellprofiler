#!/usr/bin/env python3
"""
Extract a single channel from RGB OME-TIFF files while preserving metadata.

This is useful for CellProfiler analysis when you only need one channel
(typically green for brightfield microscopy with fluorescence).

Usage:
    # Extract green channel from a single file
    uv run python extract_channel.py image.ome.tiff output/ --channel 1

    # Extract green channel from all OME-TIFF files in a directory
    uv run python extract_channel.py input_dir/ output_dir/ --green

    # Process multiple files
    uv run python extract_channel.py input/ output/ --green --recursive

The script preserves all OME-XML metadata including pixel sizes.
"""

import argparse
import sys
import xml.etree.ElementTree as ET
from pathlib import Path

import numpy as np
import tifffile


def extract_channel_from_tiff(
    input_path: Path,
    output_dir: Path,
    channel: int = 1,  # 0=R, 1=G, 2=B
    compression: str = "deflate",
    tile_size: int = 512,
    dry_run: bool = False,
    verbose: bool = False,
) -> tuple[bool, str]:
    """
    Extract a single channel from an RGB TIFF file.

    Preserves OME-XML metadata (pixel sizes, objective info, etc.)
    """
    channel_names = {0: "red", 1: "green", 2: "blue"}
    channel_name = channel_names.get(channel, f"ch{channel}")

    output_name = f"{input_path.stem}_{channel_name}.ome.tiff"
    output_path = output_dir / output_name

    if dry_run:
        return True, f"Would extract channel {channel} ({channel_name}): {input_path} -> {output_path}"

    try:
        with tifffile.TiffFile(input_path) as tif:
            # Check if it's RGB
            page = tif.pages[0]
            shape = page.shape

            if len(shape) < 3 or shape[-1] not in [3, 4]:
                return False, f"Not an RGB image (shape: {shape})"

            if verbose:
                print(f"  Input shape: {shape}, dtype: {page.dtype}")

            # Get OME metadata if present
            ome_xml = None
            if tif.ome_metadata:
                ome_xml = tif.ome_metadata
                if verbose:
                    print(f"  Found OME metadata")

            # Read the image data
            # For very large files, we'd want to do this tile by tile,
            # but for now let's read the whole thing
            data = page.asarray()

            # Extract the channel
            if len(data.shape) == 3:  # HxWxC
                channel_data = data[:, :, channel]
            elif len(data.shape) == 4:  # ZxHxWxC or TxHxWxC
                channel_data = data[:, :, :, channel]
            else:
                return False, f"Unexpected shape: {data.shape}"

            if verbose:
                print(f"  Extracted channel shape: {channel_data.shape}")

            # Update OME-XML metadata to reflect single channel
            metadata = {}
            if ome_xml:
                try:
                    metadata = update_ome_for_single_channel(ome_xml, channel, channel_name)
                except Exception as e:
                    if verbose:
                        print(f"  Warning: Could not update OME-XML: {e}")
                    # Continue without updated metadata

            # Create output directory
            output_dir.mkdir(parents=True, exist_ok=True)

            # Write the output
            write_kwargs = {
                "bigtiff": True,
                "tile": (tile_size, tile_size),
                "compression": compression,
                "photometric": "minisblack",  # Grayscale
            }

            # Add resolution if available
            if "resolution" in metadata:
                write_kwargs["resolution"] = metadata["resolution"]
                write_kwargs["resolutionunit"] = metadata["resolutionunit"]

            # Add OME-XML description if available
            if "description" in metadata:
                write_kwargs["description"] = metadata["description"]

            tifffile.imwrite(output_path, channel_data, **write_kwargs)

            size_mb = output_path.stat().st_size / (1024 * 1024)
            return True, f"Created {output_path} ({size_mb:.1f} MB)"

    except Exception as e:
        return False, f"Error: {e}"


def update_ome_for_single_channel(ome_xml: str, channel_idx: int, channel_name: str) -> dict:
    """
    Parse OME-XML and extract relevant metadata for single-channel output.
    Returns dict with resolution info for tifffile.
    """
    metadata = {}

    try:
        # Find XML start
        xml_start = ome_xml.find("<?xml")
        if xml_start < 0:
            xml_start = ome_xml.find("<OME")
        if xml_start >= 0:
            ome_xml = ome_xml[xml_start:]

        root = ET.fromstring(ome_xml)
        ns = {"ome": "http://www.openmicroscopy.org/Schemas/OME/2016-06"}

        # Find Pixels element
        pixels = root.find(".//ome:Pixels", ns)
        if pixels is None:
            pixels = root.find(".//Pixels")

        if pixels is not None:
            # Get physical sizes
            phys_x = pixels.get("PhysicalSizeX")
            phys_y = pixels.get("PhysicalSizeY")
            unit_x = pixels.get("PhysicalSizeXUnit", "µm")

            if phys_x and phys_y:
                # Convert to resolution (pixels per unit)
                try:
                    px = float(phys_x)
                    py = float(phys_y)

                    if unit_x in ["µm", "um", "micron"]:
                        # TIFF resolution is in pixels per cm
                        # 1 µm = 0.0001 cm
                        res_x = 1.0 / (px * 0.0001)  # pixels per cm
                        res_y = 1.0 / (py * 0.0001)
                        metadata["resolution"] = (res_x, res_y)
                        metadata["resolutionunit"] = tifffile.RESUNIT.CENTIMETER

                except (ValueError, ZeroDivisionError):
                    pass

            # Create minimal OME-XML for single channel
            # This preserves pixel size info that CellProfiler can read
            # Note: Use "um" instead of "µm" for ASCII compatibility
            unit_ascii = "um" if unit_x in ["µm", "um", "micron"] else unit_x
            desc = f"""<?xml version="1.0" encoding="UTF-8"?>
<OME xmlns="http://www.openmicroscopy.org/Schemas/OME/2016-06">
  <Image ID="Image:0" Name="{channel_name}">
    <Pixels DimensionOrder="XYCZT" ID="Pixels:0"
            PhysicalSizeX="{phys_x}" PhysicalSizeXUnit="{unit_ascii}"
            PhysicalSizeY="{phys_y}" PhysicalSizeYUnit="{unit_ascii}"
            SizeC="1" SizeT="1" SizeZ="1"
            Type="{pixels.get('Type', 'uint8')}">
      <Channel ID="Channel:0:0" Name="{channel_name.capitalize()}" SamplesPerPixel="1"/>
    </Pixels>
  </Image>
</OME>"""
            metadata["description"] = desc

    except ET.ParseError:
        pass

    return metadata


def find_tiff_files(root_dir: Path, recursive: bool = True) -> list[Path]:
    """Find all TIFF files in directory."""
    patterns = ["*.ome.tiff", "*.ome.tif", "*.tiff", "*.tif"]
    files = []

    for pattern in patterns:
        if recursive:
            files.extend(root_dir.glob(f"**/{pattern}"))
        else:
            files.extend(root_dir.glob(pattern))

    return sorted(set(files))


def main():
    parser = argparse.ArgumentParser(
        description="Extract a single channel from RGB TIFF files",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )
    parser.add_argument(
        "input_path",
        type=Path,
        help="Input TIFF file or directory"
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
        "--green", "-g",
        action="store_true",
        help="Extract green channel (same as --channel 1)"
    )
    parser.add_argument(
        "--red",
        action="store_true",
        help="Extract red channel (same as --channel 0)"
    )
    parser.add_argument(
        "--blue",
        action="store_true",
        help="Extract blue channel (same as --channel 2)"
    )
    parser.add_argument(
        "--compression",
        choices=["deflate", "lzw", "zlib", "none"],
        default="deflate",
        help="Compression type (default: deflate)"
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
        help="Recursively search for TIFF files"
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

    args = parser.parse_args()

    # Handle channel shortcuts
    channel = args.channel
    if args.green:
        channel = 1
    elif args.red:
        channel = 0
    elif args.blue:
        channel = 2

    channel_names = {0: "Red", 1: "Green", 2: "Blue"}
    print(f"Extracting {channel_names.get(channel, f'channel {channel}')} channel")

    # Find input files
    if args.input_path.is_file():
        files = [args.input_path]
    elif args.input_path.is_dir():
        files = find_tiff_files(args.input_path, args.recursive)
    else:
        print(f"Error: {args.input_path} not found")
        sys.exit(1)

    if not files:
        print("No TIFF files found")
        sys.exit(0)

    print(f"Found {len(files)} file(s)")
    if args.dry_run:
        print("DRY RUN - no files will be created")
    print()

    successful = 0
    failed = 0

    compression = args.compression if args.compression != "none" else None

    for i, f in enumerate(files, 1):
        print(f"[{i}/{len(files)}] {f.name}")

        success, msg = extract_channel_from_tiff(
            f,
            args.output_dir,
            channel=channel,
            compression=compression,
            tile_size=args.tile_size,
            dry_run=args.dry_run,
            verbose=args.verbose,
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
