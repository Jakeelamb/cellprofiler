#!/usr/bin/env python3
"""
Verify that OME-TIFF files contain proper metadata for CellProfiler analysis.

This script checks for:
- Physical pixel sizes (required for accurate measurements)
- Objective information
- Channel names
- Image dimensions

Usage:
    python verify_metadata.py <file.ome.tiff>
    python verify_metadata.py <directory>  # Check all OME-TIFF files
"""

import argparse
import subprocess
import sys
import xml.etree.ElementTree as ET
from pathlib import Path

# Bio-Formats showinf path
SHOWINF = Path.home() / "bin" / "bftools" / "bftools" / "showinf"

# OME namespace
OME_NS = {"ome": "http://www.openmicroscopy.org/Schemas/OME/2016-06"}


def get_ome_xml(file_path: Path) -> str | None:
    """Extract OME-XML from file using Bio-Formats."""
    try:
        result = subprocess.run(
            [str(SHOWINF), "-nopix", "-omexml-only", str(file_path)],
            capture_output=True,
            text=True,
            timeout=60
        )
        if result.returncode == 0:
            return result.stdout
        return None
    except Exception as e:
        print(f"Error extracting OME-XML: {e}")
        return None


def try_tifffile(file_path: Path) -> dict | None:
    """Try to read metadata using tifffile."""
    try:
        # Use uv run to access tifffile in the virtual environment
        result = subprocess.run(
            ["uv", "run", "python", "-c", f"""
import tifffile
import json
with tifffile.TiffFile('{file_path}') as tif:
    page = tif.pages[0]
    meta = {{}}
    # Get resolution from TIFF tags
    if hasattr(page, 'tags'):
        tags = page.tags
        if 'XResolution' in tags:
            xres = tags['XResolution'].value
            if isinstance(xres, tuple):
                meta['x_resolution'] = xres[0] / xres[1] if xres[1] else 0
        if 'YResolution' in tags:
            yres = tags['YResolution'].value
            if isinstance(yres, tuple):
                meta['y_resolution'] = yres[0] / yres[1] if yres[1] else 0
        if 'ResolutionUnit' in tags:
            unit = tags['ResolutionUnit'].value
            meta['resolution_unit'] = str(unit)

    # Get OME-XML if present
    if tif.ome_metadata:
        meta['has_ome'] = True
        meta['ome_preview'] = tif.ome_metadata[:500]
    else:
        meta['has_ome'] = False

    meta['shape'] = list(page.shape)
    meta['dtype'] = str(page.dtype)
    print(json.dumps(meta))
"""],
            capture_output=True,
            text=True,
            timeout=30
        )
        if result.returncode == 0:
            import json
            return json.loads(result.stdout.strip())
    except Exception as e:
        print(f"tifffile error: {e}")
    return None


def parse_ome_metadata(ome_xml: str) -> dict:
    """Parse OME-XML and extract key metadata."""
    metadata = {
        "physical_size_x": None,
        "physical_size_y": None,
        "physical_size_z": None,
        "unit_x": None,
        "unit_y": None,
        "objective": None,
        "magnification": None,
        "channels": [],
        "dimensions": {},
        "valid_for_cellprofiler": False,
        "warnings": [],
    }

    try:
        # Find the XML start (skip any preamble)
        xml_start = ome_xml.find("<?xml")
        if xml_start < 0:
            xml_start = ome_xml.find("<OME")
        if xml_start < 0:
            metadata["warnings"].append("No valid OME-XML found")
            return metadata

        ome_xml = ome_xml[xml_start:]

        # Parse XML
        root = ET.fromstring(ome_xml)

        # Handle namespace
        ns = {"ome": "http://www.openmicroscopy.org/Schemas/OME/2016-06"}

        # Find Pixels element (contains physical sizes)
        pixels = root.find(".//ome:Pixels", ns)
        if pixels is None:
            # Try without namespace
            pixels = root.find(".//Pixels")

        if pixels is not None:
            metadata["physical_size_x"] = pixels.get("PhysicalSizeX")
            metadata["physical_size_y"] = pixels.get("PhysicalSizeY")
            metadata["physical_size_z"] = pixels.get("PhysicalSizeZ")
            metadata["unit_x"] = pixels.get("PhysicalSizeXUnit", "µm")
            metadata["unit_y"] = pixels.get("PhysicalSizeYUnit", "µm")
            metadata["dimensions"] = {
                "width": pixels.get("SizeX"),
                "height": pixels.get("SizeY"),
                "channels": pixels.get("SizeC"),
                "z_slices": pixels.get("SizeZ"),
                "timepoints": pixels.get("SizeT"),
                "type": pixels.get("Type"),
            }

        # Find Objective
        objective = root.find(".//ome:Objective", ns)
        if objective is None:
            objective = root.find(".//Objective")

        if objective is not None:
            metadata["objective"] = objective.get("Model")
            metadata["magnification"] = objective.get("NominalMagnification")

        # Find Channels
        for channel in root.findall(".//ome:Channel", ns) or root.findall(".//Channel"):
            ch_name = channel.get("Name", "Unknown")
            metadata["channels"].append(ch_name)

        # Validate for CellProfiler
        if metadata["physical_size_x"] and metadata["physical_size_y"]:
            try:
                px = float(metadata["physical_size_x"])
                py = float(metadata["physical_size_y"])
                if px > 0 and py > 0:
                    metadata["valid_for_cellprofiler"] = True
            except ValueError:
                metadata["warnings"].append("Could not parse pixel sizes as numbers")

        if not metadata["valid_for_cellprofiler"]:
            metadata["warnings"].append("Missing or invalid physical pixel sizes")

    except ET.ParseError as e:
        metadata["warnings"].append(f"XML parse error: {e}")

    return metadata


def check_file(file_path: Path) -> dict:
    """Check a single file for metadata."""
    result = {
        "file": str(file_path),
        "exists": file_path.exists(),
        "size_mb": 0,
        "metadata": None,
        "status": "unknown"
    }

    if not file_path.exists():
        result["status"] = "file_not_found"
        return result

    result["size_mb"] = file_path.stat().st_size / (1024 * 1024)

    # Try Bio-Formats first
    ome_xml = get_ome_xml(file_path)
    if ome_xml:
        result["metadata"] = parse_ome_metadata(ome_xml)
        if result["metadata"]["valid_for_cellprofiler"]:
            result["status"] = "valid"
        else:
            result["status"] = "missing_metadata"
    else:
        # Fallback to tifffile
        tiff_meta = try_tifffile(file_path)
        if tiff_meta:
            result["metadata"] = {
                "from_tifffile": True,
                "has_ome": tiff_meta.get("has_ome", False),
                "shape": tiff_meta.get("shape"),
                "dtype": tiff_meta.get("dtype"),
                "x_resolution": tiff_meta.get("x_resolution"),
                "y_resolution": tiff_meta.get("y_resolution"),
            }
            result["status"] = "tiff_only"
        else:
            result["status"] = "cannot_read"

    return result


def print_result(result: dict, verbose: bool = False):
    """Print check result."""
    status_icons = {
        "valid": "✓",
        "missing_metadata": "⚠",
        "file_not_found": "✗",
        "cannot_read": "✗",
        "tiff_only": "○",
        "unknown": "?"
    }

    icon = status_icons.get(result["status"], "?")
    print(f"{icon} {result['file']} ({result['size_mb']:.1f} MB)")

    meta = result.get("metadata")
    if meta:
        if meta.get("from_tifffile"):
            print(f"    Format: Standard TIFF (no OME metadata)")
            if meta.get("x_resolution"):
                print(f"    Resolution: {meta.get('x_resolution')} x {meta.get('y_resolution')}")
        else:
            if meta.get("physical_size_x"):
                print(f"    Pixel size: {meta['physical_size_x']} x {meta['physical_size_y']} {meta.get('unit_x', 'µm')}")
            if meta.get("objective"):
                print(f"    Objective: {meta['objective']} ({meta.get('magnification')}x)")
            if meta.get("channels"):
                print(f"    Channels: {', '.join(meta['channels'])}")
            if meta.get("dimensions"):
                dims = meta["dimensions"]
                print(f"    Size: {dims.get('width')} x {dims.get('height')} ({dims.get('type')})")

            if meta.get("warnings"):
                for warn in meta["warnings"]:
                    print(f"    ⚠ {warn}")

    if result["status"] == "valid":
        print("    → Ready for CellProfiler analysis")
    elif result["status"] == "missing_metadata":
        print("    → Missing pixel size - measurements will be in pixels only")


def main():
    parser = argparse.ArgumentParser(description="Verify OME-TIFF metadata for CellProfiler")
    parser.add_argument("path", type=Path, help="File or directory to check")
    parser.add_argument("--verbose", "-v", action="store_true")

    args = parser.parse_args()

    if args.path.is_file():
        files = [args.path]
    elif args.path.is_dir():
        files = list(args.path.glob("*.ome.tiff")) + list(args.path.glob("*.ome.tif"))
        files.extend(args.path.glob("*.tiff"))
        files.extend(args.path.glob("*.tif"))
        files = sorted(set(files))
    else:
        print(f"Error: {args.path} not found")
        sys.exit(1)

    if not files:
        print("No TIFF files found")
        sys.exit(0)

    valid = 0
    invalid = 0

    for f in files:
        result = check_file(f)
        print_result(result, args.verbose)
        if result["status"] == "valid":
            valid += 1
        else:
            invalid += 1

    print()
    print(f"Summary: {valid} valid, {invalid} with issues")


if __name__ == "__main__":
    main()
