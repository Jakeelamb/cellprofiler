# Microscopy Image Processing for CellProfiler

Tools for converting Olympus VSI and BigTIFF microscopy files to CellProfiler-compatible OME-TIFF format while **preserving all metadata** (pixel sizes, objective info, stage positions).

## Features

- **Metadata preservation**: Keeps physical pixel sizes (µm/pixel), objective info, channel names
- **Safe processing**: Reads from source drive (mounted read-only), writes to separate destination
- **Large file support**: Handles multi-gigabyte BigTIFF and VSI pyramid images
- **CellProfiler ready**: Output files include all metadata needed for accurate measurements

## Requirements

- Python 3.12+
- Bio-Formats command line tools (bftools)
- tifffile, numpy, imagecodecs (via uv)

### Installing Bio-Formats Tools

```bash
mkdir -p ~/bin/bftools
cd ~/bin/bftools
curl -L -o bftools.zip "https://downloads.openmicroscopy.org/bio-formats/8.0.1/artifacts/bftools.zip"
unzip bftools.zip
chmod +x bftools/bfconvert bftools/showinf
```

## Usage

### 1. Mount Source Drive (Read-Only for Safety)

```bash
# Mount external drive as read-only to prevent accidental writes
udisksctl mount -b /dev/sda1 --options ro
```

### 2. Convert VSI to OME-TIFF

```bash
# Single file
python vsi_to_ometiff.py /path/to/image.vsi /output/dir

# Entire directory (recursive)
python vsi_to_ometiff.py /path/to/input_dir /output/dir --recursive

# Preview what would be converted
python vsi_to_ometiff.py /path/to/input_dir /output/dir --dry-run

# Export specific resolution level (0 = full resolution, default)
python vsi_to_ometiff.py /path/to/image.vsi /output/dir --series 0
```

### 3. Extract Green Channel Only (Recommended for Analysis)

For CellProfiler analysis, you typically only need the green channel. This reduces file size by ~66%:

```bash
# Extract green channel from converted OME-TIFF
uv run python extract_channel.py /path/to/image.ome.tiff /output/dir --green

# Process entire directory
uv run python extract_channel.py /input/dir /output/dir --green --recursive

# Or specify channel by number (0=Red, 1=Green, 2=Blue)
uv run python extract_channel.py /path/to/image.ome.tiff /output/dir --channel 1
```

### 4. Process BigTIFF Files (Alternative to VSI)

For large OME-BigTIFF files (`.ome.btf`), extract the green channel directly with compression:

```bash
# Single BigTIFF file - achieves ~8x size reduction
uv run python btf_to_green.py /path/to/image.ome.btf /output/dir

# Process entire directory
uv run python btf_to_green.py /input/dir /output/dir --recursive

# Preview what would be processed
uv run python btf_to_green.py /input/dir /output/dir --dry-run

# Adjust compression (default: deflate level 6)
uv run python btf_to_green.py /path/to/image.ome.btf /output/dir --compression lzw
```

Typical results:
- 8.2 GB uncompressed RGB → ~1 GB compressed green channel
- Metadata preserved (pixel size, channel info)
- Output ready for CellProfiler analysis

### 5. Verify Metadata Preservation

```bash
# Check a single file
python verify_metadata.py /output/dir/image.ome.tiff

# Check all files in directory
python verify_metadata.py /output/dir
```

Expected output:
```
✓ /output/dir/image.ome.tiff (7.1 MB)
    Pixel size: 0.12 x 0.12 µm
    Objective: LUCPLFLN PH (40.0x)
    Channels: BF
    Size: 3088 x 2076 (uint8)
    → Ready for CellProfiler analysis
```

## File Format Notes

### VSI Files (Olympus/Evident)
- `.vsi` file is just an index (~100KB-1MB)
- Actual image data is in companion `_filename_/` folder
- Contains resolution pyramid (multiple series)
- Must keep `.vsi` and `_filename_/` folder together

### BigTIFF Files (.ome.btf)
- Large uncompressed RGB images (8+ GB each)
- Contains full OME-XML metadata
- Use `btf_to_green.py` to extract green channel with compression
- Achieves ~8x size reduction while preserving metadata

### OME-TIFF Output
- Single file containing image + all metadata
- Uses OME-XML standard for metadata
- Compressed (deflate/LZW) by default
- BigTIFF format for large files

## Directory Structure

```
microscopy_data/              # Output directory (on local drive)
├── image1.ome.tiff          # Converted images with metadata
├── image2.ome.tiff
└── ...

/run/media/jake/easystore/   # Source drive (mounted read-only)
├── image1.vsi
├── _image1_/                # VSI companion data
│   └── stack1/
│       └── frame_t_0.ets    # Actual image data
└── ...
```

## Safety Features

The scripts include safety checks to prevent accidental overwrites:
- Refuses to write to the same filesystem as source
- Source drive should be mounted read-only
- All processing creates new files in the destination directory

## CellProfiler Integration

The output OME-TIFF files are directly compatible with CellProfiler:
1. Open CellProfiler
2. Images module → drag and drop OME-TIFF files
3. Metadata will be automatically extracted
4. Measurements will be in real units (µm) instead of pixels

## Dependencies

```toml
[project]
dependencies = [
    "imagecodecs>=2025.8.2",
    "numpy>=2.3.4",
    "tifffile>=2025.10.16",
    "tqdm>=4.66.0",
    "PyYAML>=6.0",
    "psutil>=5.9.0",
]
```

Install with: `uv sync`
