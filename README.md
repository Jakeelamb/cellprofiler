# CellProfiler Tools

Tools for converting Olympus VSI and BigTIFF microscopy files to CellProfiler-compatible OME-TIFF format while **preserving all metadata** (pixel sizes, objective info, stage positions).

## Features

- **Metadata preservation**: Keeps physical pixel sizes (um/pixel), objective info, channel names
- **Safe processing**: Reads from source drive (mounted read-only), writes to separate destination
- **Large file support**: Handles multi-gigabyte BigTIFF and VSI pyramid images
- **CellProfiler ready**: Output files include all metadata needed for accurate measurements

## Project Structure

```
cellprofiler-tools/
├── src/cellprofiler_tools/      # Python package
│   ├── converters/              # Image format converters
│   │   ├── btf_to_green.py      # BigTIFF green channel extraction
│   │   ├── extract_channel.py   # Generic channel extraction
│   │   └── vsi_to_ometiff.py    # VSI to OME-TIFF conversion
│   ├── batch/                   # Batch processing
│   │   ├── process.py           # Sequential batch processing
│   │   └── process_parallel.py  # Parallel processing with RAM monitoring
│   └── analysis/                # Analysis tools
│       ├── calculate_iod.py     # IOD calculation for genome size
│       └── verify_metadata.py   # Metadata validation
├── pipelines/                   # CellProfiler pipelines
│   ├── oil_immersion.cppipe     # For nucleus-only images
│   └── brightfield.cppipe       # For nucleus + cytoplasm images
├── scripts/                     # Shell scripts
│   ├── process_all.sh           # Parallel batch processor
│   ├── run_batch.sh             # Conservative batch runner
│   └── mount_easystore.sh       # Drive mounting helper
└── docs/                        # Documentation
    ├── PIPELINE_DOCUMENTATION.md
    └── CONVERSION_REPORT.md
```

## Installation

```bash
# Clone the repository
git clone <repo-url>
cd cellprofiler-tools

# Install dependencies
uv sync

# Install in development mode (optional)
uv pip install -e .
```

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
python src/cellprofiler_tools/converters/vsi_to_ometiff.py /path/to/image.vsi /output/dir

# Entire directory (recursive)
python src/cellprofiler_tools/converters/vsi_to_ometiff.py /path/to/input_dir /output/dir --recursive

# Preview what would be converted
python src/cellprofiler_tools/converters/vsi_to_ometiff.py /path/to/input_dir /output/dir --dry-run
```

### 3. Extract Green Channel (Recommended for Analysis)

For CellProfiler analysis, you typically only need the green channel. This reduces file size significantly:

```bash
# Extract green channel from converted OME-TIFF
python src/cellprofiler_tools/converters/extract_channel.py /path/to/image.ome.tiff /output/dir --green

# Process entire directory
python src/cellprofiler_tools/converters/extract_channel.py /input/dir /output/dir --green --recursive
```

### 4. Process BigTIFF Files

For large OME-BigTIFF files (`.ome.btf`), extract the green channel directly with compression:

```bash
# Single BigTIFF file - achieves ~8x size reduction
python src/cellprofiler_tools/converters/btf_to_green.py /path/to/image.ome.btf /output/dir

# Process entire directory
python src/cellprofiler_tools/converters/btf_to_green.py /input/dir /output/dir --recursive

# Preview what would be processed
python src/cellprofiler_tools/converters/btf_to_green.py /input/dir /output/dir --dry-run
```

### 5. Batch Processing

```bash
# Sequential processing
python src/cellprofiler_tools/batch/process.py /input/dir /output/dir

# Parallel processing with RAM monitoring
python src/cellprofiler_tools/batch/process_parallel.py /input/dir /output/dir --workers 4
```

### 6. Verify Metadata Preservation

```bash
# Check a single file
python src/cellprofiler_tools/analysis/verify_metadata.py /output/dir/image.ome.tiff

# Check all files in directory
python src/cellprofiler_tools/analysis/verify_metadata.py /output/dir
```

### 7. Calculate IOD (Post-CellProfiler Analysis)

After running CellProfiler analysis, calculate true IOD for genome size estimation:

```bash
python src/cellprofiler_tools/analysis/calculate_iod.py /path/to/cellprofiler_output.csv /path/to/images/
```

## CellProfiler Pipelines

Two analysis pipelines are included in `pipelines/`:

- **oil_immersion.cppipe**: For nucleus-only images (immersion oil microscopy)
- **brightfield.cppipe**: For nucleus + cytoplasm images (brightfield microscopy)

See `docs/PIPELINE_DOCUMENTATION.md` for detailed pipeline documentation.

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

## Safety Features

The scripts include safety checks to prevent accidental overwrites:
- Refuses to write to the same filesystem as source
- Source drive should be mounted read-only
- All processing creates new files in the destination directory

## Requirements

- Python 3.12+
- Bio-Formats command line tools (bftools)
- Dependencies managed via uv (see pyproject.toml)
