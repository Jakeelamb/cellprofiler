# BTF Processing Pipeline for HPC Clusters

A robust pipeline for processing large BTF (Bio-Formats TIFF) files on HPC clusters, designed specifically for CellProfiler analysis. This pipeline extracts green channels and creates memory-efficient tiles from large microscopy images.

## Features

- **Memory-efficient processing**: Handles extremely large BTF files without loading entire images into memory
- **HPC-optimized**: Designed for SLURM-based HPC clusters with configurable resource allocation
- **Robust file organization**: Maintains clear file naming conventions and directory structure
- **Comprehensive validation**: Built-in tools to verify processing results
- **Flexible configuration**: YAML-based configuration system with command-line overrides
- **Progress monitoring**: Real-time progress tracking and job monitoring
- **Error handling**: Comprehensive error handling and logging

## Quick Start

### 1. Setup Environment

```bash
# Setup processing directories
python pipeline_manager.py setup --input-dir /path/to/btf/files --output-dir /path/to/output

# Or with custom config
python pipeline_manager.py setup --input-dir /path/to/btf/files --output-dir /path/to/output --config my_config.yaml
```

### 2. Submit Job to HPC Cluster

```bash
# Submit with default configuration
python pipeline_manager.py submit --config config.yaml

# Monitor job progress
python pipeline_manager.py monitor --job-id 12345
```

### 3. Validate Results

```bash
# Validate processed output
python pipeline_manager.py validate --output-dir /path/to/output
```

## File Structure

```
cellprofiler_test/
├── hpc_btf_processor.py      # Main processing pipeline
├── run_btf_processing.sbatch # SLURM job script
├── pipeline_manager.py       # Management utilities
├── validate_output.py        # Output validation tool
├── config_example.yaml       # Example configuration
├── pyproject.toml           # Python dependencies
└── README.md                # This file
```

## Configuration

The pipeline uses YAML configuration files. See `config_example.yaml` for all available options:

```yaml
# Basic configuration
input_dir: "./input"          # Directory containing BTF files
output_dir: "./output"        # Directory for processed tiles
tile_size: 2048              # Size of each tile (pixels)
green_channel: 1             # Index of green channel (0=Red, 1=Green, 2=Blue)

# Logging
log_level: "INFO"            # DEBUG, INFO, WARNING, ERROR
log_file: "btf_processing.log"
```

## Usage Examples

### Command Line Interface

```bash
# Process with command line arguments
python hpc_btf_processor.py \
    --input-dir /data/btf_files \
    --output-dir /data/processed \
    --tile-size 2048 \
    --green-channel 1

# Process with configuration file
python hpc_btf_processor.py --config config.yaml
```

### SLURM Job Submission

```bash
# Submit job to cluster
sbatch run_btf_processing.sbatch

# With custom parameters
INPUT_DIR=/data/btf_files OUTPUT_DIR=/data/processed sbatch run_btf_processing.sbatch
```

### Pipeline Management

```bash
# Check job status
python pipeline_manager.py status --job-id 12345

# Monitor job progress
python pipeline_manager.py monitor --job-id 12345 --interval 30

# Show disk usage
python pipeline_manager.py stats --output-dir /data/processed

# Clean up old files
python pipeline_manager.py cleanup --days 7
```

## Output Organization

The pipeline creates a well-organized output structure:

```
output/
├── Process_374_raw/                    # Source file directory
│   ├── green_tiles/                   # Tiles directory
│   │   ├── Process_374_raw_green_tile_0000_0_0_2048x2048.tif
│   │   ├── Process_374_raw_green_tile_0001_0_2048_2048x2048.tif
│   │   └── ...
│   └── Process_374_raw_metadata.json  # Processing metadata
├── processing_summary.json            # Overall processing summary
└── btf_processing.log                # Processing logs
```

### File Naming Convention

Tiles follow this naming pattern:
```
{source_filename}_green_tile_{tile_id:04d}_{y}_{x}_{height}x{width}.tif
```

Example: `Process_374_raw_green_tile_0000_0_0_2048x2048.tif`

## HPC Resource Requirements

### Recommended SLURM Configuration

```bash
#SBATCH --cpus-per-task=8      # 8 CPU cores
#SBATCH --mem=64G              # 64 GB RAM
#SBATCH --time=24:00:00        # 24 hours max runtime
#SBATCH -N 1                   # Single node
```

### Memory Requirements

- **Minimum**: 32 GB RAM for 2048x2048 tiles
- **Recommended**: 64 GB RAM for optimal performance
- **Large files**: 128+ GB RAM for files > 10 GB

### Storage Requirements

- **Input**: Original BTF file size
- **Output**: ~25% of original size (green channel only)
- **Temporary**: Additional 10-20% during processing

## Troubleshooting

### Common Issues

1. **Memory errors**: Reduce `tile_size` or increase allocated memory
2. **File not found**: Check input directory path and file permissions
3. **Permission denied**: Ensure write access to output directory
4. **Job timeout**: Increase `--time` parameter in SLURM script

### Validation

The pipeline includes comprehensive validation:

```bash
# Validate all output files
python validate_output.py --output-dir /path/to/output

# Check specific issues
python validate_output.py --output-dir /path/to/output --log-level DEBUG
```

### Logging

All operations are logged to:
- Console output (stdout)
- Log files (`btf_processing.log`)
- SLURM output files (`btf_processing_*.out`)

## Dependencies

- Python 3.12+
- imagecodecs >= 2025.8.2
- numpy >= 2.3.4
- tifffile >= 2025.10.16
- tqdm >= 4.66.0
- PyYAML >= 6.0

## Performance Tips

1. **Optimize tile size**: Larger tiles reduce file count but increase memory usage
2. **Use fast storage**: Place input/output on high-speed storage (SSD, parallel filesystem)
3. **Monitor memory**: Watch memory usage during processing
4. **Batch processing**: Process multiple files in single job when possible

## Support

For issues or questions:
1. Check the log files for error messages
2. Run validation to identify specific problems
3. Verify file permissions and paths
4. Check SLURM job logs for cluster-specific issues

## License

This project is part of the CellProfiler analysis pipeline.
