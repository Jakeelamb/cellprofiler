#!/bin/bash
#SBATCH --job-name=btf_2pass_processor
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=your.email@domain.com
#SBATCH -N 1
#SBATCH -o btf_2pass_processing_%j.out
#SBATCH -e btf_2pass_processing_%j.err

# Two-Pass BTF File Processing Job
# ================================
# This script processes BTF files using a two-pass approach:
# 1. Extract green channel to a new file (reduces memory usage by ~70%)
# 2. Tile the smaller green-only file

# Job information
echo "=========================================="
echo "Two-Pass BTF Processing Job Started"
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $SLURM_NODELIST"
echo "CPUs: $SLURM_CPUS_PER_TASK"
echo "Memory: $SLURM_MEM_PER_NODE MB"
echo "Start time: $(date)"
echo "=========================================="

# Set up environment (if modules are available)
if command -v module &> /dev/null; then
    module purge
    module load python/3.12
    module load gcc/11.2.0
else
    echo "Module system not available, using system Python"
fi

# Install UV if not available
if ! command -v uv &> /dev/null; then
    echo "Installing UV..."
    curl -LsSf https://astral.sh/uv/install.sh | sh
    export PATH="$HOME/.cargo/bin:$PATH"
fi

# Set up working directory
WORK_DIR="/home/$USER/cellprofiler"
cd "$WORK_DIR"

# Install dependencies with UV
echo "Installing dependencies with UV..."
uv sync

# Create necessary directories
mkdir -p input output logs

# Configuration
INPUT_DIR="${INPUT_DIR:-./input}"
OUTPUT_DIR="${OUTPUT_DIR:-./output}"
TILE_SIZE="${TILE_SIZE:-2048}"
GREEN_CHANNEL="${GREEN_CHANNEL:-1}"
LOG_LEVEL="${LOG_LEVEL:-INFO}"
CLEANUP_INTERMEDIATE="${CLEANUP_INTERMEDIATE:-true}"

# Validate input directory
if [ ! -d "$INPUT_DIR" ]; then
    echo "ERROR: Input directory $INPUT_DIR does not exist!"
    exit 1
fi

# Check for BTF files
BTF_COUNT=$(find "$INPUT_DIR" -name "*.btf" -o -name "*.ome.btf" | wc -l)
if [ "$BTF_COUNT" -eq 0 ]; then
    echo "ERROR: No BTF files found in $INPUT_DIR"
    exit 1
fi

echo "Found $BTF_COUNT BTF files to process"

# Run the two-pass processing pipeline
echo "Starting two-pass BTF processing pipeline..."
uv run python hpc_btf_processor_2pass.py \
    --input-dir "$INPUT_DIR" \
    --output-dir "$OUTPUT_DIR" \
    --tile-size "$TILE_SIZE" \
    --green-channel "$GREEN_CHANNEL" \
    --log-level "$LOG_LEVEL" \
    --cleanup-intermediate "$CLEANUP_INTERMEDIATE"

# Check exit status
if [ $? -eq 0 ]; then
    echo "=========================================="
    echo "Two-Pass BTF Processing Job Completed Successfully"
    echo "End time: $(date)"
    echo "Output directory: $OUTPUT_DIR"
    echo "=========================================="
    
    # Display summary
    if [ -f "$OUTPUT_DIR/two_pass_processing_summary.json" ]; then
        echo "Processing Summary:"
        cat "$OUTPUT_DIR/two_pass_processing_summary.json" | uv run python -m json.tool
    fi
    
    # Show output structure
    echo "Output structure:"
    find "$OUTPUT_DIR" -type f -name "*.tif" | head -10
    echo "..."
    echo "Total tiles created: $(find "$OUTPUT_DIR" -type f -name "*.tif" | wc -l)"
    
else
    echo "=========================================="
    echo "Two-Pass BTF Processing Job Failed"
    echo "End time: $(date)"
    echo "Check logs for details"
    echo "=========================================="
    exit 1
fi

echo "Job completed at $(date)"
