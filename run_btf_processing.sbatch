#!/bin/bash
#SBATCH --job-name=btf_processor
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH -N 1
#SBATCH -o btf_processing_%j.out
#SBATCH -e btf_processing_%j.err

# BTF File Processing Job
# ======================
# This script processes BTF files on an HPC cluster, extracting green channels
# and creating memory-efficient tiles for CellProfiler analysis.

# Job information
echo "=========================================="
echo "BTF Processing Job Started"
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
CORRUPTED_OME_WORKAROUND="${CORRUPTED_OME_WORKAROUND:-true}"

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

# Run the processing pipeline
echo "Starting BTF processing pipeline..."
uv run python hpc_btf_processor.py \
    --input-dir "$INPUT_DIR" \
    --output-dir "$OUTPUT_DIR" \
    --tile-size "$TILE_SIZE" \
    --green-channel "$GREEN_CHANNEL" \
    --log-level "$LOG_LEVEL" \
    --corrupted-ome-workaround "$CORRUPTED_OME_WORKAROUND"

# Check exit status
if [ $? -eq 0 ]; then
    echo "=========================================="
    echo "BTF Processing Job Completed Successfully"
    echo "End time: $(date)"
    echo "Output directory: $OUTPUT_DIR"
    echo "=========================================="
    
    # Display summary
    if [ -f "$OUTPUT_DIR/processing_summary.json" ]; then
        echo "Processing Summary:"
        cat "$OUTPUT_DIR/processing_summary.json" | uv run python -m json.tool
    fi
    
    # Show output structure
    echo "Output structure:"
    find "$OUTPUT_DIR" -type f -name "*.tif" | head -10
    echo "..."
    echo "Total tiles created: $(find "$OUTPUT_DIR" -type f -name "*.tif" | wc -l)"
    
else
    echo "=========================================="
    echo "BTF Processing Job Failed"
    echo "End time: $(date)"
    echo "Check logs for details"
    echo "=========================================="
    exit 1
fi

echo "Job completed at $(date)"
