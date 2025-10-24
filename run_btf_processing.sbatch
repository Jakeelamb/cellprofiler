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

# Set up environment
module purge
module load python/3.12
module load gcc/11.2.0

# Create virtual environment if it doesn't exist
VENV_DIR="/home/$USER/btf_processing_env"
if [ ! -d "$VENV_DIR" ]; then
    echo "Creating virtual environment..."
    python -m venv "$VENV_DIR"
fi

# Activate virtual environment
source "$VENV_DIR/bin/activate"

# Install/update dependencies
echo "Installing dependencies..."
pip install --upgrade pip
pip install imagecodecs>=2025.8.2 numpy>=2.3.4 tifffile>=2025.10.16 tqdm>=4.66.0 PyYAML>=6.0

# Set up working directory
WORK_DIR="/home/$USER/Projects/cellprofiler_test"
cd "$WORK_DIR"

# Create necessary directories
mkdir -p input output logs

# Configuration
INPUT_DIR="${INPUT_DIR:-./input}"
OUTPUT_DIR="${OUTPUT_DIR:-./output}"
TILE_SIZE="${TILE_SIZE:-2048}"
GREEN_CHANNEL="${GREEN_CHANNEL:-1}"
LOG_LEVEL="${LOG_LEVEL:-INFO}"

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
python hpc_btf_processor.py \
    --input-dir "$INPUT_DIR" \
    --output-dir "$OUTPUT_DIR" \
    --tile-size "$TILE_SIZE" \
    --green-channel "$GREEN_CHANNEL" \
    --log-level "$LOG_LEVEL"

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
        cat "$OUTPUT_DIR/processing_summary.json" | python -m json.tool
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

# Deactivate virtual environment
deactivate

echo "Job completed at $(date)"
