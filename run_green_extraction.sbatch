#!/bin/bash
#SBATCH --job-name=green_extractor
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=your.email@domain.com
#SBATCH -N 1
#SBATCH -o green_extraction_%j.out
#SBATCH -e green_extraction_%j.err

# Green Channel Extraction Job
# ============================
# This script extracts only the green channel from BTF files,
# creating clean TIFF files without OME metadata corruption.

# Job information
echo "=========================================="
echo "Green Channel Extraction Job Started"
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

# Start memory monitoring in background
echo "Starting memory monitoring..."
uv run python monitor_memory.py --process-name extract_green_only.py --log-file memory_monitor.log &
MONITOR_PID=$!

# Run the green channel extraction
echo "Starting green channel extraction..."
uv run python extract_green_only.py \
    --input-dir "$INPUT_DIR" \
    --output-dir "$OUTPUT_DIR" \
    --green-channel "$GREEN_CHANNEL" \
    --log-level "$LOG_LEVEL"

# Stop memory monitoring
echo "Stopping memory monitoring..."
kill $MONITOR_PID 2>/dev/null || true

# Check exit status
if [ $? -eq 0 ]; then
    echo "=========================================="
    echo "Green Channel Extraction Job Completed Successfully"
    echo "End time: $(date)"
    echo "Output directory: $OUTPUT_DIR"
    echo "=========================================="
    
    # Display summary
    if [ -f "$OUTPUT_DIR/green_extraction_summary.json" ]; then
        echo "Extraction Summary:"
        cat "$OUTPUT_DIR/green_extraction_summary.json" | uv run python -m json.tool
    fi
    
    # Show output files
    echo "Green channel files created:"
    find "$OUTPUT_DIR" -name "*_green.tif" -exec ls -lh {} \;
    
else
    echo "=========================================="
    echo "Green Channel Extraction Job Failed"
    echo "End time: $(date)"
    echo "Check logs for details"
    echo "=========================================="
    exit 1
fi

echo "Job completed at $(date)"
