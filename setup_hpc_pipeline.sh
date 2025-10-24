#!/bin/bash
# Setup script for BTF HPC Pipeline
# =================================
# This script sets up the BTF processing pipeline on an HPC cluster

set -e  # Exit on any error

echo "=========================================="
echo "BTF HPC Pipeline Setup"
echo "=========================================="

# Check if we're in the right directory
if [ ! -f "hpc_btf_processor.py" ]; then
    echo "❌ Error: Please run this script from the project directory"
    exit 1
fi

# Make scripts executable
echo "Making scripts executable..."
chmod +x hpc_btf_processor.py
chmod +x validate_output.py
chmod +x pipeline_manager.py
chmod +x run_btf_processing.sbatch

# Create necessary directories
echo "Creating directories..."
mkdir -p input output logs

# Create default configuration if it doesn't exist
if [ ! -f "config.yaml" ]; then
    echo "Creating default configuration..."
    cp config_example.yaml config.yaml
    echo "✅ Created config.yaml from template"
else
    echo "✅ Configuration file already exists"
fi

# Check Python environment
echo "Checking Python environment..."
if command -v python3 &> /dev/null; then
    PYTHON_VERSION=$(python3 --version | cut -d' ' -f2)
    echo "✅ Python $PYTHON_VERSION found"
else
    echo "❌ Python 3 not found. Please install Python 3.12+"
    exit 1
fi

# Check if we're on a cluster with SLURM
if command -v sbatch &> /dev/null; then
    echo "✅ SLURM detected - HPC environment ready"
    SLURM_AVAILABLE=true
else
    echo "⚠️  SLURM not detected - local processing only"
    SLURM_AVAILABLE=false
fi

# Install Python dependencies
echo "Installing Python dependencies..."
if command -v pip3 &> /dev/null; then
    pip3 install --user imagecodecs numpy tifffile tqdm PyYAML
    echo "✅ Dependencies installed"
else
    echo "❌ pip3 not found. Please install pip"
    exit 1
fi

# Test the pipeline
echo "Testing pipeline components..."
python3 -c "import tifffile, numpy, tqdm, yaml; print('✅ All imports successful')"

# Create example usage script
cat > run_example.sh << 'EOF'
#!/bin/bash
# Example usage script for BTF processing

echo "BTF Processing Example"
echo "====================="

# Setup (if not already done)
echo "1. Setting up environment..."
python3 pipeline_manager.py setup --input-dir ./input --output-dir ./output

# Copy your BTF files to ./input directory
echo "2. Please copy your BTF files to ./input directory"
echo "   Example: cp /path/to/your/files/*.btf ./input/"

# Process files
echo "3. Processing files..."
python3 hpc_btf_processor.py --config config.yaml

# Validate results
echo "4. Validating results..."
python3 validate_output.py --output-dir ./output

echo "✅ Example completed!"
EOF

chmod +x run_example.sh

# Display setup summary
echo ""
echo "=========================================="
echo "Setup Complete!"
echo "=========================================="
echo ""
echo "Next steps:"
echo "1. Copy your BTF files to ./input directory"
echo "2. Edit config.yaml if needed"
echo "3. Run the pipeline:"
echo "   - Local: python3 hpc_btf_processor.py --config config.yaml"
echo "   - HPC:   sbatch run_btf_processing.sbatch"
echo "   - Or use: ./run_example.sh"
echo ""
echo "Available commands:"
echo "  python3 pipeline_manager.py --help"
echo "  python3 hpc_btf_processor.py --help"
echo "  python3 validate_output.py --help"
echo ""

if [ "$SLURM_AVAILABLE" = true ]; then
    echo "HPC Environment:"
    echo "  - SLURM job script: run_btf_processing.sbatch"
    echo "  - Submit job: sbatch run_btf_processing.sbatch"
    echo "  - Monitor: python3 pipeline_manager.py monitor --job-id <JOB_ID>"
fi

echo ""
echo "For more information, see README.md"
echo "=========================================="
