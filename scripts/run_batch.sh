#!/bin/bash
# Robust batch processing script - runs reliably unattended
# Processes BTF files then VSI files, skipping already completed ones

set -e  # Exit on error (but we handle errors per-file)

LOG_FILE="batch_run_$(date +%Y%m%d_%H%M%S).log"
INPUT_DIR="/run/media/jake/easystore"
OUTPUT_DIR="data"
WORKERS=8  # Conservative for reliability - USB I/O is the bottleneck anyway

log() {
    echo "[$(date '+%H:%M:%S')] $1" | tee -a "$LOG_FILE"
}

log "=== Starting batch processing ==="
log "Input: $INPUT_DIR"
log "Output: $OUTPUT_DIR"
log "Workers: $WORKERS"
log "Log file: $LOG_FILE"

mkdir -p "$OUTPUT_DIR"

# Count files
BTF_COUNT=$(find "$INPUT_DIR" -name "*.btf" -o -name "*.ome.btf" 2>/dev/null | wc -l)
VSI_COUNT=$(find "$INPUT_DIR" -name "*.vsi" -o -name "*.VSI" 2>/dev/null | wc -l)
log "Found $BTF_COUNT BTF files and $VSI_COUNT VSI files"

# Function to process a single BTF file
process_btf() {
    local btf_file="$1"
    local basename=$(basename "$btf_file" .btf)
    basename=${basename%.ome}  # Remove .ome if present
    local output_file="$OUTPUT_DIR/${basename}_green.ome.tiff"

    # Skip if already exists
    if [[ -f "$output_file" ]] && [[ $(stat -c%s "$output_file") -gt 1000 ]]; then
        echo "SKIP: $basename (exists)"
        return 0
    fi

    echo "PROCESSING: $basename"
    if uv run python btf_to_green.py "$btf_file" "$OUTPUT_DIR" --channel 1 2>/dev/null; then
        if [[ -f "$output_file" ]]; then
            local size=$(du -h "$output_file" | cut -f1)
            echo "OK: $basename ($size)"
            return 0
        fi
    fi
    echo "FAILED: $basename"
    return 1
}

export -f process_btf
export OUTPUT_DIR

# Process BTF files
log "=== Processing BTF files ==="
COMPLETED_BTF=0
FAILED_BTF=0

find "$INPUT_DIR" -name "*.btf" -o -name "*.ome.btf" 2>/dev/null | sort | while read btf_file; do
    basename=$(basename "$btf_file" .btf)
    basename=${basename%.ome}
    output_file="$OUTPUT_DIR/${basename}_green.ome.tiff"

    if [[ -f "$output_file" ]] && [[ $(stat -c%s "$output_file") -gt 1000 ]]; then
        log "SKIP: $basename"
        continue
    fi

    log "Processing: $basename"
    if timeout 1800 uv run python btf_to_green.py "$btf_file" "$OUTPUT_DIR" --channel 1 >> "$LOG_FILE" 2>&1; then
        if [[ -f "$output_file" ]]; then
            size=$(du -h "$output_file" | cut -f1)
            log "OK: $basename ($size)"
            ((COMPLETED_BTF++)) || true
        else
            log "FAILED: $basename (no output)"
            ((FAILED_BTF++)) || true
        fi
    else
        log "FAILED: $basename (error)"
        ((FAILED_BTF++)) || true
    fi

    # Log progress and memory every 5 files
    completed=$(ls "$OUTPUT_DIR"/*.tiff 2>/dev/null | wc -l)
    mem=$(free -h | grep Mem | awk '{print $4}')
    log "Progress: $completed files completed, $mem free RAM"
done

# Process VSI files (if bfconvert exists)
BFCONVERT="$HOME/bin/bftools/bftools/bfconvert"
if [[ -x "$BFCONVERT" ]]; then
    log "=== Processing VSI files ==="

    TEMP_DIR=$(mktemp -d)
    trap "rm -rf $TEMP_DIR" EXIT

    find "$INPUT_DIR" -name "*.vsi" -o -name "*.VSI" 2>/dev/null | sort | while read vsi_file; do
        basename=$(basename "$vsi_file" .vsi)
        basename=${basename%.VSI}
        output_file="$OUTPUT_DIR/${basename}_green.ome.tiff"

        if [[ -f "$output_file" ]] && [[ $(stat -c%s "$output_file") -gt 1000 ]]; then
            log "SKIP: $basename"
            continue
        fi

        log "Processing VSI: $basename"
        temp_tiff="$TEMP_DIR/${basename}.ome.tiff"

        # Convert VSI to TIFF
        if timeout 3600 "$BFCONVERT" -bigtiff -compression LZW -tilex 512 -tiley 512 -series 0 "$vsi_file" "$temp_tiff" >> "$LOG_FILE" 2>&1; then
            # Extract green channel
            if timeout 1800 uv run python btf_to_green.py "$temp_tiff" "$OUTPUT_DIR" --channel 1 >> "$LOG_FILE" 2>&1; then
                if [[ -f "$output_file" ]]; then
                    size=$(du -h "$output_file" | cut -f1)
                    log "OK: $basename ($size)"
                else
                    log "FAILED: $basename (no output after green extraction)"
                fi
            else
                log "FAILED: $basename (green extraction error)"
            fi
            rm -f "$temp_tiff"
        else
            log "FAILED: $basename (bfconvert error)"
        fi

        # Log progress
        completed=$(ls "$OUTPUT_DIR"/*.tiff 2>/dev/null | wc -l)
        mem=$(free -h | grep Mem | awk '{print $4}')
        log "Progress: $completed files completed, $mem free RAM"
    done
else
    log "WARNING: bfconvert not found, skipping VSI files"
fi

# Final summary
TOTAL_COMPLETED=$(ls "$OUTPUT_DIR"/*.tiff 2>/dev/null | wc -l)
log "=== COMPLETE ==="
log "Total files created: $TOTAL_COMPLETED"
log "Output directory: $OUTPUT_DIR"
