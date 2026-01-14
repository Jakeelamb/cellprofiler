#!/bin/bash
# Simple reliable batch processor - processes 4 files at a time
# Run with: nohup ./process_all.sh &

LOG="process_$(date +%Y%m%d_%H%M%S).log"
INPUT="/run/media/jake/easystore"
OUTPUT="data"
JOBS=4  # Number of parallel jobs

exec > >(tee -a "$LOG") 2>&1

echo "=== Started $(date) ==="
echo "Log: $LOG"
echo "Jobs: $JOBS parallel"

mkdir -p "$OUTPUT"

process_file() {
    local f="$1"
    local base=$(basename "$f")
    local stem="${base%.btf}"
    stem="${stem%.ome}"
    local out="$OUTPUT/${stem}_green.ome.tiff"

    if [[ -f "$out" ]] && [[ $(stat -c%s "$out" 2>/dev/null || echo 0) -gt 1000 ]]; then
        echo "[SKIP] $base"
        return 0
    fi

    echo "[START] $base"
    timeout 10800 uv run python btf_to_green.py "$f" "$OUTPUT" --channel 1 2>&1 | grep -v "OME series" || true
    # Clean up any orphaned temp files for this file
    rm -f "$OUTPUT/.tmp_${stem}"*.dat 2>/dev/null
    if [[ -f "$out" ]] && [[ $(stat -c%s "$out" 2>/dev/null || echo 0) -gt 1000000 ]]; then
        echo "[DONE] $base ($(du -h "$out" | cut -f1))"
    else
        echo "[FAIL] $base - no output or too small"
    fi
}
export -f process_file
export OUTPUT

echo "--- BTF FILES ---"
find "$INPUT" -type f \( -name "*.btf" -o -name "*.ome.btf" \) 2>/dev/null | sort | \
    xargs -P $JOBS -I {} bash -c 'process_file "$@"' _ {}

echo "--- VSI FILES ---"
BFCONVERT="$HOME/bin/bftools/bftools/bfconvert"
if [[ -x "$BFCONVERT" ]]; then
    TEMP=$(mktemp -d)
    trap "rm -rf $TEMP" EXIT

    for f in $(find "$INPUT" -type f \( -name "*.vsi" -o -name "*.VSI" \) 2>/dev/null | sort); do
        base=$(basename "$f")
        stem="${base%.vsi}"
        stem="${stem%.VSI}"
        out="$OUTPUT/${stem}_green.ome.tiff"

        if [[ -f "$out" ]] && [[ $(stat -c%s "$out" 2>/dev/null || echo 0) -gt 1000 ]]; then
            echo "[SKIP] $base"
            continue
        fi

        echo "[START VSI] $base"
        tmp="$TEMP/${stem}.ome.tiff"
        if timeout 3600 "$BFCONVERT" -bigtiff -compression LZW -tilex 512 -tiley 512 -series 0 "$f" "$tmp" 2>&1 | tail -5; then
            if timeout 10800 uv run python btf_to_green.py "$tmp" "$OUTPUT" --channel 1 2>&1 | grep -v "OME series"; then
                if [[ -f "$out" ]]; then
                    echo "[DONE] $base ($(du -h "$out" | cut -f1))"
                fi
            fi
            rm -f "$tmp"
        else
            echo "[FAIL VSI] $base"
        fi
    done
fi

echo "=== Completed $(date) ==="
echo "Total files: $(ls "$OUTPUT"/*.tiff 2>/dev/null | wc -l)"
