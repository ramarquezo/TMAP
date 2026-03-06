#!/bin/bash
# ============================================================
# MACS2 Narrow Peak Calling Loop
# Cell lines: UCD04, UCD12, UCD65
# ============================================================

set -euo pipefail

# Base path (update if mount point differs on your HPC)
BASE="/mnt/Z/ATAC_PPG_plus_MF022UCD65_MF032UCD12/RAMO_Analysis"

# Cell lines to process
CELL_LINES=("SAMPLE1" "SAMPLE2" "SAMPLE3") # Change cell lines

for CELL in "${CELL_LINES[@]}"; do

    echo "========================================"
    echo " Processing cell line: $CELL"
    echo "========================================"

    # Input BAM directory
    BAM_DIR="${BASE}/${CELL}_analysis_03-06-2025/${CELL}_results_nf-core/bwa/merged_library"

    # Output narrow peak directory
    OUTDIR="${BASE}/${CELL}_analysis_03-06-2025/${CELL}_results_nf-core/bwa/merged_library/macs2/narrow_peak"

    # Create output directory if it doesn't exist
    mkdir -p "$OUTDIR"

    echo "  BAM dir:  $BAM_DIR"
    echo "  Out dir:  $OUTDIR"
    echo ""

    # Loop over all BAM files for this cell line
    for BAM in "$BAM_DIR"/*.bam; do

        # Skip if no BAM files found
        [ -e "$BAM" ] || { echo "  No BAM files found in $BAM_DIR"; continue; }

        # Get sample name from filename (strip path and .bam extension)
        SAMPLE=$(basename "$BAM" .bam)

        echo "  Running MACS2 on: $SAMPLE"

        macs2 callpeak \
            -t "$BAM" \
            -f BAMPE \
            -n "${CELL}_${SAMPLE}" \
            --outdir "$OUTDIR" \
            -g hs \
            --nomodel \
            --shift -100 \
            --extsize 200 \
            --keep-dup all \
            -q 0.05 \
            2>&1 | tee "$OUTDIR/${CELL}_${SAMPLE}_macs2.log"

        echo "  Done: ${CELL}_${SAMPLE}"
        echo ""

    done

done

echo "========================================"
echo " ALL DONE!"
echo " Narrow peaks saved per cell line in:"
for CELL in "${CELL_LINES[@]}"; do
    echo "  ${BASE}/${CELL}_analysis_03-06-2025/${CELL}_results_nf-core/bwa/merged_library/macs2/narrow_peak"
done
echo "========================================"
