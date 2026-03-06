#!/bin/bash
# ============================================================
# Create Consensus Narrow Peaks from MACS2 output
# Merges all narrowPeak files across:
#   - 3 cell lines (UCD04, UCD12, UCD65)
#   - 2 conditions (CM, CONTROL)
#   - 3 replicates each
# ============================================================

set -euo pipefail

BASE="/mnt/Z/ATAC_PPG_plus_MF022UCD65_MF032UCD12/RAMO_Analysis"
CELL_LINES=("SAMPLE1" "SAMPLE2" "SAMPLE3") # Change cell lines

OUTDIR="/mnt/Z/ATAC_PPG_plus_MF022UCD65_MF032UCD12/RAMO_Analysis/Plot_AP-1_sites/consensus_narrow_peaks"
mkdir -p "$OUTDIR"

MERGED_ALL="${OUTDIR}/all_narrowpeaks_merged.bed"
CONSENSUS="${OUTDIR}/consensus_peaks_narrow.bed"

echo "=============================================="
echo " Creating consensus narrow peaks"
echo "=============================================="

# ============================================================
# STEP 1: Collect and concatenate all narrowPeak files
# ============================================================

echo ""
echo "[1/4] Collecting all narrowPeak files..."

# Remove previous merged file if exists
rm -f "$MERGED_ALL"

PEAK_COUNT=0
for CELL in "${CELL_LINES[@]}"; do

    PEAK_DIR="${BASE}/${CELL}_analysis_03-06-2025/${CELL}_results_nf-core/bwa/merged_library/macs2/narrow_peak"

    echo "  Looking in: $PEAK_DIR"

    for PEAK in "$PEAK_DIR"/*.narrowPeak; do

        [ -e "$PEAK" ] || { echo "  WARNING: No narrowPeak files found in $PEAK_DIR"; continue; }

        echo "  Found: $(basename $PEAK)"

        # Extract chr, start, end, score, qvalue columns only
        # narrowPeak format: chr start end name score strand signalValue pValue qValue summit
        awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10}' "$PEAK" \
            >> "$MERGED_ALL"

        PEAK_COUNT=$((PEAK_COUNT + 1))
    done
done

echo ""
echo "  Total narrowPeak files collected: $PEAK_COUNT"
echo "  Total peaks before merging: $(wc -l < $MERGED_ALL)"

# ============================================================
# STEP 2: Sort the combined peak file
# ============================================================

echo ""
echo "[2/4] Sorting peaks..."

sort -k1,1 -k2,2n "$MERGED_ALL" > "${OUTDIR}/all_narrowpeaks_sorted.bed"

# ============================================================
# STEP 3: Merge overlapping peaks (bedtools merge)
# A peak is kept if present in at least 1 sample
# ============================================================

echo ""
echo "[3/4] Merging overlapping peaks with bedtools..."

bedtools merge \
    -i "${OUTDIR}/all_narrowpeaks_sorted.bed" \
    -d 0 \
    > "$CONSENSUS"

echo "  Total consensus peaks after merging: $(wc -l < $CONSENSUS)"

# ============================================================
# STEP 4: Filter out blacklisted regions (hg38)
# ============================================================

echo ""
echo "[4/4] Filtering ENCODE blacklist regions..."

BLACKLIST="/home/ramarquezo/hg38-blacklist.v2.bed"

# Use the same blacklist from your nf-core run
if [ -f "$BLACKLIST" ]; then
    bedtools intersect \
        -a "$CONSENSUS" \
        -b "$BLACKLIST" \
        -v \
        > "${OUTDIR}/consensus_peaks_narrow_filtered.bed"

    echo "  Peaks after blacklist filtering: $(wc -l < ${OUTDIR}/consensus_peaks_narrow_filtered.bed)"
    FINAL="${OUTDIR}/consensus_peaks_narrow_filtered.bed"
else
    echo "  WARNING: Blacklist not found at $BLACKLIST"
    echo "  Skipping blacklist filtering — using unfiltered consensus"
    FINAL="$CONSENSUS"
fi

# ============================================================
# SUMMARY
# ============================================================

echo ""
echo "=============================================="
echo " DONE!"
echo " Final consensus peak file:"
echo " $FINAL"
echo ""
echo " Use this path in 02_scan_ap1_motifs.sh:"
echo " CONSENSUS_PEAKS=\"$FINAL\""
echo "=============================================="
