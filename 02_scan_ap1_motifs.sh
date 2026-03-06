#!/bin/bash
# ============================================================
# STEP 2: Scan promoter sequences for AP-1 motifs using FIMO
# Motif: FOS::JUN (JASPAR MA0099.3)
# ============================================================

set -euo pipefail

# ============================================================
# UPDATE THESE PATHS FOR EACH CELL LINE
# ============================================================
CELL="SAMPLE1"   # <-- # Change CELL="SAMPLE1" / "SAMPLE2" / "SAMPLE3" in each script

BASE="/mnt/Z/ATAC_PPG_plus_MF022UCD65_MF032UCD12/RAMO_Analysis"
ANALYSIS_DIR="${BASE}/Plot_AP-1_sites/${CELL}_AP-1"

OUTDIR="${ANALYSIS_DIR}/results/motifs"
PROMOTER_FA="${ANALYSIS_DIR}/results/promoter_sequences.fa"
MOTIF_DB="${ANALYSIS_DIR}/motifs/AP1_FOSJUN.meme"

CONSENSUS_PEAKS="${BASE}/Plot_AP-1_sites/consensus_narrow_peaks/consensus_peaks_narrow_filtered.bed"

# Per-condition peak files (merged replicates)
CM_PEAKS="${ANALYSIS_DIR}/results/${CELL}_CM_peaks.bed"
CONTROL_PEAKS="${ANALYSIS_DIR}/results/${CELL}_CONTROL_peaks.bed"

# narrowPeak directory for this cell line
NARROW_PEAK_DIR="${BASE}/${CELL}_analysis_03-06-2025/${CELL}_results_nf-core/bwa/merged_library/macs2/narrow_peak"

mkdir -p "$OUTDIR" "${ANALYSIS_DIR}/motifs"

# ============================================================
# DOWNLOAD AP-1 MOTIF FROM JASPAR
# FOS::JUN  MA0099.3  (canonical AP-1 heterodimer)
# Also include FOSL2::JUN, FOSL1::JUNB for broader coverage
# ============================================================

echo "Downloading AP-1 motifs from JASPAR..."

# Create MEME format motif file for FIMO
# MA0099.3 FOS::JUN - canonical AP-1 motif (TGASTCA)
cat > "$MOTIF_DB" << 'MEME_EOF'
MEME version 4

ALPHABET= ACGT

strands: + -

Background letter frequencies
A 0.25 C 0.25 G 0.25 T 0.25

MOTIF MA0099.3 FOS::JUN
letter-probability matrix: alength= 4 w= 11 nsites= 7 E= 0
0.014493 0.014493 0.014493 0.956522
0.014493 0.014493 0.956522 0.014493
0.956522 0.014493 0.014493 0.014493
0.014493 0.014493 0.014493 0.956522
0.014493 0.014493 0.014493 0.956522
0.014493 0.014493 0.014493 0.956522
0.014493 0.956522 0.014493 0.014493
0.956522 0.014493 0.014493 0.014493
0.014493 0.014493 0.956522 0.014493
0.014493 0.014493 0.014493 0.956522
0.956522 0.014493 0.014493 0.014493

MOTIF MA0476.1 FOS::JUND
letter-probability matrix: alength= 4 w= 10 nsites= 20 E= 0
0.050000 0.050000 0.050000 0.850000
0.050000 0.050000 0.850000 0.050000
0.850000 0.050000 0.050000 0.050000
0.050000 0.050000 0.050000 0.850000
0.050000 0.050000 0.050000 0.850000
0.050000 0.050000 0.050000 0.850000
0.050000 0.850000 0.050000 0.050000
0.850000 0.050000 0.050000 0.050000
0.050000 0.050000 0.850000 0.050000
0.850000 0.050000 0.050000 0.050000
MEME_EOF

echo "Motif file created: $MOTIF_DB"

# ============================================================
# ALTERNATIVELY: Download full JASPAR 2024 CORE database
# Uncomment to use full database instead
# ============================================================
# wget -O ./motifs/JASPAR2024_CORE_vertebrates.meme \
#   "https://jaspar.elixir.no/download/data/2024/CORE/JASPAR2024_CORE_vertebrates_non-redundant_pfms_meme.txt"
# MOTIF_DB="./motifs/JASPAR2024_CORE_vertebrates.meme"

# ============================================================
# RUN FIMO MOTIF SCANNING
# ============================================================

echo "Running FIMO motif scan..."

fimo \
    --oc "$OUTDIR/fimo_output" \
    --thresh 1e-3 \
    "$MOTIF_DB" \
    "$PROMOTER_FA"

echo "FIMO complete. Results in: $OUTDIR/fimo_output/"

# ============================================================
# CONVERT FIMO OUTPUT TO BED FORMAT
# FIXED: FIMO outputs chr directly in $3, no coordinate parsing needed
# FIXED: Filter only valid chr lines to remove garbage lines
# ============================================================

echo "Converting FIMO output to BED..."

# fimo.tsv columns: motif_id, motif_alt_id, sequence_name, start, stop, strand, score, p-value, q-value, matched_sequence
awk 'NR>1 && $1!="#" && $3~/^chr/ {
    print $3"\t"$4"\t"$5"\t"$1"_"$3"\t"$8"\t"$6
}' "$OUTDIR/fimo_output/fimo.tsv" \
    | sort -k1,1 -k2,2n \
    > "$OUTDIR/AP1_motif_sites.bed"

echo "AP-1 motif BED file: $OUTDIR/AP1_motif_sites.bed"
echo "Total AP-1 motif sites: $(wc -l < $OUTDIR/AP1_motif_sites.bed)"
echo ""
echo "Motifs per gene:"
cut -f1 "$OUTDIR/AP1_motif_sites.bed" | sort | uniq -c

# ============================================================
# MERGE REPLICATES PER CONDITION AND INTERSECT
# ============================================================

echo ""
echo "Merging replicates per condition..."

# Merge CM replicates
cat "${NARROW_PEAK_DIR}/${CELL}_CM_REP1.mLb.clN.sorted_peaks.narrowPeak" \
    "${NARROW_PEAK_DIR}/${CELL}_CM_REP2.mLb.clN.sorted_peaks.narrowPeak" \
    "${NARROW_PEAK_DIR}/${CELL}_CM_REP3.mLb.clN.sorted_peaks.narrowPeak" \
    | sort -k1,1 -k2,2n | bedtools merge -i - > "$CM_PEAKS"

# Merge CONTROL replicates
cat "${NARROW_PEAK_DIR}/${CELL}_CONTROL_REP1.mLb.clN.sorted_peaks.narrowPeak" \
    "${NARROW_PEAK_DIR}/${CELL}_CONTROL_REP2.mLb.clN.sorted_peaks.narrowPeak" \
    "${NARROW_PEAK_DIR}/${CELL}_CONTROL_REP3.mLb.clN.sorted_peaks.narrowPeak" \
    | sort -k1,1 -k2,2n | bedtools merge -i - > "$CONTROL_PEAKS"

echo "CM peaks: $(wc -l < $CM_PEAKS)"
echo "CONTROL peaks: $(wc -l < $CONTROL_PEAKS)"

# ============================================================
# INTERSECT AP-1 MOTIFS WITH PEAKS PER CONDITION
# ============================================================

echo ""
echo "=== AP-1 motifs in open chromatin per gene ==="
echo ""
echo "--- CM (Conditioned Media) ---"
bedtools intersect \
    -a "$OUTDIR/AP1_motif_sites.bed" \
    -b "$CM_PEAKS" \
    -wa > "$OUTDIR/AP1_motifs_in_CM_peaks.bed"
cut -f1 "$OUTDIR/AP1_motifs_in_CM_peaks.bed" | sort | uniq -c

echo ""
echo "--- CONTROL (Starvation) ---"
bedtools intersect \
    -a "$OUTDIR/AP1_motif_sites.bed" \
    -b "$CONTROL_PEAKS" \
    -wa > "$OUTDIR/AP1_motifs_in_CONTROL_peaks.bed"
cut -f1 "$OUTDIR/AP1_motifs_in_CONTROL_peaks.bed" | sort | uniq -c

# ============================================================
# INTERSECT WITH CONSENSUS PEAKS
# ============================================================

echo ""
echo "Intersecting with consensus peaks..."
if [ -f "$CONSENSUS_PEAKS" ]; then
    bedtools intersect \
        -a "$OUTDIR/AP1_motif_sites.bed" \
        -b "$CONSENSUS_PEAKS" \
        -wa -wb \
        > "$OUTDIR/AP1_motifs_in_ATAC_peaks.bed"
    echo "AP-1 motifs overlapping consensus ATAC peaks: $(wc -l < $OUTDIR/AP1_motifs_in_ATAC_peaks.bed)"
else
    echo "WARNING: Consensus peaks file not found at $CONSENSUS_PEAKS"
fi

echo ""
echo "DONE. Next step: run 03_compute_atac_signal.py"
