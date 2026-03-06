#!/bin/bash
# ============================================================
# STEP 1: Extract promoter regions for genes of interest
# Genes: GeneA, GeneB, GeneC, GeneD
# Genome: hg38
# Window: TSS ± 2kb (adjustable)
# ============================================================

set -euo pipefail

# ============================================================
# UPDATE THIS FOR EACH CELL LINE
# ============================================================
CELL="SAMPLE1"   # <-- # Change CELL="SAMPLE1" / "SAMPLE2" / "SAMPLE3" in each script

# --- CONFIGURATION ---
BASE="/mnt/Z/ATAC_PPG_plus_MF022UCD65_MF032UCD12/RAMO_Analysis"
OUTDIR="${BASE}/Plot_AP-1_sites/${CELL}_AP-1/results"
WINDOW=2000       # bp upstream and downstream of TSS

GENOME_FASTA="${BASE}/genome/genome.fa"
GENOME_SIZES="${BASE}/genome/genome.fa.sizes"

mkdir -p "$OUTDIR"

echo "=============================================="
echo " STEP 1: Extract promoter regions"
echo " Cell line: $CELL"
echo " Window: TSS ± ${WINDOW}bp"
echo "=============================================="

# --- GENE COORDINATES (hg38, MANE Select TSS positions) ---
# Verified from UCSC GENCODE V49 MANE Select transcripts
# Format: chr  TSS  TSS+1  gene_name  .  strand

cat > "$OUTDIR/genes_of_interest.bed" << 'EOF'
chr12	12891562	12891563	GeneA	.	+
chr10	102395705	102395706	GeneB	.	+
chr19	45001464	45001465	GeneC	.	+
chr1	150876599	150876600	GeneD	.	-
EOF

echo "Gene TSS coordinates written."

# --- EXTEND TO PROMOTER WINDOWS (TSS ± WINDOW) ---
bedtools slop \
    -i "$OUTDIR/genes_of_interest.bed" \
    -g "$GENOME_SIZES" \
    -b $WINDOW \
    > "$OUTDIR/promoter_windows_${WINDOW}bp.bed"

echo "Promoter windows (TSS ± ${WINDOW}bp) saved to:"
echo "  $OUTDIR/promoter_windows_${WINDOW}bp.bed"

# --- EXTRACT FASTA SEQUENCES FOR MOTIF SCANNING ---
if [ -f "$GENOME_FASTA" ]; then
    bedtools getfasta \
        -fi "$GENOME_FASTA" \
        -bed "$OUTDIR/promoter_windows_${WINDOW}bp.bed" \
        -name \
        -fo "$OUTDIR/promoter_sequences.fa"
    echo "FASTA sequences extracted to: $OUTDIR/promoter_sequences.fa"
else
    echo "WARNING: Genome FASTA not found at $GENOME_FASTA"
    echo "Update GENOME_FASTA path and re-run, OR download with:"
    echo "  wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz"
    echo "  gunzip hg38.fa.gz && samtools faidx hg38.fa"
fi

echo ""
echo "DONE. Next step: run 02_scan_ap1_motifs.sh"
echo "(Remember to set CELL=\"${CELL}\" in 02_scan_ap1_motifs.sh too!)"
