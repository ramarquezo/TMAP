#!/bin/bash
# ============================================================
# MASTER RUN SCRIPT
# AP-1 ATAC-seq analysis: Starvation vs Conditioned Media
# ============================================================

set -euo pipefail

echo "=============================================="
echo " AP-1 Motif Analysis at Gene Promoters"
echo " GPRC5A | NFKB2 | RELB | ARNT"
echo "=============================================="

# --- STEP 0: Setup environment ---
echo "[0/4] Activating conda environment..."
#conda activate ap1_atac

# --- STEP 1: Get promoter sequences ---
echo "[1/4] Extracting promoter regions..."
bash 01_get_promoters.sh

# --- STEP 2: Scan AP-1 motifs ---
echo "[2/4] Scanning for AP-1 motifs (FIMO)..."
bash 02_scan_ap1_motifs.sh

# --- STEP 3: Quantify ATAC signal ---
echo "[3/4] Quantifying ATAC signal per condition..."
python3 03_compute_atac_signal.py

# --- STEP 4: Generate figure ---
echo "[4/4] Generating publication figure..."
python3 04_plot_tracks.py

echo ""
echo "=============================================="
echo " DONE! Check figures/ directory for outputs"
echo "=============================================="


# ============================================================
# EXPECTED nf-core ATAC-seq OUTPUT STRUCTURE
# Use this to find your files
# ============================================================
: '
results/
├── bigwig/
│   ├── Starvation_rep1.mLb.clN.bigWig      ← normalized bigwigs
│   ├── Starvation_rep2.mLb.clN.bigWig
│   ├── Starvation_rep3.mLb.clN.bigWig
│   ├── CM_rep1.mLb.clN.bigWig
│   ├── CM_rep2.mLb.clN.bigWig
│   └── CM_rep3.mLb.clN.bigWig
├── bwa/
│   └── mergedLibrary/
│       ├── Starvation_rep1.mLb.clN.sorted.bam   ← BAM files
│       └── ...
├── macs2/
│   ├── Starvation_rep1/
│   │   └── *_peaks.narrowPeak
│   └── consensus/
│       └── consensus_peaks.mLb.clN.bed     ← consensus peaks
└── deeptools/
    └── plotFingerprint/                    ← QC plots
'


# ============================================================
# TROUBLESHOOTING
# ============================================================
: '
Q: FIMO not found
A: conda install -c bioconda meme

Q: pyBigWig import error
A: pip install pyBigWig

Q: bedtools not found
A: conda install -c bioconda bedtools

Q: No AP-1 motifs found for a gene
A: Try widening the window in 01_get_promoters.sh (WINDOW=5000)
   Or lower the FIMO threshold in 02_scan_ap1_motifs.sh (--thresh 1e-3)

Q: BigWig signal is all zeros
A: Check your bigwig file paths in 04_plot_tracks.py
   Verify chrom naming (chr1 vs 1) matches your bigwig files:
   bigWigInfo your_file.bigWig | head

Q: TSS coordinates seem off
A: Verify with:
   grep "GPRC5A" results/promoters/genes_of_interest.bed
   Check against UCSC or Ensembl for your gene of interest
'
