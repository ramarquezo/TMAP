#!/usr/bin/env python3
"""
STEP 3: Quantify ATAC-seq signal at AP-1 motif sites
Per condition: Starvation vs Conditioned Media
Uses deepTools computeMatrix - only at motifs within open chromatin peaks
"""

import os
import subprocess
import pandas as pd
import numpy as np
from pathlib import Path

# ============================================================
# UPDATE THIS FOR EACH CELL LINE
# ============================================================
CELL="SAMPLE1"   # <-- # Change CELL="SAMPLE1" / "SAMPLE2" / "SAMPLE3" in each script

# ============================================================
# PATHS - auto-built from CELL variable
# ============================================================
BASE = "/mnt/Z/ATAC_PPG_plus_MF022UCD65_MF032UCD12/RAMO_Analysis"
CELL_BW_DIR = f"{BASE}/{CELL}_analysis_03-06-2025/{CELL}_results_nf-core/bwa/merged_library/bigwig"
ANALYSIS_DIR = f"{BASE}/Plot_AP-1_sites/{CELL}_AP-1"

# BigWig files
BIGWIG_FILES = {
    "Starvation": [
        f"{CELL_BW_DIR}/CONTROL_REP1.mLb.clN.bigWig",
        f"{CELL_BW_DIR}/CONTROL_REP2.mLb.clN.bigWig",
        f"{CELL_BW_DIR}/CONTROL_REP3.mLb.clN.bigWig",
    ],
    "CM": [
        f"{CELL_BW_DIR}/CM_REP1.mLb.clN.bigWig",
        f"{CELL_BW_DIR}/CM_REP2.mLb.clN.bigWig",
        f"{CELL_BW_DIR}/CM_REP3.mLb.clN.bigWig",
    ]
}

# Input files
MOTIF_BED        = f"{ANALYSIS_DIR}/results/motifs/AP1_motif_sites.bed"
CM_PEAKS         = f"{ANALYSIS_DIR}/results/{CELL}_CM_peaks.bed"
CONTROL_PEAKS    = f"{ANALYSIS_DIR}/results/{CELL}_CONTROL_peaks.bed"

# Filtered motif files (motifs within open chromatin only)
MOTIF_IN_CM      = f"{ANALYSIS_DIR}/results/motifs/AP1_motifs_in_CM_peaks.bed"
MOTIF_IN_CONTROL = f"{ANALYSIS_DIR}/results/motifs/AP1_motifs_in_CONTROL_peaks.bed"
MOTIF_IN_EITHER  = f"{ANALYSIS_DIR}/results/motifs/AP1_motifs_in_either_peaks.bed"

# Output directory
OUTDIR = Path(f"{ANALYSIS_DIR}/results/signal_quantification")
OUTDIR.mkdir(parents=True, exist_ok=True)

# Genes of interest and their chromosomes
GENES = {
    "GPRC5A": "chr12",
    "NFKB2":  "chr10",
    "RELB":   "chr19",
    "ARNT":   "chr1"
}

print("="*60)
print(f"STEP 3: Quantifying ATAC signal at AP-1 motif sites")
print(f"Cell line: {CELL}")
print("="*60)

# ============================================================
# STEP 3A: Filter motifs to open chromatin regions
# ============================================================

print("\n[1/3] Filtering AP-1 motifs to open chromatin regions...")

# Motifs in CM peaks
os.system(f"bedtools intersect -a {MOTIF_BED} -b {CM_PEAKS} -wa > {MOTIF_IN_CM}")
print(f"  Motifs in CM peaks:      {sum(1 for _ in open(MOTIF_IN_CM))}")

# Motifs in CONTROL peaks
os.system(f"bedtools intersect -a {MOTIF_BED} -b {CONTROL_PEAKS} -wa > {MOTIF_IN_CONTROL}")
print(f"  Motifs in CONTROL peaks: {sum(1 for _ in open(MOTIF_IN_CONTROL))}")

# Motifs in EITHER condition (union) - use for computeMatrix
os.system(f"cat {MOTIF_IN_CM} {MOTIF_IN_CONTROL} | sort -k1,1 -k2,2n | uniq > {MOTIF_IN_EITHER}")
print(f"  Motifs in either peaks:  {sum(1 for _ in open(MOTIF_IN_EITHER))}")

# Per gene breakdown
print("\n  Open chromatin AP-1 motifs per gene:")
print("  Gene      CM    CONTROL")
print("  --------  ----  -------")
for gene, chrom in GENES.items():
    cm_count = int(os.popen(f"grep -c '^{chrom}' {MOTIF_IN_CM} 2>/dev/null || echo 0").read().strip())
    ctrl_count = int(os.popen(f"grep -c '^{chrom}' {MOTIF_IN_CONTROL} 2>/dev/null || echo 0").read().strip())
    print(f"  {gene:<8}  {cm_count:<4}  {ctrl_count}")

# ============================================================
# STEP 3B: Compute matrix using ONLY open chromatin motifs
# ============================================================

print("\n[2/3] Computing signal matrix at open chromatin AP-1 motifs...")

def run_compute_matrix(gene, chrom, bigwig_files, motif_bed, outdir, window=500):
    """Run deepTools computeMatrix for open chromatin AP-1 motifs of one gene"""

    # Filter to this gene's chromosome
    gene_motif_bed = outdir / f"{gene}_AP1_motifs_open.bed"
    os.system(f"grep '^{chrom}\\b' {motif_bed} > {gene_motif_bed}")

    if not os.path.exists(gene_motif_bed) or os.path.getsize(gene_motif_bed) == 0:
        print(f"  WARNING: No open chromatin AP-1 motifs for {gene}")
        return None

    n_motifs = sum(1 for _ in open(gene_motif_bed))
    print(f"\n  {gene}: {n_motifs} motifs in open chromatin")

    all_bw = bigwig_files["Starvation"] + bigwig_files["CM"]
    bw_str = " ".join(all_bw)
    matrix_file = outdir / f"{gene}_matrix_open.gz"

    cmd = (f"computeMatrix reference-point "
           f"--referencePoint center "
           f"--regionsFileName {gene_motif_bed} "
           f"--scoreFileName {bw_str} "
           f"--outFileName {matrix_file} "
           f"--beforeRegionStartLength {window} "
           f"--afterRegionStartLength {window} "
           f"--missingDataAsZero "
           f"--numberOfProcessors 4")

    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)

    if result.returncode != 0:
        print(f"  ERROR: {result.stderr[:300]}")
        return None

    return matrix_file


def extract_mean_signals(matrix_file, n_starv, n_cm):
    """Extract mean ATAC signal from deepTools matrix"""
    import gzip, json

    with gzip.open(matrix_file, 'rb') as f:
        header = json.loads(f.readline().decode().strip('@').strip())
        data = np.loadtxt(f, comments='@', usecols=range(6, 606))

    if data.ndim == 1:
        data = data.reshape(1, -1)

    n_bins = data.shape[1] // (n_starv + n_cm)

    signals = {}
    for condition, indices in [("Starvation", range(n_starv)),
                                ("CM", range(n_starv, n_starv + n_cm))]:
        sample_data = []
        for s in indices:
            col_start = s * n_bins
            col_end = col_start + n_bins
            sample_data.append(np.nanmean(data[:, col_start:col_end]))
        signals[condition] = np.mean(sample_data)

    return signals


# ============================================================
# STEP 3C: Run for all genes
# ============================================================

print("\n[3/3] Extracting signal values per gene per condition...")

results = []

for gene, chrom in GENES.items():

    matrix_file = run_compute_matrix(
        gene, chrom, BIGWIG_FILES, MOTIF_IN_EITHER, OUTDIR
    )

    if matrix_file and os.path.exists(matrix_file):
        try:
            signals = extract_mean_signals(
                matrix_file,
                n_starv=len(BIGWIG_FILES["Starvation"]),
                n_cm=len(BIGWIG_FILES["CM"])
            )
            log2fc = np.log2((signals["CM"] + 1e-6) / (signals["Starvation"] + 1e-6))
            results.append({
                "Gene":                    gene,
                "Cell_line":               CELL,
                "Starvation_mean":         round(signals["Starvation"], 4),
                "CM_mean":                 round(signals["CM"], 4),
                "log2FC_CM_vs_Starvation": round(log2fc, 4)
            })
            print(f"  {gene}: Starvation={signals['Starvation']:.4f}  "
                  f"CM={signals['CM']:.4f}  log2FC={log2fc:.4f}")
        except Exception as e:
            print(f"  {gene}: Could not extract signals — {e}")
    else:
        print(f"  {gene}: Skipped — no open chromatin motifs found")

# ============================================================
# SAVE RESULTS
# ============================================================

if results:
    df = pd.DataFrame(results)
    out_file = OUTDIR / f"AP1_signal_summary_{CELL}.tsv"
    df.to_csv(out_file, sep="\t", index=False)
    print(f"\nSignal summary saved to: {out_file}")
    print("\n" + "="*60)
    print("SUMMARY:")
    print("="*60)
    print(df.to_string(index=False))
else:
    print("\nWARNING: No results generated!")

print(f"\nDONE. Next step: run 04_plot_tracks.py")
print(f"(Remember to set CELL = '{CELL}' in 04_plot_tracks.py too!)")
