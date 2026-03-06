#!/usr/bin/env python3
"""
STEP 4: Publication-quality figure
Layout per gene (4 columns):
  Row 0: Overlaid ATAC tracks (CM + Starvation, transparent fills)
  Row 1: AP-1 motif accessibility (CM=red, Starv=blue, both=grey)
  Row 2: Gene annotation + TSS arrow
  Row 3: Barplot ATAC log2FC + RNA-seq log2FC side by side
  Row 4: Heatmap spanning all 4 genes (CM - Starvation)
"""

import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.gridspec as gridspec
from matplotlib.colors import TwoSlopeNorm
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

# ============================================================
# UPDATE THIS FOR EACH CELL LINE
# ============================================================
CELL="SAMPLE1"   # <-- # Change CELL="SAMPLE1" / "SAMPLE2" / "SAMPLE3" in each script

# ============================================================
# RNA-seq log2FC values (CM vs Starvation)
# Use None for cell lines without RNA-seq
# ============================================================
RNA_LOG2FC = {
    "UCD04": {"GeneA": 1.20,   "GeneB": 0.53,   "GeneC": 1.34,   "GeneD": -0.33},  # <-- UPDATE
    "UCD12": {"GeneA": 1.20,   "GeneB": 0.53,   "GeneC": 1.14,   "GeneD": -0.24},  # <-- UPDATE
    "UCD65": {"GeneA": 1.19,   "GeneB": 0.52,   "GeneC": 1.54,   "GeneD": -0.42},  # <-- UPDATE
}

# ============================================================
# PATHS
# ============================================================
BASE         = "/mnt/Z/ATAC_PPG_plus_MF022UCD65_MF032UCD12/RAMO_Analysis"
ANALYSIS_DIR = f"{BASE}/Plot_AP-1_sites/{CELL}_AP-1"
CELL_BW_DIR  = f"{BASE}/{CELL}_analysis_03-06-2025/{CELL}_results_nf-core/bwa/merged_library/bigwig"

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

SIGNAL_SUMMARY = f"{ANALYSIS_DIR}/results/signal_quantification/AP1_signal_summary_{CELL}.tsv"
MOTIF_BED      = f"{ANALYSIS_DIR}/results/motifs/AP1_motif_sites.bed"
MOTIF_IN_CM    = f"{ANALYSIS_DIR}/results/motifs/AP1_motifs_in_CM_peaks.bed"
MOTIF_IN_CTRL  = f"{ANALYSIS_DIR}/results/motifs/AP1_motifs_in_CONTROL_peaks.bed"

OUTDIR = Path(f"{ANALYSIS_DIR}/figures")
OUTDIR.mkdir(exist_ok=True)

# Gene metadata
GENES = {
    "GeneA": {"chr": "chr12", "tss": 12891562,  "strand": "+", "window": 3000}, # <- Update
    "GeneB":  {"chr": "chr10", "tss": 102395705, "strand": "+", "window": 3000}, # <- Update
    "GeneC":   {"chr": "chr19", "tss": 45001464,  "strand": "+", "window": 3000}, # <- Update
    "GeneD":   {"chr": "chr1",  "tss": 150876599, "strand": "-", "window": 3000}, # <- Update
}

COLORS = {
    "Starvation": "#4878CF",
    "CM":         "#D65F5F",
    "MOTIF_CM":   "#D65F5F",
    "MOTIF_CTRL": "#4878CF",
    "MOTIF_BOTH": "#95a5a6",
    "MOTIF_NONE": "#ecf0f1",
    "UP":         "#2ecc71",
    "DOWN":       "#e74c3c",
    "ATAC_bar":   "#3498db",
    "RNA_bar":    "#e67e22",
}

print("="*60)
print(f"STEP 4: Generating figure for {CELL}")
print("="*60)

# ============================================================
# HELPERS
# ============================================================
def get_bigwig_signal(bw_file, chrom, start, end, bins=300):
    try:
        import pyBigWig
        bw   = pyBigWig.open(bw_file)
        vals = bw.stats(chrom, start, end, type="mean", nBins=bins)
        bw.close()
        return np.array([v if v is not None else 0 for v in vals])
    except:
        return np.zeros(bins)

def get_mean_signal(bw_files, chrom, start, end, bins=300):
    sigs = [get_bigwig_signal(f, chrom, start, end, bins)
            for f in bw_files if os.path.exists(f)]
    return np.mean(sigs, axis=0) if sigs else np.zeros(bins)

def get_std_signal(bw_files, chrom, start, end, bins=300):
    sigs = [get_bigwig_signal(f, chrom, start, end, bins)
            for f in bw_files if os.path.exists(f)]
    return np.std(sigs, axis=0) if len(sigs) > 1 else np.zeros(bins)

def load_motifs(bed_file, chrom, start, end):
    try:
        df   = pd.read_csv(bed_file, sep="\t", header=None,
                           names=["chr","start","end","name","score","strand"])
        mask = (df["chr"]==chrom) & (df["start"]>=start) & (df["end"]<=end)
        return df[mask].to_dict("records")
    except:
        return []

def motif_center(m, tss):
    return ((m["start"] + m["end"]) / 2) - tss

def motif_in_list(mc, mlist, tss, tol=10):
    return any(abs(motif_center(m, tss) - mc) < tol for m in mlist)

# ============================================================
# LOAD DATA
# ============================================================
print("\nLoading ATAC signal summary...")
df_signal = pd.read_csv(SIGNAL_SUMMARY, sep="\t")
print(df_signal.to_string(index=False))

genes_list = list(GENES.keys())
n_genes    = len(genes_list)
BINS       = 300

# Pre-load all signals
print("\nLoading BigWig signals...")
data = {}
for gene, info in GENES.items():
    chrom  = info["chr"]
    tss    = info["tss"]
    start  = tss - info["window"]
    end    = tss + info["window"]
    data[gene] = {
        "ctrl_mean":   get_mean_signal(BIGWIG_FILES["Starvation"], chrom, start, end),
        "ctrl_std":    get_std_signal( BIGWIG_FILES["Starvation"], chrom, start, end),
        "cm_mean":     get_mean_signal(BIGWIG_FILES["CM"],         chrom, start, end),
        "cm_std":      get_std_signal( BIGWIG_FILES["CM"],         chrom, start, end),
        "all_motifs":  load_motifs(MOTIF_BED,    chrom, start, end),
        "cm_motifs":   load_motifs(MOTIF_IN_CM,  chrom, start, end),
        "ctrl_motifs": load_motifs(MOTIF_IN_CTRL,chrom, start, end),
    }
    print(f"  {gene}: done")

# ============================================================
# FIGURE LAYOUT
# 5 rows x 4 cols
# Row 0: overlaid tracks          height=4
# Row 1: motif accessibility      height=1
# Row 2: gene annotation          height=1
# Row 3: barplots                 height=3
# Row 4: heatmap (spans all cols) height=2
# ============================================================
fig = plt.figure(figsize=(18, 15))
fig.patch.set_facecolor('white')

gs = gridspec.GridSpec(
    5, n_genes,
    height_ratios=[4, 1, 1, 3, 2],
    hspace=0.20,
    wspace=0.30,
    left=0.07, right=0.95,
    top=0.93,  bottom=0.08
)

heatmap_matrix = np.zeros((n_genes, BINS))

for col, gene in enumerate(genes_list):
    info   = GENES[gene]
    tss    = info["tss"]
    window = info["window"]
    x_rel  = np.linspace(-window, window, BINS)

    d          = data[gene]
    ctrl_mean  = d["ctrl_mean"]
    ctrl_std   = d["ctrl_std"]
    cm_mean    = d["cm_mean"]
    cm_std     = d["cm_std"]
    all_motifs = d["all_motifs"]
    cm_motifs  = d["cm_motifs"]
    ctrl_motifs= d["ctrl_motifs"]

    ymax = max(
    (ctrl_mean + ctrl_std).max(),
    (cm_mean + cm_std).max()
    ) * 1.15 or 1

    # Store for heatmap
    heatmap_matrix[col] = cm_mean - ctrl_mean

    # ---- ROW 0: Overlaid ATAC tracks ----
    ax0 = fig.add_subplot(gs[0, col])

    # Starvation: blue transparent fill + line
    ax0.fill_between(x_rel, ctrl_mean,
                     alpha=0.40, color=COLORS["Starvation"], linewidth=0)
    ax0.fill_between(x_rel,
                     np.maximum(ctrl_mean - ctrl_std, 0),
                     ctrl_mean + ctrl_std,
                     alpha=0.06, color=COLORS["Starvation"])
    ax0.plot(x_rel, ctrl_mean,
             color=COLORS["Starvation"], linewidth=1.0, alpha=0.95,
             label="Starvation")

    # CM: red transparent fill + line on top
    ax0.fill_between(x_rel, cm_mean,
                     alpha=0.40, color=COLORS["CM"], linewidth=0)
    ax0.fill_between(x_rel,
                     np.maximum(cm_mean - cm_std, 0),
                     cm_mean + cm_std,
                     alpha=0.06, color=COLORS["CM"])
    ax0.plot(x_rel, cm_mean,
             color=COLORS["CM"], linewidth=1.0, alpha=0.95,
             label="Cond. Media")

    ax0.axvline(0, color='black', lw=0.8, ls='--', alpha=0.4)
    ax0.set_xlim(-window, window)
    ax0.set_ylim(0, ymax)
    ax0.spines[['top','right','bottom']].set_visible(False)
    ax0.set_xticks([-3000, -2000, -1000, 0, 1000, 2000, 3000])
    ax0.set_xticklabels(["-3kb", "", "", "TSS", "", "", "+3kb"], fontsize=6)
    ax0.spines['bottom'].set_visible(True)
    ax0.tick_params(labelsize=6)
    ax0.set_title(gene, fontsize=12, fontweight='bold', pad=4)

    # Legend only on first plot
    if col == 0:
        ax0.set_ylabel("ATAC Signal\n(mean ± SD)", fontsize=8)
        ax0.legend(fontsize=7, loc='upper right',
                   frameon=False, handlelength=1)
    else:
        ax0.set_yticklabels([])

    # ---- ROW 1: Motif accessibility ----
    ax1 = fig.add_subplot(gs[1, col])
    ax1.set_xlim(-window, window)
    ax1.set_ylim(0, 3)   # makes lines appear shorter (only bottom third visible)
    ax1.tick_params(bottom=False, labelbottom=False)  # removes small ticks
    ax1.set_xticklabels([])
    ax1.set_yticks([])
    ax1.spines[['top','right','left','bottom']].set_visible(False)
    ax1.axvline(0, color='black', lw=0.5, ls='--', alpha=0.3)

    for m in all_motifs:
        mc      = motif_center(m, tss)
        in_cm   = motif_in_list(mc, cm_motifs,   tss)
        in_ctrl = motif_in_list(mc, ctrl_motifs,  tss)
        if in_cm and in_ctrl:
            color, lw = COLORS["MOTIF_BOTH"], 2.0
        elif in_cm:
            color, lw = COLORS["MOTIF_CM"],   2.5
        elif in_ctrl:
            color, lw = COLORS["MOTIF_CTRL"], 2.5
        else:
            color, lw = COLORS["MOTIF_NONE"], 1.0
        ax1.axvline(mc, color=color, lw=lw, alpha=0.9)

    if col == 0:
        ax1.set_ylabel("AP-1\nmotifs", fontsize=7)
    ax1.text(window * 0.6, 0.65,
             f"n={len(all_motifs)}", fontsize=6,
             color='#8e44ad', style='italic')

    # ---- ROW 2: Gene annotation ----
    ax2 = fig.add_subplot(gs[2, col])
    ax2.set_xlim(-window, window)
    ax2.set_ylim(-0.5, 1.5)
    ax2.set_yticks([])
    ax2.spines[['top','right','left','bottom']].set_visible(False)
    ax2.plot([-window*0.6, window*0.6], [0.5, 0.5], color='#2c3e50', lw=2)
    arrow_dir = 1 if info["strand"] == "+" else -1
    ax2.annotate('', xy=(arrow_dir*250, 0.5), xytext=(0, 0.5),
                 arrowprops=dict(arrowstyle='->', color='#2c3e50', lw=1.5))
#    ax2.text(0, 1.2, gene, ha='center', fontsize=8,
#             fontweight='bold', color='#2c3e50')
#    ax2.text(0, -0.3, "TSS", ha='center', fontsize=6, color='gray')
    ax2.set_xticks([])
    ax2.set_xticklabels([])

    # ---- ROW 3: Barplot ATAC + RNA log2FC ----
    ax3 = fig.add_subplot(gs[3, col])

    row     = df_signal[df_signal["Gene"] == gene]
    atac_fc = float(row["log2FC_CM_vs_Starvation"]) if not row.empty else 0
    rna_fc  = RNA_LOG2FC[CELL].get(gene, None)

    if rna_fc is not None:
        bar_x    = [0.25, 0.75]
        bar_vals = [atac_fc, rna_fc]
        bar_cols = [COLORS["ATAC_bar"],
                    COLORS["UP"] if rna_fc >= 0 else COLORS["DOWN"]]
        bar_lbls = ["ATAC", "RNA"]
    else:
        bar_x    = [0.5]
        bar_vals = [atac_fc]
        bar_cols = [COLORS["ATAC_bar"]]
        bar_lbls = ["ATAC"]

    all_vals = [v for v in bar_vals if v is not None]
    ymin_bar = min(0, min(all_vals)) * 1.35
    ymax_bar = max(0, max(all_vals)) * 1.35
    if ymin_bar == 0: ymin_bar = -0.1
    if ymax_bar == 0: ymax_bar =  0.1

    bars = ax3.bar(bar_x, bar_vals, width=0.35,
                   color=bar_cols, edgecolor='black', linewidth=0.7)
    ax3.axhline(0, color='black', lw=0.8)
    ax3.set_xlim(0, 1)
    ax3.set_ylim(ymin_bar, ymax_bar)
    ax3.set_xticks(bar_x)
    ax3.set_xticklabels(bar_lbls, fontsize=9)
    ax3.tick_params(labelsize=7)
    ax3.spines[['top','right']].set_visible(False)

    if col == 0:
        ax3.set_ylabel("log₂FC\n(CM / Starvation)", fontsize=8)
    else:
        ax3.set_yticklabels([])

    # Value labels on bars
    for bar, val in zip(bars, bar_vals):
        va     = 'bottom' if val >= 0 else 'top'
        offset = (ymax_bar - ymin_bar) * 0.03
        offset = offset if val >= 0 else -offset
        ax3.text(bar.get_x() + bar.get_width()/2,
                 val + offset, f"{val:.2f}",
                 ha='center', va=va, fontsize=7, fontweight='bold')

# ---- ROW 4: Heatmap spanning all columns ----
ax4 = fig.add_subplot(gs[4, :])

norm = TwoSlopeNorm(
    vmin=min(heatmap_matrix.min(), -0.1),
    vcenter=0,
    vmax=max(heatmap_matrix.max(),  0.1)
)
im = ax4.imshow(heatmap_matrix, aspect='auto',
                cmap='RdBu_r', norm=norm,
                interpolation='bilinear')

# Add AP-1 motif tick marks above each gene row
for row_idx, gene in enumerate(genes_list):
    info    = GENES[gene]
    tss     = info["tss"]
    window  = info["window"]
    all_motifs = data[gene]["all_motifs"]

    for m in all_motifs:
        mc = motif_center(m, tss)  # position relative to TSS in bp

        # Convert bp position to pixel/bin coordinate
        x_bin = (mc + window) / (2 * window) * (BINS - 1)

        # Draw small tick mark just above the gene row
        ax4.plot(
            [x_bin, x_bin],
            [row_idx - 0.5, row_idx - 0.42],  # just above the row
            color='black',
            lw=1.0,
            clip_on=False
        )

# Add white separator lines between gene rows
for row_idx in range(1, n_genes):
    ax4.axhline(row_idx - 0.5, color='white', linewidth=2.0, zorder=5)

ax4.set_yticks(range(n_genes))
ax4.set_yticklabels(genes_list, fontsize=9, fontweight='bold')
ax4.set_xticks(np.linspace(0, BINS-1, 7))
ax4.set_xticklabels(["-3kb","-2kb","-1kb","TSS","+1kb","+2kb","+3kb"],
                    fontsize=8)
ax4.set_xlabel("Position relative to TSS (bp)", fontsize=9)
ax4.set_title("ATAC-seq Δ signal (CM − Starvation) at AP-1 motif regions",
              fontsize=9, pad=4)

cbar = plt.colorbar(im, ax=ax4, orientation='vertical',
                    fraction=0.015, pad=0.01)
cbar.set_label("Δ Signal\n(CM−Starv)", fontsize=7)
cbar.ax.tick_params(labelsize=6)

# ============================================================
# LEGEND
# ============================================================
legend_elements = [
    mpatches.Patch(color=COLORS["Starvation"], alpha=0.7, label="Starvation"),
    mpatches.Patch(color=COLORS["CM"],         alpha=0.7, label="Cond. Media"),
    mpatches.Patch(color=COLORS["MOTIF_CM"],   label="AP-1 motif: CM-accessible"),
    mpatches.Patch(color=COLORS["MOTIF_CTRL"], label="AP-1 motif: Starv-accessible"),
    mpatches.Patch(color=COLORS["MOTIF_BOTH"], label="AP-1 motif: both conditions"),
    mpatches.Patch(color=COLORS["ATAC_bar"],   label="ATAC log₂FC"),
    mpatches.Patch(color=COLORS["RNA_bar"],    label="RNA-seq log₂FC"),
]
fig.legend(handles=legend_elements, loc='lower center',
           ncol=7, fontsize=7, frameon=True,
           bbox_to_anchor=(0.5, 0.001))

fig.suptitle(
    f"AP-1 Chromatin Accessibility at Target Gene Promoters — {CELL}\n"
    "ATAC-seq: Starvation vs Conditioned Media",
    fontsize=12, fontweight='bold', y=0.985
)

# ============================================================
# SAVE
# ============================================================
out_pdf = OUTDIR / f"AP1_figure_{CELL}_v2.pdf"
out_png = OUTDIR / f"AP1_figure_{CELL}_v2.png"

plt.savefig(out_pdf, dpi=300, bbox_inches='tight', facecolor='white')
plt.savefig(out_png, dpi=150, bbox_inches='tight', facecolor='white')

print(f"\nFigure saved:")
print(f"  PDF: {out_pdf}")
print(f"  PNG: {out_png}")
print("\nDONE!")
