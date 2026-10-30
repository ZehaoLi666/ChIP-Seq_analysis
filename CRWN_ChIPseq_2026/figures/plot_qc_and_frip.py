#!/usr/bin/env python3
"""
QC figures: replicate/control correlation structure, and FRiP + peak reproducibility.

Figure 1 (correlation) applies the decision rule that controls must NOT
correlate with IP samples more strongly than true replicates do. The per-
replicate verdict is computed, not eyeballed: for each IP, is its WORST
sibling correlation still above its BEST control correlation?

Figure 2 (FRiP) plots fraction of reads in peaks per replicate x control x
peak mode, plus the 3-of-3 vs union peak counts.

Inputs are the deepTools correlation matrices and the FRiP/reproducibility
tables emitted by 02_CRWN_qc_bigwig.sh and 03_CRWN_peaks.sh.

Usage:
  python plot_qc_and_frip.py --qc-dir deeptools --frip stats/CRWN_FRiP.tsv
"""

import argparse

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

ORDER = ["HA_D1", "HA_D2", "HA_D3", "IgG", "Input"]
SHORT = {"HA_D1": "rep1", "HA_D2": "rep2", "HA_D3": "rep3",
         "IgG": "IgG", "Input": "Input"}
IPS = ["HA_D1", "HA_D2", "HA_D3"]
CONTROLS = ["IgG", "Input"]

# 3-of-3 reproducible vs union peak counts, from 03_CRWN_peaks.sh
REPRO = pd.DataFrame({
    "control": ["IGG", "IGG", "INPUT", "INPUT"],
    "mode": ["narrow", "broad", "narrow", "broad"],
    "reproducible_3of3": [39, 111, 167, 255],
    "union": [1401, 2430, 1489, 2342],
})


def apply_style():
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"],
        "font.size": 7, "axes.titlesize": 8, "axes.labelsize": 7,
        "xtick.labelsize": 6, "ytick.labelsize": 6, "axes.linewidth": 0.6,
        "savefig.dpi": 300, "pdf.fonttype": 42, "ps.fonttype": 42,
    })


def open_frame(ax):
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    ax.spines["left"].set_position(("outward", 3))
    ax.spines["bottom"].set_position(("outward", 3))


def read_matrix(path):
    """deepTools writes a comment header line and quotes the labels."""
    df = pd.read_csv(path, sep="\t", skiprows=1, index_col=0)
    df.index = [i.strip("'") for i in df.index]
    df.columns = [c.strip("'") for c in df.columns]
    return df.loc[ORDER, ORDER]


def diagnose(matrix, label):
    """Per-replicate: closer to its siblings, or to a control?"""
    out = []
    for s in IPS:
        sib = [matrix.loc[s, o] for o in IPS if o != s]
        ctl = [matrix.loc[s, o] for o in CONTROLS]
        out.append({"matrix": label, "replicate": s,
                    "sibling_min": round(min(sib), 4),
                    "control_max": round(max(ctl), 4),
                    "passes_rule": bool(min(sib) > max(ctl))})
    return pd.DataFrame(out)


def plot_correlation(bam_p, bw_p, out_stem="CRWN_chipseq_QC_correlation"):
    labels = [SHORT[o] for o in ORDER]
    sib_min = [min(bw_p.loc[s, o] for o in IPS if o != s) for s in IPS]
    ctl_max = [max(bw_p.loc[s, o] for o in CONTROLS) for s in IPS]

    fig = plt.figure(figsize=(8.6, 3.2))
    ax_a = fig.add_axes([0.055, 0.24, 0.245, 0.60])
    ax_b = fig.add_axes([0.345, 0.24, 0.245, 0.60])
    cax = fig.add_axes([0.605, 0.24, 0.011, 0.60])
    ax_c = fig.add_axes([0.735, 0.24, 0.235, 0.60])

    for ax, matrix, title, show_y in ((ax_a, bam_p, "10-kb read counts", True),
                                      (ax_b, bw_p, "CPM-normalized tracks", False)):
        v = matrix.values.astype(float)
        im = ax.imshow(v, cmap="viridis", vmin=0.80, vmax=1.0)
        for i in range(len(ORDER)):
            for j in range(len(ORDER)):
                ax.text(j, i, f"{v[i, j]:.2f}", ha="center", va="center",
                        fontsize=6, color="white" if v[i, j] < 0.93 else "black")
        ax.set_xticks(range(len(ORDER))); ax.set_yticks(range(len(ORDER)))
        ax.set_xticklabels(labels, rotation=45, ha="right")
        ax.set_yticklabels(labels if show_y else [])
        ax.set_title(title, loc="left")
        for spine in ax.spines.values():
            spine.set_visible(False)
        ax.tick_params(length=0)

    cb = fig.colorbar(im, cax=cax, ticks=[0.80, 0.90, 1.0])
    cb.set_label("Pearson r", fontsize=6)
    cb.ax.tick_params(labelsize=6)

    ys = np.arange(len(IPS))
    ax_c.hlines(y=ys, xmin=np.minimum(sib_min, ctl_max),
                xmax=np.maximum(sib_min, ctl_max), color="#cccccc", lw=1.2)
    ax_c.scatter(sib_min, ys, s=40, color="#1f77b4", zorder=3,
                 label="worst sibling replicate")
    ax_c.scatter(ctl_max, ys, s=40, color="#d95f02", marker="D", zorder=3,
                 label="best-matching control")
    ax_c.set_yticks(ys); ax_c.set_yticklabels([SHORT[s] for s in IPS])
    ax_c.set_xlabel("Pearson r (CPM tracks)")
    ax_c.set_title("Controls rival replicates", loc="left")
    ax_c.set_xlim(0.87, 1.015); ax_c.margins(y=0.34)
    ax_c.legend(frameon=False, fontsize=6, loc="lower left",
                handletextpad=0.4, borderpad=0.2, labelspacing=0.55)
    open_frame(ax_c)

    for ax, letter in ((ax_a, "a"), (ax_b, "b"), (ax_c, "c")):
        ax.text(-0.17, 1.14, letter, transform=ax.transAxes,
                fontsize=9, fontweight="bold", va="top")

    fig.savefig(f"{out_stem}.png", dpi=300)
    fig.savefig(f"{out_stem}.pdf")


def plot_frip(frip, out_stem="CRWN_chipseq_peaks_FRiP"):
    reps = ["CRWN_HA_D1", "CRWN_HA_D2", "CRWN_HA_D3"]
    combos = [("IGG", "narrow"), ("IGG", "broad"),
              ("INPUT", "narrow"), ("INPUT", "broad")]
    cols = {"IGG": "#d95f02", "INPUT": "#1f77b4"}

    fig = plt.figure(figsize=(7.8, 3.1))
    ax_a = fig.add_axes([0.075, 0.20, 0.35, 0.62])
    ax_b = fig.add_axes([0.585, 0.20, 0.34, 0.62])

    xs = np.arange(len(reps))
    width = 0.19
    for k, (ctl, mode) in enumerate(combos):
        sub = frip[(frip.control == ctl) & (frip["mode"] == mode)]
        sub = sub.set_index("IP").loc[reps]
        ax_a.bar(xs + (k - 1.5) * width, sub.FRiP_pct, width=width,
                 color=cols[ctl], alpha=1.0 if mode == "broad" else 0.55,
                 edgecolor=cols[ctl], linewidth=0.6,
                 hatch="" if mode == "broad" else "///")
    ax_a.axhline(1.0, color="#888888", lw=0.8, ls=":", zorder=0)
    ax_a.annotate("1% = weak-IP floor", xy=(2.42, 1.0), xytext=(2.42, 3.4),
                  fontsize=6, color="#555555", ha="right",
                  arrowprops=dict(arrowstyle="-", color="#888888", lw=0.5))
    ax_a.set_xticks(xs); ax_a.set_xticklabels(["rep1", "rep2", "rep3"])
    ax_a.set_ylabel("FRiP (% of reads in peaks)")
    ax_a.set_title("Enrichment falls sharply rep1 \u2192 rep3", loc="left")
    ax_a.set_ylim(0, 21.5)
    ax_a.legend(handles=[
        Patch(facecolor=cols["IGG"], alpha=0.55, hatch="///",
              edgecolor=cols["IGG"], label="vs IgG, narrow"),
        Patch(facecolor=cols["IGG"], edgecolor=cols["IGG"], label="vs IgG, broad"),
        Patch(facecolor=cols["INPUT"], alpha=0.55, hatch="///",
              edgecolor=cols["INPUT"], label="vs Input, narrow"),
        Patch(facecolor=cols["INPUT"], edgecolor=cols["INPUT"],
              label="vs Input, broad")],
        frameon=False, fontsize=6, loc="upper right", handlelength=1.1,
        labelspacing=0.3)
    open_frame(ax_a)

    repro = REPRO.copy()
    repro["pct_of_union"] = (100 * repro.reproducible_3of3 / repro.union).round(1)
    ys = np.arange(len(repro))
    labels = [f"{c.lower()}, {m}" for c, m in zip(repro.control, repro["mode"])]
    ax_b.barh(ys, repro.union, color="#e4e4e4", edgecolor="#bbbbbb", height=0.6)
    ax_b.barh(ys, repro.reproducible_3of3, color="#2c7fb8", height=0.6)
    for i, (u, r3, pct) in enumerate(zip(repro.union, repro.reproducible_3of3,
                                         repro.pct_of_union)):
        ax_b.text(u + 45, i, f"{r3}/{u}  ({pct}%)", va="center", fontsize=6)
    ax_b.set_yticks(ys); ax_b.set_yticklabels(labels)
    ax_b.set_xlabel("Peaks")
    ax_b.set_title("Few peaks reproduce across replicates", loc="left")
    ax_b.set_xlim(0, 3900); ax_b.set_xticks([0, 1000, 2000, 3000])
    ax_b.invert_yaxis()
    ax_b.legend(handles=[Patch(facecolor="#2c7fb8", label="in all 3 replicates"),
                         Patch(facecolor="#e4e4e4", edgecolor="#bbbbbb",
                               label="union of 3 replicates")],
                frameon=False, fontsize=6, loc="upper right", handlelength=1.1,
                labelspacing=0.4, bbox_to_anchor=(1.02, 0.62))
    open_frame(ax_b)

    for ax, letter in ((ax_a, "a"), (ax_b, "b")):
        ax.text(-0.17, 1.12, letter, transform=ax.transAxes,
                fontsize=9, fontweight="bold", va="top")

    fig.savefig(f"{out_stem}.png", dpi=300)
    fig.savefig(f"{out_stem}.pdf")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--qc-dir", default="deeptools")
    ap.add_argument("--frip", default="stats/CRWN_FRiP.tsv")
    args = ap.parse_args()
    apply_style()

    bam_p = read_matrix(f"{args.qc_dir}/mapping_pearson_correlation_matrix.tsv")
    bam_s = read_matrix(f"{args.qc_dir}/mapping_spearman_correlation_matrix.tsv")
    bw_p = read_matrix(f"{args.qc_dir}/bigwig_CPM_pearson_correlation_matrix.tsv")

    verdict = pd.concat([diagnose(bam_p, "BAM Pearson"),
                         diagnose(bam_s, "BAM Spearman"),
                         diagnose(bw_p, "BigWig Pearson")], ignore_index=True)
    verdict.to_csv("CRWN_replicate_vs_control_rule.csv", index=False)
    print(verdict.to_string(index=False))

    for matrix, name in ((bam_p, "bam_pearson"), (bam_s, "bam_spearman"),
                         (bw_p, "bigwig_CPM_pearson")):
        matrix.round(4).to_csv(f"CRWN_correlation_{name}.csv")

    plot_correlation(bam_p, bw_p)

    frip = pd.read_csv(args.frip, sep="\t")
    if "FRiP_pct" not in frip.columns:
        frip["FRiP_pct"] = (100 * frip["FRiP"]).round(2)
    frip.to_csv("CRWN_FRiP_peaks_summary.csv", index=False)
    plot_frip(frip)
    print("\nmean FRiP by replicate:")
    print(frip.groupby("IP").FRiP_pct.mean().round(2).to_string())


if __name__ == "__main__":
    main()
