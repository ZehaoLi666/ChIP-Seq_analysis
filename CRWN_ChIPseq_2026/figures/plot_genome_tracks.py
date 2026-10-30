#!/usr/bin/env python3
"""
IGV-style stacked coverage tracks across the concatenated P. falciparum genome.

Eight rows: five CPM tracks (3 IP replicates, IgG, Input) on one shared
y-scale so heights are comparable, then three IP-minus-control tracks.
Also checks, numerically, whether the strongest bin per chromosome coincides
with the AT-derived centromere -- the check that turned the "spikes look
periodic" impression into a measurement.

Usage:
  python plot_genome_tracks.py --bw-dir . --centromeres Pf3D7_v48_centromeres.csv
"""

import argparse
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pyBigWig

BIN_BP = 5000

CPM_TRACKS = [("CRWN_HA_D1.CPM.bw", "HA rep1", "#1f77b4"),
              ("CRWN_HA_D2.CPM.bw", "HA rep2", "#1f77b4"),
              ("CRWN_HA_D3.CPM.bw", "HA rep3", "#1f77b4"),
              ("CRWN_IGG.CPM.bw", "IgG", "#d95f02"),
              ("CRWN_INPUT.CPM.bw", "Input", "#7f7f7f")]
SUB_TRACKS = [("CRWN_HA_D1_minus_CRWN_INPUT.subtract.bw", "rep1 \u2212 Input", "#2c7fb8"),
              ("CRWN_HA_D2_minus_CRWN_INPUT.subtract.bw", "rep2 \u2212 Input", "#2c7fb8"),
              ("CRWN_HA_D3_minus_CRWN_INPUT.subtract.bw", "rep3 \u2212 Input", "#2c7fb8")]


def apply_style():
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"],
        "font.size": 7, "axes.titlesize": 8, "axes.labelsize": 7,
        "xtick.labelsize": 6, "ytick.labelsize": 6, "axes.linewidth": 0.6,
        "savefig.dpi": 300, "pdf.fonttype": 42, "ps.fonttype": 42,
    })


def nuclear_chroms(chrom_dict):
    nuc = {k: v for k, v in chrom_dict.items()
           if k.startswith("Pf3D7_") and "API" not in k}
    return dict(sorted(nuc.items(), key=lambda kv: int(kv[0].split("_")[1])))


def read_track(path, nuc, bin_bp=BIN_BP):
    bw = pyBigWig.open(path)
    out = {}
    for chrom, length in nuc.items():
        n = max(1, length // bin_bp)
        vals = bw.stats(chrom, 0, length, nBins=n, type="mean")
        out[chrom] = np.array([0.0 if v is None else float(v) for v in vals])
    bw.close()
    return out


def centromere_check(tracks, nuc, cen_mid, bin_bp=BIN_BP):
    """Distance from each chromosome's strongest bin to its centromere."""
    rows = []
    for name, prof in tracks.items():
        for chrom in nuc:
            y = prof[chrom]
            i = int(np.argmax(y))
            pos = i * bin_bp + bin_bp / 2
            rows.append({"track": name, "chrom": chrom,
                         "max_pos": int(pos), "max_val": round(float(y[i]), 3),
                         "dist_to_cen_kb": round((pos - cen_mid[chrom]) / 1000, 1)})
    return pd.DataFrame(rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bw-dir", default=".")
    ap.add_argument("--centromeres", required=True)
    ap.add_argument("--out", default="CRWN_genome_coverage_tracks")
    args = ap.parse_args()
    apply_style()

    first = os.path.join(args.bw_dir, CPM_TRACKS[0][0])
    bw = pyBigWig.open(first)
    nuc = nuclear_chroms(bw.chroms())
    bw.close()

    order = CPM_TRACKS + SUB_TRACKS
    tracks = {label: read_track(os.path.join(args.bw_dir, fn), nuc)
              for fn, label, _ in order}

    cen = pd.read_csv(args.centromeres)
    cen_mid = {r.chrom: (r.start + r.end) / 2.0
               for r in cen.itertuples() if r.chrom in nuc}

    # concatenated genome coordinates
    offs, cum = {}, 0
    for chrom, length in nuc.items():
        offs[chrom] = cum
        cum += length
    genome_mb = cum / 1e6
    bounds = [(offs[c] / 1e6, (offs[c] + nuc[c]) / 1e6) for c in nuc]
    cen_x = [(cen_mid[c] + offs[c]) / 1e6 for c in nuc]

    cpm_labels = [lab for _, lab, _ in CPM_TRACKS]
    sub_labels = [lab for _, lab, _ in SUB_TRACKS]
    all_cpm = np.concatenate([np.concatenate([tracks[l][c] for c in nuc])
                              for l in cpm_labels])
    all_sub = np.concatenate([np.concatenate([tracks[l][c] for c in nuc])
                              for l in sub_labels])
    y_cpm = float(np.percentile(all_cpm, 99.8))
    y_sub = float(np.percentile(np.abs(all_sub), 99.8))

    colour = {lab: col for _, lab, col in order}
    n = len(order)
    fig = plt.figure(figsize=(9.0, 8.1))
    top, height, gap = 0.930, 0.078, 0.0075
    axes = [fig.add_axes([0.115, top - (i + 1) * (height + gap), 0.865, height])
            for i in range(n)]

    for i, (_, label, _) in enumerate(order):
        ax = axes[i]
        is_sub = label in sub_labels
        for chrom in nuc:
            y = tracks[label][chrom]
            x = (offs[chrom] + np.arange(len(y)) * BIN_BP + BIN_BP / 2) / 1e6
            ax.fill_between(x, 0, y, color=colour[label], lw=0, step="mid")
        for j, (a, _b) in enumerate(bounds):
            if j:
                ax.axvline(a, color="#ffffff", lw=1.6, zorder=5)
                ax.axvline(a, color="#c8c8c8", lw=0.45, zorder=6)
        if is_sub:
            ax.axhline(0, color="#444444", lw=0.5, zorder=4)
            ax.set_ylim(-y_sub * 0.5, y_sub * 1.12)
            ax.set_yticks([0, 4])
        else:
            ax.set_ylim(0, y_cpm)
            ax.set_yticks([0, 13])
        ax.set_xlim(0, genome_mb)
        for cx in cen_x:
            ax.plot(cx, ax.get_ylim()[1], marker="v", ms=2.5, color="#c1272d",
                    clip_on=False, zorder=8)
        ax.tick_params(labelsize=6, length=2)
        ax.set_xticks([])
        for side in ("top", "right", "bottom"):
            ax.spines[side].set_visible(False)
        ax.spines["left"].set_linewidth(0.6)
        ax.text(-0.012, 0.5, label, transform=ax.transAxes, ha="right",
                va="center", fontsize=7)

    fig.text(0.985, top + 0.004, "CPM", ha="right", fontsize=6, color="#666666")
    fig.text(0.985, top - len(cpm_labels) * (height + gap) + 0.004,
             "IP \u2212 Input (CPM)", ha="right", fontsize=6, color="#666666")

    ax_last = axes[-1]
    ax_last.set_xticks([(a + b) / 2 for a, b in bounds])
    ax_last.set_xticklabels([str(i) for i in range(1, len(nuc) + 1)], fontsize=6.5)
    ax_last.tick_params(axis="x", length=0, pad=2)
    fig.text(0.5475, 0.212, "Chromosome", ha="center", fontsize=8)
    fig.text(0.115, 0.978,
             "CRWN-HA ChIP-seq coverage across the P. falciparum 3D7 genome",
             fontsize=9.5, ha="left")
    fig.text(0.115, 0.960,
             "5-kb bins. Red triangles mark centromeres.",
             fontsize=6.5, ha="left", color="#555555")

    fig.savefig(f"{args.out}.png", dpi=300)
    fig.savefig(f"{args.out}.pdf")

    check = centromere_check({l: tracks[l] for l in sub_labels}, nuc, cen_mid)
    check.to_csv("CRWN_centromere_peak_coincidence.csv", index=False)
    for label in sub_labels:
        d = check[check.track == label].dist_to_cen_kb.abs()
        print(f"{label:16s} within10kb={(d <= 10).sum()}/{len(d)} "
              f"median|d|={d.median():.1f}kb")


if __name__ == "__main__":
    main()
