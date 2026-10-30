#!/usr/bin/env python3
"""
Generate CRWN_Troph_coverage_per_chromosome.png / .pdf

Per-chromosome IGV-style coverage tracks for the CRWN-HA trophozoite (16 hpi)
ChIP-seq, drawn as 14 stacked rows (one per nuclear chromosome) on a shared
absolute-Mb x-axis and a shared y-axis.

Signal = mean of the two Rep1 IP-minus-IgG subtraction tracks (C1, C2),
with negative values set to zero AFTER averaging, then summarized in 5-kb bins.
Averaging before clamping matters: a bin where one replicate is +2 and the
other -2 averages to 0 rather than to +1.

Inputs
------
  --bw1 / --bw2   the two CPM-normalized IP-minus-IgG BigWig tracks
  --centromeres   Pf3D7_v48_centromeres.csv (columns: chrom,start,end,...)
                  Centromere positions are INFERRED from local AT content, not
                  taken from the PlasmoDB v48 annotation (which has no
                  centromere features). The chromosome 10 position is
                  unresolved by that method and is drawn as an open triangle.

Usage
-----
  python plot_CRWN_Troph_coverage.py \
      --bw1 CRWN_Troph_Rep1_C1_minus_IgG_16hpi_Rep1.subtract.CPM.bw \
      --bw2 CRWN_Troph_Rep1_C2_minus_IgG_16hpi_Rep1.subtract.CPM.bw \
      --centromeres Pf3D7_v48_centromeres.csv \
      --out CRWN_Troph_coverage_per_chromosome

Requires: pyBigWig, numpy, pandas, matplotlib
"""

import argparse

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
import pyBigWig

# ---------------------------------------------------------------- parameters
BIN_BP = 5000          # display bin size
CAP_PERCENTILE = 99.9  # y-axis cap; bins above are marked, not flattened
END_WINDOW_BP = 100000  # subtelomeric window for the end/core ratio
TRACK_COLOR = "#238b45"
CAP_MARK_COLOR = "#00441b"
CEN_COLOR = "#c1272d"
UNRESOLVED_CHROM = "Pf3D7_10_v3"  # AT-content method fails here; open marker


def apply_style():
    """Publication defaults. Mirrors the figure-style skill's base rcParams so
    the script runs standalone."""
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"],
        "font.size": 7,
        "axes.titlesize": 8,
        "axes.labelsize": 7,
        "xtick.labelsize": 6,
        "ytick.labelsize": 6,
        "axes.linewidth": 0.6,
        "xtick.major.width": 0.6,
        "ytick.major.width": 0.6,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "figure.dpi": 150,
        "savefig.dpi": 300,
        "savefig.bbox": None,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    })


def nuclear_chroms(chrom_dict):
    """14 nuclear chromosomes in numeric order; drops apicoplast and mito."""
    nuc = {k: v for k, v in chrom_dict.items()
           if k.startswith("Pf3D7_") and "API" not in k}
    return dict(sorted(nuc.items(), key=lambda kv: int(kv[0].split("_")[1])))


def binned_mean(bw, chrom, length, bin_bp):
    """Mean signal per bin; BigWig gaps (None) become 0."""
    n = length // bin_bp
    vals = bw.stats(chrom, 0, n * bin_bp, nBins=n, type="mean")
    return np.array([0.0 if v is None else float(v) for v in vals])


def load_profiles(bw1_path, bw2_path, bin_bp=BIN_BP):
    """Average the two replicate tracks, then clamp negatives to zero."""
    bw1, bw2 = pyBigWig.open(bw1_path), pyBigWig.open(bw2_path)
    if bw1.chroms() != bw2.chroms():
        raise ValueError("BigWig headers differ; tracks are not on one assembly")
    nuc = nuclear_chroms(bw1.chroms())

    raw, prof = {}, {}
    for chrom, length in nuc.items():
        v1 = binned_mean(bw1, chrom, length, bin_bp)
        v2 = binned_mean(bw2, chrom, length, bin_bp)
        raw[chrom] = (v1 + v2) / 2.0          # average FIRST
        prof[chrom] = np.maximum(raw[chrom], 0)  # then clamp
    bw1.close()
    bw2.close()
    return nuc, raw, prof


def centromere_midpoints(csv_path, nuc):
    cen = pd.read_csv(csv_path)
    missing = set(nuc) - set(cen["chrom"])
    if missing:
        raise ValueError(f"centromere table is missing {sorted(missing)}")
    return {r.chrom: (r.start + r.end) / 2.0
            for r in cen.itertuples() if r.chrom in nuc}


def per_chromosome_stats(nuc, raw, prof, cen_mid, bin_bp=BIN_BP):
    """Means, medians, and the two length-independent enrichment ratios."""
    w = END_WINDOW_BP // bin_bp
    rows = []
    for chrom, length in nuc.items():
        a = prof[chrom]
        ends = np.concatenate([a[:w], a[-w:]])
        core = a[w:-w]
        ci = int(cen_mid[chrom] // bin_bp)
        cen_win = a[max(0, ci - 1):ci + 2]
        rows.append({
            "chromosome": chrom,
            "chr": int(chrom.split("_")[1]),
            "length_bp": length,
            "mean": float(a.mean()),
            "median": float(np.median(a)),
            "end_100kb_mean": float(ends.mean()),
            "core_mean": float(core.mean()),
            "end_over_core": float(ends.mean() / core.mean()),
            "cen_mean": float(cen_win.mean()),
            "cen_over_core": float(cen_win.mean() / core.mean()),
            "frac_bins_negative_before_clamp": float((raw[chrom] < 0).mean()),
        })
    return pd.DataFrame(rows)


def plot(nuc, prof, cen_mid, out_stem):
    all_vals = np.concatenate([prof[c] for c in nuc])
    y_cap = float(np.percentile(all_vals, CAP_PERCENTILE))
    max_mb = max(nuc.values()) / 1e6

    # explicit axes placement: 14 rows need predictable geometry, and
    # constrained_layout/tight_layout fight a stack this dense
    top, left, width = 0.930, 0.155, 0.815
    row_h, row_gap = 0.0545, 0.0075

    fig = plt.figure(figsize=(7.2, 5.6))
    axes = []
    for i, (chrom, length) in enumerate(nuc.items()):
        ax = fig.add_axes([left, top - (i + 1) * (row_h + row_gap), width, row_h])
        y = prof[chrom]
        x = (np.arange(len(y)) * BIN_BP + BIN_BP / 2) / 1e6

        ax.fill_between(x, 0, np.minimum(y, y_cap), color=TRACK_COLOR,
                        lw=0, step="mid")

        # bins exceeding the cap are flagged, never silently flattened
        over = y > y_cap
        if over.any():
            ax.plot(x[over], np.full(over.sum(), y_cap * 0.985), marker="^",
                    ms=1.8, ls="", color=CAP_MARK_COLOR, zorder=7)

        # centromere marker; hollow where the AT-content call is unresolved
        hollow = (chrom == UNRESOLVED_CHROM)
        ax.plot(cen_mid[chrom] / 1e6, y_cap, marker="v", ms=3.2,
                markerfacecolor=("none" if hollow else CEN_COLOR),
                markeredgecolor=CEN_COLOR, markeredgewidth=0.8,
                color=CEN_COLOR, clip_on=False, zorder=8)

        # baseline spans only the real chromosome, so row length reads as length
        ax.plot([0, length / 1e6], [0, 0], color="#999999", lw=0.7, zorder=2)

        ax.set_xlim(0, max_mb * 1.005)
        ax.set_ylim(0, y_cap)
        ax.set_xticks([])
        ax.set_yticks([])
        for side in ("top", "right", "bottom"):
            ax.spines[side].set_visible(False)
        ax.spines["left"].set_linewidth(0.5)

        ax.text(-0.055, 0.5, f"chr{i + 1:02d}", transform=ax.transAxes,
                ha="right", va="center", fontsize=7)
        ax.text(length / 1e6 + 0.04, 0.5 * y_cap, f"{length / 1e6:.2f} Mb",
                fontsize=5.5, color="#999999", va="center", ha="left")
        axes.append(ax)

    # one y-scale on the first row only; repeating it 14x is noise
    axes[0].set_yticks([0, 10])
    axes[0].tick_params(axis="y", labelsize=6, length=1.8, pad=1)

    fig.text(0.030, top - 7 * (row_h + row_gap), "Coverage (IP \u2212 IgG, CPM)",
             rotation=90, va="center", ha="center", fontsize=8)

    ax_last = axes[-1]
    ticks = np.arange(0, 3.6, 0.5)
    ax_last.set_xticks(ticks)
    ax_last.set_xticklabels([f"{v:.1f}" for v in ticks], fontsize=7)
    ax_last.tick_params(axis="x", length=2, pad=2)
    ax_last.spines["bottom"].set_visible(True)
    ax_last.spines["bottom"].set_linewidth(0.5)
    fig.text(left + width / 2, top - 14 * (row_h + row_gap) - 0.058,
             "Position along chromosome (Mb)", ha="center", fontsize=8.5)

    fig.text(0.030, 0.982, "CRWN-HA trophozoite (16 hpi) coverage per chromosome",
             fontsize=9.5, ha="left")
    fig.text(0.030, 0.963,
             "Mean of Rep1 C1 and C2 IP\u2212IgG tracks, negatives set to 0, "
             "5-kb bins. Shared axis 0\u2013%.1f CPM." % y_cap,
             fontsize=6.3, ha="left", color="#555555")

    # Marker key drawn with real markers, not unicode glyphs: the triangle
    # code points are absent from several common sans fonts and silently
    # render as missing-character boxes.
    key_y = 0.947
    fig.add_artist(Line2D([0.034], [key_y], marker="^", ms=2.6, ls="",
                          color=CAP_MARK_COLOR, transform=fig.transFigure))
    fig.text(0.046, key_y, "bin above axis cap", fontsize=6.3, ha="left",
             va="center", color="#555555")
    fig.add_artist(Line2D([0.215], [key_y], marker="v", ms=3.2, ls="",
                          color=CEN_COLOR, transform=fig.transFigure))
    fig.text(0.227, key_y, "centromere (AT-inferred)", fontsize=6.3, ha="left",
             va="center", color="#8a1f22")
    fig.add_artist(Line2D([0.406], [key_y], marker="v", ms=3.2, ls="",
                          markerfacecolor="none", markeredgecolor=CEN_COLOR,
                          markeredgewidth=0.8, color=CEN_COLOR,
                          transform=fig.transFigure))
    fig.text(0.418, key_y, "chr10 unresolved by this method", fontsize=6.3,
             ha="left", va="center", color="#8a1f22")

    fig.savefig(f"{out_stem}.png", dpi=300)
    fig.savefig(f"{out_stem}.pdf")
    return fig, y_cap, int((all_vals > y_cap).sum())


def check_no_text_overlaps(fig):
    """Geometric legibility check: flag colliding text, ignoring neighbouring
    tick labels (which abut by design)."""
    renderer = fig.canvas.get_renderer()
    items = [(t, t.get_window_extent(renderer))
             for t in fig.findobj(matplotlib.text.Text)
             if t.get_text().strip() and t.get_visible()]
    ticks = set()
    for ax in fig.axes:
        ticks |= set(ax.get_xticklabels() + ax.get_yticklabels())
    return [(a.get_text(), b.get_text())
            for i, (a, box_a) in enumerate(items)
            for b, box_b in items[i + 1:]
            if box_a.overlaps(box_b) and not (a in ticks and b in ticks)]


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bw1", required=True, help="Rep1 C1 IP-minus-IgG BigWig")
    ap.add_argument("--bw2", required=True, help="Rep1 C2 IP-minus-IgG BigWig")
    ap.add_argument("--centromeres", required=True,
                    help="Pf3D7_v48_centromeres.csv")
    ap.add_argument("--out", default="CRWN_Troph_coverage_per_chromosome",
                    help="output stem (.png and .pdf are written)")
    ap.add_argument("--stats-csv", default="CRWN_Troph_per_chromosome_coverage.csv")
    args = ap.parse_args()

    apply_style()
    nuc, raw, prof = load_profiles(args.bw1, args.bw2)
    cen_mid = centromere_midpoints(args.centromeres, nuc)

    stats = per_chromosome_stats(nuc, raw, prof, cen_mid)
    stats.round(4).to_csv(args.stats_csv, index=False)

    fig, y_cap, n_capped = plot(nuc, prof, cen_mid, args.out)

    overlaps = check_no_text_overlaps(fig)
    genome_mean = float(np.concatenate([prof[c] for c in nuc]).mean())
    print(f"chromosomes={len(nuc)} bins={sum(len(prof[c]) for c in nuc)}")
    print(f"genome_mean={genome_mean:.4f} y_cap={y_cap:.2f} capped_bins={n_capped}")
    print(f"mean_end_over_core={stats.end_over_core.mean():.2f} "
          f"mean_cen_over_core={stats.cen_over_core.mean():.2f}")
    print(f"text_overlaps={len(overlaps)}")


if __name__ == "__main__":
    main()
