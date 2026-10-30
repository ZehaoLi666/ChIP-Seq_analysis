#!/usr/bin/env python3
"""
Binomial enrichment tests + figure for CRWN-HA peak localization.

Consumes the counts printed by 05_heterochromatin_enrichment.sh and tests
each peak set against a genome-fraction background:
  - subtelomeric    : peaks within 50/100/200 kb of a chromosome end
  - VSA genes       : peaks overlapping var / rifin / stevor loci
  - centromeres     : peaks within +/-5 kb of an AT-derived centromere

Counts are hard-coded below exactly as the shell stage reported them, so the
tests are reproducible without re-running bedtools. Re-run 05_ and update
these if the peak sets change.

Caveats carried into the output:
  * overlap is counted per PEAK, not per base pair, which slightly favours
    small features;
  * the background is a plain genome fraction, NOT mappability-corrected --
    repetitive subtelomeres have lower effective mappability, which makes
    these enrichments conservative;
  * a permutation test against shuffled peaks would be the stronger control.

Outputs: CRWN_subtelomeric_enrichment_test.csv
         CRWN_VSA_gene_overlap_test.csv
         CRWN_centromere_overlap_test.csv
         CRWN_heterochromatin_enrichment.png / .pdf
"""

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

GENOME_BP = 23332839  # 14 nuclear chromosomes, PlasmoDB v48

# ---- observed counts, from 05_heterochromatin_enrichment.sh ----------------
SETS = ["all peaks (vs IgG)", "all peaks (vs Input)",
        "3-of-3 (vs IgG)", "3-of-3 (vs Input)"]
N_PEAKS = {"all peaks (vs IgG)": 2430, "all peaks (vs Input)": 2342,
           "3-of-3 (vs IgG)": 111,     "3-of-3 (vs Input)": 255}

# fraction of peaks within 50 / 100 / 200 kb of a chromosome end
SUBTEL_OBS = {
    "all peaks (vs IgG)":   {50: 0.110, 100: 0.185, 200: 0.282},
    "all peaks (vs Input)": {50: 0.127, 100: 0.193, 200: 0.289},
    "3-of-3 (vs IgG)":      {50: 0.027, 100: 0.072, 200: 0.135},
    "3-of-3 (vs Input)":    {50: 0.012, 100: 0.024, 200: 0.075},
}
SUBTEL_EXP = {50: 0.0617, 100: 0.1217, 200: 0.2417}

VSA_OVERLAP = {"all peaks (vs IgG)": 358, "all peaks (vs Input)": 325,
               "3-of-3 (vs IgG)": 1,      "3-of-3 (vs Input)": 1}
VSA_EXP = 0.0353          # 823,693 bp of var+rifin+stevor / genome

CEN_OVERLAP = {"all peaks (vs IgG)": 32, "all peaks (vs Input)": 30,
               "3-of-3 (vs IgG)": 23,    "3-of-3 (vs Input)": 29}
CEN_EXP = 174800 / GENOME_BP  # 14 centromeres +/- 5 kb = 0.75% of genome


def binom_table(obs_counts, expected_frac, label_col="peak_set"):
    rows = []
    for s in SETS:
        n = N_PEAKS[s]
        k = obs_counts[s]
        p = stats.binomtest(k, n, expected_frac, alternative="two-sided").pvalue
        rows.append({label_col: s, "n_peaks": n, "overlapping": k,
                     "observed_frac": round(k / n, 4),
                     "expected_frac": round(expected_frac, 4),
                     "enrichment": round((k / n) / expected_frac, 2),
                     "binom_p": f"{p:.2e}"})
    return pd.DataFrame(rows)


def subtelomeric_table():
    rows = []
    for s in SETS:
        n = N_PEAKS[s]
        for kb in (50, 100, 200):
            frac = SUBTEL_OBS[s][kb]
            k = int(round(frac * n))
            exp = SUBTEL_EXP[kb]
            p = stats.binomtest(k, n, exp, alternative="two-sided").pvalue
            rows.append({"set": s, "window_kb": kb, "observed_frac": frac,
                         "expected_frac": exp,
                         "enrichment": round(frac / exp, 2),
                         "binom_p": f"{p:.2e}"})
    return pd.DataFrame(rows)


def apply_style():
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"],
        "font.size": 7, "axes.titlesize": 8, "axes.labelsize": 7,
        "xtick.labelsize": 6, "ytick.labelsize": 6,
        "axes.linewidth": 0.6, "axes.spines.top": False,
        "axes.spines.right": False, "savefig.dpi": 300,
        "pdf.fonttype": 42, "ps.fonttype": 42,
    })


def open_frame(ax):
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    ax.spines["left"].set_position(("outward", 3))
    ax.spines["bottom"].set_position(("outward", 3))


def plot(subtel, vsa, out_stem="CRWN_heterochromatin_enrichment"):
    short = ["union\nvs IgG", "union\nvs Input", "3-of-3\nvs IgG", "3-of-3\nvs Input"]
    pivot = subtel.pivot(index="set", columns="window_kb",
                         values="enrichment").loc[SETS]
    vsa_e = vsa.set_index("peak_set").loc[SETS, "enrichment"].values

    fig = plt.figure(figsize=(7.4, 3.2))
    ax_a = fig.add_axes([0.085, 0.24, 0.38, 0.58])
    ax_b = fig.add_axes([0.605, 0.24, 0.345, 0.58])

    xs = np.arange(4)
    width = 0.25
    shades = {50: "#08519c", 100: "#4292c6", 200: "#9ecae1"}
    for i, kb in enumerate((50, 100, 200)):
        ax_a.bar(xs + (i - 1) * width, pivot[kb].values, width=width,
                 color=shades[kb], label=f"{kb} kb")
    ax_a.axhline(1.0, color="#333333", lw=0.9)
    ax_a.text(3.42, 1.06, "no enrichment", fontsize=6, color="#333333", ha="right")
    ax_a.set_xticks(xs); ax_a.set_xticklabels(short)
    ax_a.set_ylabel("Enrichment vs genome background")
    ax_a.set_title("Distance to chromosome end", loc="left")
    ax_a.set_ylim(0, 2.5)
    ax_a.legend(frameon=False, fontsize=6, title="window", title_fontsize=6,
                loc="upper right", handlelength=1.1, labelspacing=0.3)
    open_frame(ax_a)

    ax_b.bar(xs, vsa_e, width=0.55,
             color=["#c1272d", "#c1272d", "#bbbbbb", "#bbbbbb"])
    ax_b.axhline(1.0, color="#333333", lw=0.9)
    for i, v in enumerate(vsa_e):
        ax_b.text(i, v + 0.12, f"{v:.1f}\u00d7", ha="center", fontsize=6.5)
    ax_b.set_xticks(xs); ax_b.set_xticklabels(short)
    ax_b.set_ylabel("Enrichment vs genome background")
    ax_b.set_title("Overlap with $\\it{var}$/rifin/stevor genes", loc="left")
    ax_b.set_ylim(0, 5.0)
    open_frame(ax_b)

    for ax, letter in ((ax_a, "a"), (ax_b, "b")):
        ax.text(-0.17, 1.12, letter, transform=ax.transAxes,
                fontsize=9, fontweight="bold", va="top")

    fig.savefig(f"{out_stem}.png", dpi=300)
    fig.savefig(f"{out_stem}.pdf")
    return fig


def main():
    apply_style()
    subtel = subtelomeric_table()
    vsa = binom_table(VSA_OVERLAP, VSA_EXP)
    cen = binom_table(CEN_OVERLAP, CEN_EXP)

    subtel.to_csv("CRWN_subtelomeric_enrichment_test.csv", index=False)
    vsa.to_csv("CRWN_VSA_gene_overlap_test.csv", index=False)
    cen.to_csv("CRWN_centromere_overlap_test.csv", index=False)
    plot(subtel, vsa)

    print(subtel.to_string(index=False))
    print()
    print(vsa.to_string(index=False))
    print()
    print(cen.to_string(index=False))


if __name__ == "__main__":
    main()
