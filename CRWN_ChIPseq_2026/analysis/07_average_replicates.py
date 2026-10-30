#!/usr/bin/env python3
"""
Average two IP-minus-control BigWig tracks and summarize coverage per chromosome.

Produces the averaged tracks used for the D1+D2 coverage figure:
  <out>_avg.bw      averaged, negatives retained
  <out>_avg_pos.bw  averaged, negatives clamped to zero

ORDER MATTERS: averaging happens FIRST, clamping SECOND. A bin where one
replicate is +2 and the other -2 then averages to 0 rather than to +1;
clamping before averaging would inflate the baseline.

Per-chromosome statistics include the two length-independent ratios
(centromere/arm and end/core). The plain per-chromosome MEAN is reported but
is dominated by chromosome length -- see the note in the printed output.

Usage:
  python 07_average_replicates.py \
      --bw1 CRWN_HA_D1_minus_CRWN_INPUT.subtract.bw \
      --bw2 CRWN_HA_D2_minus_CRWN_INPUT.subtract.bw \
      --centromeres Pf3D7_v48_centromeres.csv \
      --out CRWN_HA_rep1rep2_minus_Input
"""

import argparse

import numpy as np
import pandas as pd
import pyBigWig
from scipy import stats as st

STEP = 10            # native resolution of the output BigWig, bp
DISPLAY_BIN = 5000   # bin size for the per-chromosome summary
END_WINDOW = 100000  # subtelomeric window for the end/core ratio


def nuclear_chroms(chrom_dict):
    nuc = {k: v for k, v in chrom_dict.items()
           if k.startswith("Pf3D7_") and "API" not in k}
    return dict(sorted(nuc.items(), key=lambda kv: int(kv[0].split("_")[1])))


def average_tracks(bw1_path, bw2_path, step=STEP):
    bw1, bw2 = pyBigWig.open(bw1_path), pyBigWig.open(bw2_path)
    hdr = bw1.chroms()
    if hdr != bw2.chroms():
        raise ValueError("BigWig headers differ; not the same assembly")
    avg = {}
    for chrom, length in hdr.items():
        n = length // step
        v1 = np.nan_to_num(np.array(
            bw1.stats(chrom, 0, n * step, nBins=n, type="mean"), dtype=float))
        v2 = np.nan_to_num(np.array(
            bw2.stats(chrom, 0, n * step, nBins=n, type="mean"), dtype=float))
        avg[chrom] = ((v1 + v2) / 2.0).astype(np.float32)   # average FIRST
    bw1.close(); bw2.close()
    return hdr, avg


def write_bigwig(path, hdr, values, step=STEP):
    out = pyBigWig.open(path, "w")
    out.addHeader([(c, hdr[c]) for c in hdr])
    for chrom in hdr:
        out.addEntries(chrom, 0,
                       values=values[chrom].astype(float).tolist(),
                       span=step, step=step)
    out.close()


def rebin(arr, step, bin_bp):
    """Collapse step-bp values into bin_bp means (drops a ragged tail bin)."""
    k = bin_bp // step
    n = len(arr) // k
    return arr[:n * k].reshape(n, k).mean(axis=1)


def per_chromosome_stats(nuc, clipped, cen_mid, step=STEP, bin_bp=DISPLAY_BIN):
    rows = []
    for chrom, length in nuc.items():
        a = clipped[chrom]
        ci0, ci1 = int(cen_mid[chrom] // step) - 125, int(cen_mid[chrom] // step) + 125
        ci0 = max(0, ci0)
        mask = np.ones(len(a), bool)
        mask[ci0:ci1] = False
        w = END_WINDOW // step
        ends = np.concatenate([a[:w], a[-w:]])
        core = a[w:-w]
        rows.append({
            "chromosome": chrom,
            "chr": int(chrom.split("_")[1]),
            "length_bp": length,
            "mean_coverage": float(a.mean()),
            "median_coverage": float(np.median(a)),
            "p99_coverage": float(np.percentile(a, 99)),
            "max_coverage": float(a.max()),
            "centromere_mean": float(a[ci0:ci1].mean()),
            "non_centromere_mean": float(a[mask].mean()),
            "end_100kb_mean": float(ends.mean()),
            "core_mean": float(core.mean()),
        })
    df = pd.DataFrame(rows)
    df["cen_vs_arm_ratio"] = (df.centromere_mean / df.non_centromere_mean).round(1)
    df["end_over_core"] = (df.end_100kb_mean / df.core_mean).round(3)
    return df


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bw1", required=True)
    ap.add_argument("--bw2", required=True)
    ap.add_argument("--centromeres", required=True)
    ap.add_argument("--out", default="CRWN_HA_rep1rep2_minus_Input")
    ap.add_argument("--stats-csv", default="CRWN_per_chromosome_coverage.csv")
    args = ap.parse_args()

    hdr, avg = average_tracks(args.bw1, args.bw2)
    clipped = {c: np.maximum(a, 0).astype(np.float32) for c, a in avg.items()}

    write_bigwig(f"{args.out}_avg.bw", hdr, avg)
    write_bigwig(f"{args.out}_avg_pos.bw", hdr, clipped)

    nuc = nuclear_chroms(hdr)
    cen = pd.read_csv(args.centromeres)
    cen_mid = {r.chrom: (r.start + r.end) / 2.0
               for r in cen.itertuples() if r.chrom in nuc}

    stats_df = per_chromosome_stats(nuc, clipped, cen_mid)
    stats_df.round(4).to_csv(args.stats_csv, index=False)

    genome_mean = float(np.concatenate([clipped[c] for c in nuc]).mean())
    neg_frac = float(np.mean([(avg[c] < 0).mean() for c in nuc]))
    r_p, p_p = st.pearsonr(stats_df.length_bp, stats_df.mean_coverage)
    r_s, p_s = st.spearmanr(stats_df.length_bp, stats_df.mean_coverage)

    print(stats_df[["chr", "length_bp", "mean_coverage", "median_coverage",
                    "centromere_mean", "cen_vs_arm_ratio",
                    "end_over_core"]].round(3).to_string(index=False))
    print(f"\ngenome_mean={genome_mean:.4f} "
          f"frac_bins_negative_before_clamp={neg_frac:.3f}")
    print(f"length vs mean coverage: Pearson r={r_p:.3f} p={p_p:.2e} | "
          f"Spearman rho={r_s:.3f} p={p_s:.2e}")
    print("NOTE: the length correlation is an artefact of averaging sparse "
          "focal signal over chromosomes of different size -- interpret "
          "cen_vs_arm_ratio and end_over_core, not mean_coverage.")


if __name__ == "__main__":
    main()
