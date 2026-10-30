# CRWN-HA ChIP-seq — data processing code

All code used to process the CRWN-HA ChIP-seq data, from raw FASTQ to the
figures and report. Two execution contexts, per lab policy:

* `hpcc/` runs on UCR HPCC. Trimming + alignment go to SLURM as an array job;
  everything downstream runs on the standing interactive CPU node via
  `srun --jobid=<ID> --overlap`. No FASTQ or BAM leaves the cluster.
* `analysis/` and `figures/` run locally on exported BigWig / peak / stats
  files.

Reference throughout: *P. falciparum* 3D7, PlasmoDB v48, Bowtie2 index
`Pfalciparum_v48`. Effective genome size for MACS2: 2.3e7.

## Execution order

| Step | Script | Where | Produces |
|---|---|---|---|
| 1 | `hpcc/01_CRWN_chipseq_map.slurm` | SLURM array (5 tasks) | trimmed FASTQ → sorted / MAPQ30 / rmdup BAMs, `stats/*` |
| 2 | `hpcc/02_CRWN_qc_bigwig.sh` | interactive node | correlation + PCA, 5 CPM BigWigs, 12 IP-vs-control tracks |
| 3 | `hpcc/03_CRWN_peaks.sh` | interactive node | 12 MACS2 peak calls, FRiP, 3-of-3 + union BEDs, TSS/TES profiles |
| 4 | `hpcc/find_centromeres.py` → `cen_final.py` | interactive node | `Pf3D7_v48_centromeres.bed` |
| 5 | `analysis/05_heterochromatin_enrichment.sh` | interactive node | peak-distance and VSA/centromere overlap counts |
| 6 | `analysis/06_enrichment_tests.py` | local | binomial tests + `CRWN_heterochromatin_enrichment.png` |
| 7 | `analysis/07_average_replicates.py` | local | averaged BigWigs + per-chromosome coverage table |
| 8 | `figures/plot_*.py` | local | coverage-track and QC figures |
| 9 | `analysis/09_mapping_report_docx.py` | local | `CRWN_ChIPseq_mapping_report.docx` |

## Centromere coordinates — read before using

PlasmoDB v48 contains **no centromere features** (only gene/mRNA/CDS/exon), and
no centromere annotation existed anywhere on the cluster. The coordinates in
`Pf3D7_v48_centromeres.bed` are **derived from sequence**, not annotation:

`find_centromeres.py` scans 2-kb windows for maximum AT content, excluding
60 kb at each chromosome end so subtelomeric AT-rich repeats cannot win;
`cen_final.py` adds the requirement that the winning window sit in a ≥5-kb
gene desert. 13 of 14 chromosomes give a clean call (97–98% AT against a
~80% genome background, ~2.5 kb wide, consistent with the published PfCENH3
domain size).

**Chromosome 10 is not resolved by this method.** Its best 2-kb window reaches
only 91.5% AT with several near-ties, and a 4-kb scan matched to the PfCENH3
domain size gives 0.896 against a 0.892 runner-up (`cen10_4kb.py`). The region
where the coverage maximum falls (~930–940 kb) is **not** AT-rich (0.83–0.86),
has no assembly gap, and is flanked by ordinary metabolic genes
(`chr10_recheck.py`, `chr10_probe.py`) — so it is not a centromere by this
criterion either. Figures mark chr10 with a hollow triangle. Resolving it
requires published CenH3 coordinates or a CenH3 ChIP.

## Conventions worth preserving

* **Average before clamping.** `07_average_replicates.py` and
  `plot_CRWN_Troph_coverage.py` average replicates first and only then set
  negatives to zero, so opposing values between replicates cancel instead of
  inflating the baseline.
* **Per-chromosome means are length-confounded.** Mean coverage correlates with
  chromosome length at r = −0.95 because a sparse focal signal gets diluted
  over longer chromosomes. Interpret the length-independent ratios
  (`cen_vs_arm_ratio`, `end_over_core`), not `mean_coverage`.
* **Capped bins are marked, not flattened.** Figures that cap the y-axis draw a
  marker on every bin above the cap.
* **No unicode triangle glyphs in figure text.** Those code points are missing
  from several common sans fonts and render as blank boxes; marker keys are
  drawn with real markers instead.

## HPCC environment quirks baked into the scripts

* The `deeptools`, `gatk`, and `cutadapt` modulefiles are broken by a conda
  `--stack` incompatibility — tools are invoked by absolute path.
* `fastp` and bare `samtools` need explicit `LD_LIBRARY_PATH` entries
  (htslib 1.19.1, libdeflate 1.10, isa-l) or they exit 127 in seconds with an
  empty error log. Every binary is version-probed at job start so a shared
  library failure prints a message instead of failing silently.
* `bamCompare` refuses `--normalizeUsing CPM` unless paired with
  `--scaleFactorsMethod None`.
* `--regionsLabel` is not valid on `computeMatrix`; labelling happens at plot time.
* Control FASTQ filenames carry an extra plate-well token (`CRWN_IGG_C11_…`)
  that IP samples lack, so a uniform symlink loop silently produces broken
  links for exactly the controls. Every symlink is existence-tested.

## Verification status

Run against the saved artifacts in this project:

* `06_enrichment_tests.py` — reproduces the published enrichment table.
  One deliberate change: the centromere test now uses a **two-sided** binomial
  for consistency with the subtelomeric and VSA tests (originally one-sided),
  so its p-values are ~1.4× larger (e.g. 2.9e-03 vs 2.1e-03 for the IgG union
  set). Enrichment values are unchanged.
* `09_mapping_report_docx.py` — rebuilds the report (29 paragraphs, 1 table,
  1 image).
* `plot_genome_tracks.py` — reproduces the centromere coincidence result
  (13/14 chromosomes within 10 kb, median 1.9–4.2 kb across replicates).
* `plot_CRWN_Troph_coverage.py` — reproduces the trophozoite figure.
* `07_average_replicates.py` — runs ~7 min; not re-executed in the final pass.
  Its logic is unchanged from the run that produced the saved BigWigs.

## Dependencies

HPCC side: fastp 0.23.2, Bowtie2 2.5.3, samtools 1.19.2, bedtools 2.30.0,
MACS2 2.2.7.1, deepTools 3.5.1.
Local side: python ≥3.10 with `pyBigWig`, `numpy`, `pandas`, `scipy`,
`matplotlib`, `python-docx`.
