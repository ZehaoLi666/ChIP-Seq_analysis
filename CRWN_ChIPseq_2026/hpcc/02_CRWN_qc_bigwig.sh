#!/bin/bash
# CRWN ChIP-seq downstream stage 1: QC correlation/PCA + CPM BigWigs + IP-vs-control tracks
# RUN ON INTERACTIVE CPU NODE (srun --jobid=<CPU_JOBID> --overlap), not sbatch.
set -euo pipefail

PROJ=/rhome/zli529/lab/Aug2026_CRWN_ChIPseq
cd $PROJ
mkdir -p deeptools bigwig stats logs annotation

DT=/opt/linux/rocky/8.x/x86_64/pkgs/deeptools/3.5.1/bin
PKG=/opt/linux/rocky/8.x/x86_64/pkgs
export LD_LIBRARY_PATH=$PKG/htslib/1.19.1/lib:$PKG/libdeflate/1.10/lib:$PKG/isa-l/8b7c1b8/lib:${LD_LIBRARY_PATH:-}
SAMTOOLS=$PKG/samtools/1.19.2/bin/samtools

resolve () { command -v "$1" >/dev/null 2>&1 && command -v "$1" || echo "$DT/$1"; }
MULTIBAMSUMMARY=$(resolve multiBamSummary)
MULTIBIGWIGSUMMARY=$(resolve multiBigwigSummary)
PLOTCORRELATION=$(resolve plotCorrelation)
PLOTPCA=$(resolve plotPCA)
BAMCOVERAGE=$(resolve bamCoverage)
BAMCOMPARE=$(resolve bamCompare)
for t in "$MULTIBAMSUMMARY" "$MULTIBIGWIGSUMMARY" "$PLOTCORRELATION" "$PLOTPCA" "$BAMCOVERAGE" "$BAMCOMPARE"; do
  [[ -x "$t" ]] || { echo "ERROR: $t not executable"; exit 1; }
done
echo "deepTools resolved: $BAMCOVERAGE"

NPROC=2
SAMPLES=$(cat samples.txt)
for S in $SAMPLES; do
  [[ -s bam/${S}.rmdup.bam ]] || { echo "ERROR: missing bam/${S}.rmdup.bam"; exit 1; }
  $SAMTOOLS quickcheck bam/${S}.rmdup.bam || { echo "ERROR: truncated bam ${S}"; exit 1; }
done

# ---------- 1. BAM-level correlation + PCA ----------
echo "[$(date)] multiBamSummary"
$MULTIBAMSUMMARY bins \
  --bamfiles bam/CRWN_HA_D1.rmdup.bam bam/CRWN_HA_D2.rmdup.bam bam/CRWN_HA_D3.rmdup.bam \
             bam/CRWN_IGG.rmdup.bam bam/CRWN_INPUT.rmdup.bam \
  --labels HA_D1 HA_D2 HA_D3 IgG Input \
  --binSize 10000 --numberOfProcessors $NPROC \
  --outFileName deeptools/mapping_bam_summary.npz \
  --outRawCounts deeptools/mapping_bam_summary_counts.tsv
[[ -s deeptools/mapping_bam_summary.npz ]] || { echo "ERROR: bam summary missing"; exit 1; }

for M in pearson spearman; do
  $PLOTCORRELATION -in deeptools/mapping_bam_summary.npz \
    --corMethod $M --whatToPlot heatmap --skipZeros --removeOutliers \
    --plotNumbers \
    --plotTitle "${M^} correlation of CRWN ChIP-seq BAMs" \
    --plotFile deeptools/mapping_${M}_correlation_heatmap.pdf \
    --outFileCorMatrix deeptools/mapping_${M}_correlation_matrix.tsv
  [[ -s deeptools/mapping_${M}_correlation_matrix.tsv ]] || { echo "ERROR: $M matrix missing"; exit 1; }
done

$PLOTPCA -in deeptools/mapping_bam_summary.npz \
  -o deeptools/mapping_PCA.pdf \
  --outFileNameData deeptools/mapping_PCA_data.tsv \
  --plotTitle "PCA of CRWN ChIP-seq BAMs"
[[ -s deeptools/mapping_PCA.pdf ]] || { echo "ERROR: PCA missing"; exit 1; }
echo "[$(date)] BAM QC done"

# ---------- 2. CPM BigWigs, every sample ----------
for S in $SAMPLES; do
  echo "[$(date)] bamCoverage $S"
  $BAMCOVERAGE --bam bam/${S}.rmdup.bam \
    --outFileName bigwig/${S}.CPM.bw --outFileFormat bigwig \
    --normalizeUsing CPM --binSize 10 --numberOfProcessors $NPROC --extendReads
  [[ -s bigwig/${S}.CPM.bw ]] || { echo "ERROR: CPM bw failed $S"; exit 1; }
done

# ---------- 3. BigWig-level correlation + PCA ----------
echo "[$(date)] multiBigwigSummary"
$MULTIBIGWIGSUMMARY bins \
  --bwfiles bigwig/CRWN_HA_D1.CPM.bw bigwig/CRWN_HA_D2.CPM.bw bigwig/CRWN_HA_D3.CPM.bw \
            bigwig/CRWN_IGG.CPM.bw bigwig/CRWN_INPUT.CPM.bw \
  --labels HA_D1 HA_D2 HA_D3 IgG Input \
  --binSize 10000 --numberOfProcessors $NPROC \
  --outFileName deeptools/bigwig_CPM_summary.npz \
  --outRawCounts deeptools/bigwig_CPM_summary_counts.tsv
[[ -s deeptools/bigwig_CPM_summary.npz ]] || { echo "ERROR: bw summary missing"; exit 1; }

for M in pearson spearman; do
  $PLOTCORRELATION -in deeptools/bigwig_CPM_summary.npz \
    --corMethod $M --whatToPlot heatmap --skipZeros --removeOutliers --plotNumbers \
    --plotTitle "${M^} correlation of CRWN CPM BigWigs" \
    --plotFile deeptools/bigwig_CPM_${M}_correlation_heatmap.pdf \
    --outFileCorMatrix deeptools/bigwig_CPM_${M}_correlation_matrix.tsv
done
$PLOTPCA -in deeptools/bigwig_CPM_summary.npz \
  -o deeptools/bigwig_CPM_PCA.pdf \
  --outFileNameData deeptools/bigwig_CPM_PCA_data.tsv \
  --plotTitle "PCA of CRWN CPM BigWigs"
echo "[$(date)] BigWig QC done"

# ---------- 4. IP vs control tracks (both IgG and Input) ----------
# NOTE: bamCompare requires --scaleFactorsMethod None alongside --normalizeUsing CPM
for IP in CRWN_HA_D1 CRWN_HA_D2 CRWN_HA_D3; do
  for CTRL in CRWN_IGG CRWN_INPUT; do
    echo "[$(date)] bamCompare $IP vs $CTRL"
    $BAMCOMPARE --bamfile1 bam/${IP}.rmdup.bam --bamfile2 bam/${CTRL}.rmdup.bam \
      --outFileName bigwig/${IP}_vs_${CTRL}.log2ratio.bw --outFileFormat bigwig \
      --operation log2 --pseudocount 1 \
      --scaleFactorsMethod None --normalizeUsing CPM \
      --binSize 10 --numberOfProcessors $NPROC
    [[ -s bigwig/${IP}_vs_${CTRL}.log2ratio.bw ]] || { echo "ERROR: log2 failed $IP/$CTRL"; exit 1; }

    $BAMCOMPARE --bamfile1 bam/${IP}.rmdup.bam --bamfile2 bam/${CTRL}.rmdup.bam \
      --outFileName bigwig/${IP}_minus_${CTRL}.subtract.bw --outFileFormat bigwig \
      --operation subtract \
      --scaleFactorsMethod None --normalizeUsing CPM \
      --binSize 10 --numberOfProcessors $NPROC
    [[ -s bigwig/${IP}_minus_${CTRL}.subtract.bw ]] || { echo "ERROR: subtract failed $IP/$CTRL"; exit 1; }
  done
done

echo "[$(date)] STAGE2_COMPLETE"
