#!/bin/bash
# CRWN ChIP-seq downstream stage 2: peak calling (narrow+broad x IgG+Input), FRiP,
# replicate consistency, TSS/TES profiles, summary table.
# RUN ON INTERACTIVE CPU NODE (srun --jobid=<CPU_JOBID> --overlap), not sbatch.
set -euo pipefail

PROJ=/rhome/zli529/lab/Aug2026_CRWN_ChIPseq
cd $PROJ
mkdir -p peaks stats deeptools annotation logs

PKG=/opt/linux/rocky/8.x/x86_64/pkgs
DT=$PKG/deeptools/3.5.1/bin
export LD_LIBRARY_PATH=$PKG/htslib/1.19.1/lib:$PKG/libdeflate/1.10/lib:$PKG/isa-l/8b7c1b8/lib:${LD_LIBRARY_PATH:-}
SAMTOOLS=$PKG/samtools/1.19.2/bin/samtools
BEDTOOLS=$PKG/bedtools/2.30.0/bin/bedtools
GTF=/rhome/zli529/lab/PlasmoDB_Genome/PlasmoDB_v48/PlasmoDB-48_Pfalciparum3D7.gtf

resolve () { command -v "$1" >/dev/null 2>&1 && command -v "$1" || echo "$DT/$1"; }
COMPUTEMATRIX=$(resolve computeMatrix)
PLOTPROFILE=$(resolve plotProfile)
PLOTHEATMAP=$(resolve plotHeatmap)

if command -v macs2 >/dev/null 2>&1; then MACS=$(command -v macs2)
elif command -v macs3 >/dev/null 2>&1; then MACS=$(command -v macs3)
elif [[ -x $PKG/macs2/2.2.7/bin/macs2 ]]; then MACS=$PKG/macs2/2.2.7/bin/macs2
else echo "ERROR: no macs2/macs3 found"; exit 1; fi
"$MACS" --version || { echo "ERROR: macs not runnable"; exit 1; }
echo "Using peak caller: $MACS"

GSIZE=2.3e7

# ---------- 0. genes.bed from GTF ----------
if [[ ! -s annotation/genes.bed ]]; then
  awk 'BEGIN{OFS="\t"} $3=="gene" {
        id=""; if (match($0,/gene_id "[^"]+"/)) { id=substr($0,RSTART+9,RLENGTH-10) }
        if (id=="") id="gene"NR
        print $1, $4-1, $5, id, ".", $7 }' "$GTF" \
    | sort -k1,1 -k2,2n > annotation/genes.bed
fi
[[ -s annotation/genes.bed ]] || { echo "ERROR: genes.bed empty"; exit 1; }
echo "genes.bed lines: $(wc -l < annotation/genes.bed)"

# ---------- 1. peak calling: 3 IPs x 2 controls x 2 modes ----------
for IP in CRWN_HA_D1 CRWN_HA_D2 CRWN_HA_D3; do
  for CTRL in CRWN_IGG CRWN_INPUT; do
    CTAG=${CTRL#CRWN_}

    # narrow
    OUT=peaks/${IP}_vs_${CTAG}_narrow
    mkdir -p $OUT
    echo "[$(date)] MACS narrow $IP vs $CTAG"
    "$MACS" callpeak -t bam/${IP}.rmdup.bam -c bam/${CTRL}.rmdup.bam \
      -f BAMPE -g $GSIZE -n ${IP}_vs_${CTAG} --outdir $OUT -q 0.01 \
      > logs/macs_${IP}_vs_${CTAG}_narrow.log 2>&1 || echo "WARNING: narrow MACS nonzero exit $IP/$CTAG"
    NP=$OUT/${IP}_vs_${CTAG}_peaks.narrowPeak
    [[ -s $NP ]] && echo "  narrow peaks: $(wc -l < $NP)" || echo "  WARNING: no narrow peaks $IP/$CTAG"

    # broad
    OUTB=peaks/${IP}_vs_${CTAG}_broad
    mkdir -p $OUTB
    echo "[$(date)] MACS broad $IP vs $CTAG"
    "$MACS" callpeak -t bam/${IP}.rmdup.bam -c bam/${CTRL}.rmdup.bam \
      -f BAMPE -g $GSIZE -n ${IP}_vs_${CTAG} --outdir $OUTB \
      --broad --broad-cutoff 0.05 -q 0.05 \
      > logs/macs_${IP}_vs_${CTAG}_broad.log 2>&1 || echo "WARNING: broad MACS nonzero exit $IP/$CTAG"
    BP=$OUTB/${IP}_vs_${CTAG}_peaks.broadPeak
    [[ -s $BP ]] && echo "  broad peaks: $(wc -l < $BP)" || echo "  WARNING: no broad peaks $IP/$CTAG"
  done
done

# ---------- 2. FRiP for every peak set ----------
echo -e "IP\tcontrol\tmode\tpeak_reads\ttotal_reads\tFRiP\tn_peaks" > stats/CRWN_FRiP.tsv
for IP in CRWN_HA_D1 CRWN_HA_D2 CRWN_HA_D3; do
  TOTAL=$($SAMTOOLS view -c bam/${IP}.rmdup.bam)
  for CTAG in IGG INPUT; do
    for MODE in narrow broad; do
      if [[ $MODE == narrow ]]; then PK=peaks/${IP}_vs_${CTAG}_narrow/${IP}_vs_${CTAG}_peaks.narrowPeak
      else PK=peaks/${IP}_vs_${CTAG}_broad/${IP}_vs_${CTAG}_peaks.broadPeak; fi
      if [[ -s $PK ]]; then
        NPK=$(wc -l < $PK)
        PR=$($BEDTOOLS intersect -a bam/${IP}.rmdup.bam -b $PK -u 2>/dev/null | $SAMTOOLS view -c -)
        awk -v i=$IP -v c=$CTAG -v m=$MODE -v p=$PR -v t=$TOTAL -v n=$NPK \
          'BEGIN{OFS="\t"; print i,c,m,p,t,(t>0?p/t:"NA"),n}' >> stats/CRWN_FRiP.tsv
      else
        echo -e "${IP}\t${CTAG}\t${MODE}\tNA\t${TOTAL}\tNA\t0" >> stats/CRWN_FRiP.tsv
      fi
    done
  done
done
cat stats/CRWN_FRiP.tsv

# ---------- 3. replicate consistency (per control x mode) ----------
for CTAG in IGG INPUT; do
  for MODE in narrow broad; do
    if [[ $MODE == narrow ]]; then EXT=narrowPeak; else EXT=broadPeak; fi
    P1=peaks/CRWN_HA_D1_vs_${CTAG}_${MODE}/CRWN_HA_D1_vs_${CTAG}_peaks.${EXT}
    P2=peaks/CRWN_HA_D2_vs_${CTAG}_${MODE}/CRWN_HA_D2_vs_${CTAG}_peaks.${EXT}
    P3=peaks/CRWN_HA_D3_vs_${CTAG}_${MODE}/CRWN_HA_D3_vs_${CTAG}_peaks.${EXT}
    if [[ -s $P1 && -s $P2 && -s $P3 ]]; then
      # peaks present in all three replicates
      $BEDTOOLS intersect -a $P1 -b $P2 -u | $BEDTOOLS intersect -a - -b $P3 -u \
        | sort -k1,1 -k2,2n > peaks/CRWN_HA_${CTAG}_${MODE}.reproducible_3of3.bed
      cat $P1 $P2 $P3 | sort -k1,1 -k2,2n | $BEDTOOLS merge \
        > peaks/CRWN_HA_${CTAG}_${MODE}.merged_union.bed
      echo "${CTAG} ${MODE}: 3of3=$(wc -l < peaks/CRWN_HA_${CTAG}_${MODE}.reproducible_3of3.bed) union=$(wc -l < peaks/CRWN_HA_${CTAG}_${MODE}.merged_union.bed)"
    else
      echo "WARNING: incomplete peak sets for ${CTAG} ${MODE}, skipping consistency"
    fi
  done
done

# ---------- 4. TSS / TES profiles from log2 tracks ----------
BWS=""
for IP in CRWN_HA_D1 CRWN_HA_D2 CRWN_HA_D3; do
  [[ -s bigwig/${IP}_vs_CRWN_IGG.log2ratio.bw ]] && BWS="$BWS bigwig/${IP}_vs_CRWN_IGG.log2ratio.bw"
done
if [[ -n "$BWS" ]]; then
  for RP in TSS TES; do
    echo "[$(date)] computeMatrix $RP"
    $COMPUTEMATRIX reference-point --referencePoint $RP -b 1000 -a 1000 --binSize 10 \
      -R annotation/genes.bed -S $BWS --samplesLabel HA_D1 HA_D2 HA_D3 \
      --skipZeros --missingDataAsZero --numberOfProcessors 2 \
      -o deeptools/CRWN_HA_vs_IgG_${RP}_matrix.gz
    $PLOTPROFILE -m deeptools/CRWN_HA_vs_IgG_${RP}_matrix.gz \
      -out deeptools/CRWN_HA_vs_IgG_${RP}_profile.pdf \
      --outFileNameData deeptools/CRWN_HA_vs_IgG_${RP}_profile_data.tsv \
      --plotTitle "CRWN-HA (log2 IP/IgG) around $RP" \
      --yAxisLabel "log2 IP/IgG" --dpi 600
    $PLOTHEATMAP -m deeptools/CRWN_HA_vs_IgG_${RP}_matrix.gz \
      -out deeptools/CRWN_HA_vs_IgG_${RP}_heatmap.pdf \
      --plotTitle "CRWN-HA around $RP" --colorMap Reds --dpi 600
  done
fi

# ---------- 5. mapping + summary table ----------
echo -e "Sample\tTotal_pairs_in\tOverall_align_rate\tMAPQ30_dup_removed_reads\tDup_rate" > stats/CRWN_mapping_summary.tsv
for S in $(cat samples.txt); do
  ALN=$(grep -m1 'overall alignment rate' stats/${S}.bowtie2.log 2>/dev/null | awk '{print $1}' || echo NA)
  TOTIN=$(grep -m1 'reads; of these' stats/${S}.bowtie2.log 2>/dev/null | awk '{print $1}' || echo NA)
  FINAL=$(cat stats/${S}.final_read_count.txt 2>/dev/null || echo NA)
  DUP=$(grep -m1 'DUPLICATE TOTAL' stats/${S}.markdup.txt 2>/dev/null | awk '{print $NF}' || echo NA)
  echo -e "${S}\t${TOTIN}\t${ALN}\t${FINAL}\t${DUP}" >> stats/CRWN_mapping_summary.tsv
done
cat stats/CRWN_mapping_summary.tsv

echo "[$(date)] STAGE3_COMPLETE"
