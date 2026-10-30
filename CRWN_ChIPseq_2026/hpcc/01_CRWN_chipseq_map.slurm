#!/bin/bash
#SBATCH --job-name=CRWN_chip_map
#SBATCH --account=lerochlab
#SBATCH --partition=epyc
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=60G
#SBATCH --time=48:00:00
#SBATCH --output=logs/CRWN_chip_map_%A_%a.out
#SBATCH --error=logs/CRWN_chip_map_%A_%a.err
#SBATCH --array=1-5

set -euo pipefail

PROJ=/rhome/zli529/lab/Aug2026_CRWN_ChIPseq
cd $PROJ
mkdir -p clean_fastq fastqc bam stats logs

BOWTIE2_INDEX=/rhome/zli529/lab/PlasmoDB_Genome/PlasmoDB_v48/bowtie2_index/Pfalciparum_v48

PKG=/opt/linux/rocky/8.x/x86_64/pkgs
export LD_LIBRARY_PATH=$PKG/htslib/1.19.1/lib:$PKG/libdeflate/1.10/lib:$PKG/isa-l/8b7c1b8/lib:${LD_LIBRARY_PATH:-}
FASTP=$PKG/fastp/0.23.2/bin/fastp
BOWTIE2=$PKG/bowtie2/2.5.3/bowtie2
SAMTOOLS=$PKG/samtools/1.19.2/bin/samtools

# fail loudly and early if shared libs are wrong (silent 127 otherwise)
for b in "$FASTP" "$BOWTIE2" "$SAMTOOLS"; do
  "$b" --version >/dev/null 2>&1 || { echo "ERROR: $b not runnable (shared libs?)"; exit 1; }
done

SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" samples.txt)
echo "[$(date)] START ${SAMPLE} on $(hostname)"

# ---------- 1. fastp trimming ----------
"$FASTP" \
  -i raw_data/${SAMPLE}_R1.fastq.gz \
  -I raw_data/${SAMPLE}_R2.fastq.gz \
  -o clean_fastq/${SAMPLE}_R1.clean.fastq.gz \
  -O clean_fastq/${SAMPLE}_R2.clean.fastq.gz \
  -h fastqc/${SAMPLE}.fastp.html \
  -j fastqc/${SAMPLE}.fastp.json \
  -w 8 2> stats/${SAMPLE}.fastp.log

[[ -s clean_fastq/${SAMPLE}_R1.clean.fastq.gz && -s clean_fastq/${SAMPLE}_R2.clean.fastq.gz ]] \
  || { echo "ERROR: fastp output missing for ${SAMPLE}"; exit 1; }
echo "[$(date)] fastp done ${SAMPLE}"

# ---------- 2. bowtie2 alignment ----------
"$BOWTIE2" \
  -x ${BOWTIE2_INDEX} \
  -1 clean_fastq/${SAMPLE}_R1.clean.fastq.gz \
  -2 clean_fastq/${SAMPLE}_R2.clean.fastq.gz \
  -p 8 \
  --very-sensitive \
  2> stats/${SAMPLE}.bowtie2.log \
  | "$SAMTOOLS" sort -@ 4 -m 2G -T bam/${SAMPLE}.sorttmp -o bam/${SAMPLE}.sorted.bam
"$SAMTOOLS" index -@ 4 bam/${SAMPLE}.sorted.bam

[[ -s bam/${SAMPLE}.sorted.bam && -s bam/${SAMPLE}.sorted.bam.bai ]] \
  || { echo "ERROR: sorted BAM missing for ${SAMPLE}"; exit 1; }
"$SAMTOOLS" quickcheck bam/${SAMPLE}.sorted.bam || { echo "ERROR: BAM truncated ${SAMPLE}"; exit 1; }
echo "[$(date)] bowtie2 done ${SAMPLE}"

# ---------- 3. MAPQ30 filter (keep properly paired) ----------
"$SAMTOOLS" view -@ 8 -b -q 30 -f 2 -F 3844 \
  bam/${SAMPLE}.sorted.bam \
  -o bam/${SAMPLE}.mapq30.bam
"$SAMTOOLS" sort -@ 4 -m 2G -T bam/${SAMPLE}.q30tmp -o bam/${SAMPLE}.mapq30.sorted.bam bam/${SAMPLE}.mapq30.bam
"$SAMTOOLS" index -@ 4 bam/${SAMPLE}.mapq30.sorted.bam
rm -f bam/${SAMPLE}.mapq30.bam
[[ -s bam/${SAMPLE}.mapq30.sorted.bam ]] || { echo "ERROR: MAPQ BAM missing ${SAMPLE}"; exit 1; }
echo "[$(date)] mapq30 done ${SAMPLE}"

# ---------- 4. duplicate removal (paired-end) ----------
"$SAMTOOLS" sort -@ 4 -m 2G -n -T bam/${SAMPLE}.nstmp -o bam/${SAMPLE}.namesort.bam bam/${SAMPLE}.mapq30.sorted.bam
"$SAMTOOLS" fixmate -@ 4 -m bam/${SAMPLE}.namesort.bam bam/${SAMPLE}.fixmate.bam
"$SAMTOOLS" sort -@ 4 -m 2G -T bam/${SAMPLE}.pstmp -o bam/${SAMPLE}.positionsort.bam bam/${SAMPLE}.fixmate.bam
"$SAMTOOLS" markdup -@ 4 -r -s bam/${SAMPLE}.positionsort.bam bam/${SAMPLE}.rmdup.bam 2> stats/${SAMPLE}.markdup.txt
"$SAMTOOLS" index -@ 4 bam/${SAMPLE}.rmdup.bam

[[ -s bam/${SAMPLE}.rmdup.bam && -s bam/${SAMPLE}.rmdup.bam.bai ]] \
  || { echo "ERROR: rmdup BAM missing ${SAMPLE}"; exit 1; }
rm -f bam/${SAMPLE}.namesort.bam bam/${SAMPLE}.fixmate.bam bam/${SAMPLE}.positionsort.bam
echo "[$(date)] markdup done ${SAMPLE}"

# ---------- 5. stats ----------
"$SAMTOOLS" flagstat bam/${SAMPLE}.rmdup.bam > stats/${SAMPLE}.flagstat.txt
"$SAMTOOLS" idxstats bam/${SAMPLE}.rmdup.bam > stats/${SAMPLE}.idxstats.txt
"$SAMTOOLS" stats    bam/${SAMPLE}.rmdup.bam > stats/${SAMPLE}.samtools.stats.txt
"$SAMTOOLS" view -c  bam/${SAMPLE}.rmdup.bam > stats/${SAMPLE}.final_read_count.txt

# ---------- 6. cleanup temp fastq ----------
rm -f clean_fastq/${SAMPLE}_R1.clean.fastq.gz clean_fastq/${SAMPLE}_R2.clean.fastq.gz
rm -f bam/${SAMPLE}.sorttmp*.bam bam/${SAMPLE}.q30tmp*.bam bam/${SAMPLE}.nstmp*.bam bam/${SAMPLE}.pstmp*.bam

echo "[$(date)] FINISHED ${SAMPLE}"
