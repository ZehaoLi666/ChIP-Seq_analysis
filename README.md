# ChIP-Seq-analysis

Hello everyone! This is a guaidline for ChIP-seq analuysis on Linux. 

# 1. Introduction 
## 1.1 ChIP-seq 
### 1.1.1 A few differences compared to RNA_seq analysis 
- Don't need to trim bases at the beginning of each read (no bias)
- Can use a simpler aligner like bowtie2 because splicing doesn't need to be accounted for
- The IGV tracks are very important and need to make normalized tracks to carefully look at where the peaks are.
- Definitely need to do removal of reads marked as PCR duplicates (we've have some debate on this for RNA-seq)

## 1.2 Methods and Materials 

# 2. Part1: Quality control 
## 2.1 Using FastQC on cluster 
```
module load fastqc 
fastqc *.fastq   # to run all fastq files in this directory
```

# 3. Part2: Read mapping with Bowtie2
way1: Interactive sessions with srun, if this job will not run for a long time.
```
wget                              # download genome file 
bowtie2-build input.fasta reference.genome  # build reference genome next read mapping, reference.genome is output file name, you can name putput file name as you want
bowtie2 -x reference.genome -1 output-forward-paired.fq -2 output-reverse-paired.fq -S output-paired-aligned.sam  
```

way2:Job submission with sbatch, if this job will run for a long time.

```
#!/bin/bash -l

#SBATCH --nodes=6
#SBATCH --ntasks=6
#SBATCH --cpus-per-task=6
#SBATCH --mem-per-cpu=10G
#SBATCH --time=0-03:00:00  
#SBATCH --mail-user=zli529@ucr.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name="zli_bowtie"
#SBATCH --output=bowtie2.log
#SBATCH -p batch 

module load bowtie2
bowtie2 -x reference  -1 ~/chip/data/GAM_Input_S25_L004_R1_001.fastq -2 ~/chip/data/GAM_Input_S25_L004_R2_001.fastq -S GAM.sam 

```

# 4. Part3: Alignment
In this part, we can use SAMTools.
LINKS TO SAMTOOLS PAGE:

## 4.1 remove reads with mapping quality below 10. (optional)
This is step is optional. In some cases, filtering out low-quality reads may be necessary to obtain reliable results, particularly in situations where the signal-to-noise ratio is low or the sequencing depth is low. However, in other cases, retaining all reads regardless of mapping quality may be necessary to capture rare events or low-abundance features. Therefore, whether doing this step depends on research purpose and specific condition. We can cmompare the results of doing this step and not doing this step, to select the most suitbale one for your study.

If used Bowtie2, removing reads with mapping quality below 10 is a good idea, for others scores are assigned differently, can examine mapping quality distribution using BamQC.

```
samtools view -b -q 10 output-sorted-paired-aligned.bam > output-removed-sorted-paired-aligned.bam 
```

## 4.2 Sort the aligment file.
Sorting the alignment file is an important step in ChIP-seq data analysis because it allows for more efficient downstream processing and analysis. When sequencing reads are mapped to the genome, the resulting alignment file is typically unsorted, meaning that the reads are ordered arbitrarily. Sorting the alignment file by genomic coordinates places the reads in order according to their position on the genome, which enables more efficient processing of the data.

```
module load samtools 
samtools sort -o sorted.bam  input.bam
```
We chould also submit job with sbatch.

## 4.3 Filter out reads to rRNA genes
rRNA genes are highly expressed and can account for a large proportion of sequencing reads, which can reduce the effective sequencing depth for other genomic regions of interest. 

### 4.3.1 Creat a bed file with the rRNA genes 

Download the referencegenome.gff file from database. 
```
grep '\trRNA\t' referencegenome.gff > rRNA.bed
```
### 4.3.2 Use BEDtools to intersect alignment BAM file with rRNA.bed file and filter out reads to rRNA genes
```
module load bedtools 
bedtools intersect -abam $1-sorted-removed.bam -b rRNA.bed -v > $1-sorted-filteredrRNA.bam 
```

