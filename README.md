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
### 1.2.1 Methods 
#### Workflow 
#### Pachages 
#### Automatic analysis 
### 1.2.2 Materials


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
wget URL 
grep '\trRNA\t' referencegenome.gff > rRNA.bed
```
### 4.3.2 Use BEDtools to intersect alignment BAM file with rRNA.bed file and filter out reads to rRNA genes
```
module load bedtools 
bedtools intersect -abam $1-sorted-removed.bam -b rRNA.bed -v > $1-sorted-filteredrRNA.bam 
```
## 4.4 Keep only properly paired reads
ChIP-seq experiments typically involve sequencing DNA fragments that have been immunoprecipitated with a specific antibody, followed by mapping the resulting sequencing reads to the reference genome. The immunoprecipitation process can introduce biases in the sequencing library preparation, leading to non-specific binding, PCR amplification artifacts, and other sources of noise.

Filtering out improperly paired reads can help to reduce the effects of these biases and improve the accuracy of downstream analyses. Improperly paired reads are reads where one or both reads in a paired-end sequencing experiment do not align as expected based on the insert size distribution, indicating that the read pair is likely to be a PCR artifact or other source of noise.
```
samtools view -b -f 2 -F 4 output-removed-sorted-paired-aligned.bam > output-only-properly-paired-aligned.bam 
```

# 5. Part4 Visualize the alignment in IGV (genome browser)
This step in very important. We need to make TDF or BigWig file by converting BAM file. To directly compare samples on IGV, normalize by number of million mapped reads before producing TDF/BigWig. One possible way of doing this: Convert BAM to BED file, use bedtools genomecov with -d option to get read count per nucleotide in the genome, convert to WIG with custom Python script, convert WIG to TDF using Igvtools ToTDF. Can then visualize the TDF. 

## 5.1 Make referencegenome.genome file
```
samtools faidx referencegenome.fasta   # faidx creates an index file of the fasta file for quick access sequence elements
```

## 5.2 Convert BAM to BED file
```
bedtools bamtobed -i output-only-properly-paired-aligned.bam > output-only-properly-paired-aligned.bed 
``` 

## 5.3 Get read count per nucleotide in the genome
```
bedtools genomecov -d -i output-only-properly-paired-aligned.bed -g referencegenome.genome > output.txt 
```

## 5.4 Convert txt file to Wig and TDF format
Use [MMR_normalize_by_nucleotide.py](https://github.com/ZehaoLi666/ChIP-Seq_analysis/blob/main/MMR_normalize_by_nucleotide.py) to do normalization.
There are three parameters we need to input: input file name, out name file name and MMR_number.
MMR-number can get from read counts you get from the last bam file from samtools. If read counts is 5 millions, MMR_number is 5. 
```
samtool view -c output-only-properly-paired-aligned.bam
python3.8 MMR_normalize_by_nucleotide.py input.txt output.txt MMR_number 

```
use [make_wig_from_nucleotide_coverage.py](https://github.com/ZehaoLi666/ChIP-Seq_analysis/blob/main/make_wig_from_nucleotide_coverage.py) to convert txt file to Wig format
```
python3.8 make_wig_from_nucleotide_coverage.py input.txt output.wig 
```
use igvtools to convert Wig format to TDF format
```
module load igv
igvtools toTDF output-only-properly-paired-norm.wig output-only-properly-paired-norm.tdf referencegenome.fasta.fai 
```
Here are some results from IGV: 
In the example we gave, we were trying to find two proteins binding area in Plasmodium berghei's genome. 
![results of chip-seq analysis on IGV](https://github.com/ZehaoLi666/ChIP-Seq_analysis/blob/main/IGV-chip.png)
It's a little difficult to find peaks by naked eyes. Actually, this step can only provide us a general preview. We need more sophisticated statistic analysis in the next step.



# 6. Part5 Peak calling 
We can use [MACS2](https://hbctraining.github.io/Intro-to-ChIPseq/lessons/05_peak_calling_macs.html) to do this step.
For many cases, the chip-seq would be much more complicated, so if we directly visualize the alignment file in IGV, we may not get what we want. 
So we need to follow the workflow below. 
![Chip-seq analysis workflow with Masc2](https://hbctraining.github.io/Intro-to-ChIPseq/img/chip_workflow_june2017_step2.png)

Because we've already combined sense.bam and antisense.bam (R1.bam and R2.bam) together, but Macs2 needs them to be separate. Like this:
```
module load macs2    
module load gcc 
macs2 callpeak -t sample1_R1.bam -c control_R1.bam -f BAM -g 1.9e+7 -n KIN5_G1 --outdir macs2 
macs2 callpeak -t sample2_R2.bam -c control_R2.bam -f BAM -g 1.9e+7 -n KIN5_G1 --outdir macs2
# -t input bam file 
# -c control bam file
# -f type of input file 
# -g genome size 
# -n output prefix name
# --ourdir  directory of output files
```

From the previous steps, we've already combined sense.bam and antisense.bam (R1.bam and R2.bam) together, we need add paraments like this:
```
macs2 callpeak -t ~/chip/3.alignment/4_keep_paired/paired_KIN5_G1.bam  -c ~/chip/3.alignment/4_keep_paired/paired_GAM.bam -f BAM -g 1.9e+7 -n KIN5_G1 --outdir macs2 --nomodel --extsize 147     
# --nomodel, this option will skip the model building step
# --extsize 147, extsize option specifies the extension size for fragment length estimation.
```
