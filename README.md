# ChIP-Seq-analysis

Hello everyone! This is a guaidline for ChIP-seq analuysis on Linux. 

# Introduction 
## ChIP-seq 
### A few differences compared to RNA_seq analysis 
- Don't need to trim bases at the beginning of each read (no bias)
- Can use a simpler aligner like bowtie2 because splicing doesn't need to be accounted for
- The IGV tracks are very important and need to make normalized tracks to carefully look at where the peaks are.
- Definitely need to do removal of reads marked as PCR duplicates (we've have some debate on this for RNA-seq)

## Methods and Materials 

# Part1: Quality control 
## Using FastQC on cluster 
```
module load fastqc 
fastqc *.fastq   # to run all fastq files in this directory
```

# Part2: Read mapping with Bowtie2
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

# Part3: Alignment
In this part, we can use SAMTools.
LINKS TO SAMTOOLS PAGE:



