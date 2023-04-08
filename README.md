# ChIP-Seq-analysis

Hello everyone! This is a guaidline for ChIP-seq analuysis on Linux. 

# Introduction 
## ChIP-seq 
## Methods and Materials 

# Part1: Quality control 
## Using FastQC on cluster 
```
module load fastqc 
fastqc *.fastq   # to run all fastq files in this directory
```
# Part2: Read mapping with Bowtie2
```
wget                              # download genome file 
bowtie2-build input.fasta output  # build reference genome next read mapping
bowtie2 -x referencegenome -1 output-forward-paired.fq.gz -2 output-reverse-paired.fq.gz -S output-paired-aligned.sam  
```


```
sdgdsfkhgkldfsg

```
