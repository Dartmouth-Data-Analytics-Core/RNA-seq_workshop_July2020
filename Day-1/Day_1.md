## Day 1 

### Pre-processing and quality control of raw data 
- Basic quality control of FASTQ files 
- Go over a basic FASTQC report 
- Generate a quick FASTQC report 

<br>

### Basics 
- Read number required for experiment  
- Single-end vs paired-end 
- Full-length vs 3'-end assays 
- Replicates & experimental design 
- Structure of a for loop

<br>

### Read pre-processing & trimming  
#### Principles of read trimming: Downstream steps are more efficient
Read trimming is used to clean up the library of raw reads. Trimming can be used to trim adapter sequences, polyA tails, low quality bases, or reads that are too short. Mostly trimming is used to remove low quality bases, the quality score of a base (Q score) denotes the probability that a base was called correctly. The higher the Q score the more likely the base was called correctly. A Q score of 30 indicates a 99.9% chance that the base call was accurate, a Q score of 20 indicates a 99% changes of accuracy,etc. Generally 20 is the lowest filter we would recommend for quality trimming. 

#### Example command for read trimming with cutadapt 
(If you use Cutadapt, please cite DOI:10.14806/ej.17.1.200)
Usage:
    cutadapt -a ADAPTER [options] [-o output.fastq] input.fastq

For paired-end reads:
    cutadapt -a ADAPT1 -A ADAPT2 [options] -o out1.fastq -p out2.fastq in1.fastq in2.fastq
    
  ```bash
  for i in `cat prefix.list`;
  do cutadapt -a 'A{76}' -m 5 --nextseq-trim=20 -o "$i".trimmed.fq.gz "$i"_R1_001.fastq.gz >"$i".cutadapt.out;
  done
  ```



### Read alignment  
- Principles of gapped short read alignment (mentioning splice junctions)
- Perform alignment with STAR
- Read clipping 
- View & explore some reads in IGV (show difference on read distributions in full-length transcript & 3'end data)
- Quick mention of quasi-mapping with tools like salmon 

### Post-alignment quality control 
- Picard tools CollectRNASeqMetrics 
- Identiofy PCR duplicates, discuss value of checking, and controversey of removing them 
- RNA-seq QC metrics 
- Discuss a bad QC report 
- Using MultiQC to synthesize a QC report 

<br>

### Read count quantification  
- Principles of read counting 
- Quantification with HTSeq-count 
- Mention of more complex probablistic methods (e.g. RSEM)

<br>

### Exercises? (if time, could be take home) 
Just ideas...

1. All of your samples show good quality alignment %s except for one, which has an alignment rate of 55%?  
You noted when doing QC of the raw FASTQ files that this sample has many duplicated sequences.  
What might these sequences represent and how could you quickly test this?





