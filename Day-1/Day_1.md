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

<br>

### Read pre-processing & trimming  
- Principles of read trimming: makes downstream steps more efficient
- Performing read trimming with cutadapt


### Read alignment  
- Principles of gapped short read alignment (mentioning splice junctions)
- Perform alignment with STAR

### Post-alignment quality control 
- Picard tools CollectRNASeqMetrics 
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





