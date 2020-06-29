### Pre-processing and quality control of raw data 

Before you begin working with genomic data it is important to asses the quality of the data, look for possible contamination, and get an idea of what you can expect the assembled product to look like. The distribution of base qualities, the distribution of read lengths, GC bias, and duplication rates are all informative metrics for assessing raw sequencing data. 

FastQC runs several 'analysis' modules on the raw FASTQ files that allow us to evaluate the quality of that sample, and identify potential quality issues that need to be addressed in the downstream analysis steps. Many short reads, skewed GC content, and high duplication rates can indicate contamination by adapters, bacteria, etc. 

```bash
fastqc infile.fq.gz --outdir=fastqc_out
```
Lets have a look at a typical HTML report. 

[Good Illumina Data FastQC Report](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/good_sequence_short_fastqc.html)

[Bad Illumina Data FastQC Report](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/bad_sequence_fastqc.html)


<br>
 

Openning 1 HTML file at a time to do this is annoying, so we use a tool called MultiQC to aggregate these reports together. 
```bash
multiqc .
```
MultiQC is very flexible and comprehensive, and will process output from many bioinformatics programs (FastqQC, cutadapt, STAR, picard, samtools), so we can use it to do more QA/QC and collect key metrics later in the pre-processing steps. More on this later. 
