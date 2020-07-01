# Part 2 - Pre-processing and quality control (QC) of raw data 

### Learning objectives: 
- 


## Read pre-processing & trimming

### Principles of read trimming: Downstream steps are more efficient
Read trimming is used to clean up the library of raw reads. Trimming can be used to trim adapter sequences, polyA tails, low quality bases, or reads that are too short. Mostly trimming is used to remove low quality bases, generally 20 is the lowest filter we would recommend for quality trimming. 


![Read alignment](../figures/read_processing.png)


### Example command for read trimming with cutadapt 
(If you use Cutadapt, please cite DOI:10.14806/ej.17.1.200)

Usage:
    
cutadapt -a ADAPTER [options] [-o output.fastq] input.fastq

For paired-end reads:

cutadapt -a ADAPT1 -A ADAPT2 [options] -o out1.fastq -p out2.fastq in1.fastq in2.fastq

```bash
cutadapt -a 'A{76}' -m 5 --nextseq-trim=20 -o "$i".trimmed.fq.gz "$i"_R1_001.fastq.gz >"$i".cutadapt.out;
```
-a trims adapters from the 3' end
-g trims adapters from the 5' end
-b trims adapters from both ends
-m trims reads that are samller than the minimum threshold (length is measured after all quality and adapter trimming)
-q qulaity threshold for trimming bases

--nextseq-trim works like q=20 except that quality scores on Gs are ignored, this accomodates the effects of dark cycles from some illumina instruments that leave a string of high quality but in correct Gs at the 3' end of reads

-o output file
