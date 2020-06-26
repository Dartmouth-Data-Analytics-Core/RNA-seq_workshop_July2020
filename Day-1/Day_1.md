## Day 1 



### FASTQ file format

FASTQ files are arguably the workhorse format of bioinformatics. FASTQs are used to store sequence reads generated in next-generatoon sequencing (NGS) experiments. Similarly to FASTA files, FASTQ files contain a herder line, followed by the sequence read, however individual quality of base calls from the sequencer are included for each record in a FASTQ file. 

Here is what a the first record of an example FASTQ file looks like
```
@SRR1039508.1 HWI-ST177:290:C0TECACXX:1:1101:1225:2130 length=63
CATTGCTGATACCAANNNNNNNNGCATTCCTCAAGGTCTTCCTCCTTCCCTTACGGAATTACA
+
HJJJJJJJJJJJJJJ########00?GHIJJJJJJJIJJJJJJJJJJJJJJJJJHHHFFFFFD
```

**Four rows exist for each record in a FASTQ file:**
- Row 1: Header line that stores information about the read (always starts with an `@`), such as the *instrument ID*, *flowcell ID*, *lane on flowcell*, *file number*, *cluster coordinates*, *sample barcode*, etc.
- Row 2: The sequence if bases called
- Row 3: Usually just a `+` and sometimes followed by the read info. in line 1
- Row 4: Individual base qualities (must be same length as line 2

Quality scores, also known as **Phred scores**, in row 4 represent the probability that the associated base call is incorrect, which are defined by the below formaula for current Illumina machines:
```
Q = -10 x log10(P), where Q = base quality, P = probability of incorrect base call

or 

P = 10^-Q/10
```

Intuitively, this means that a base with a Phred score of `10` has a `1 in 10` chance of being an incorrectly called base, or *90%*. Likewise, a score of `20` has a `1 in 100` chance (99% accuracy), `30` a `1 in 1000` chance (99.9%) and `40` a `1 in 10,000` chance (99.99%). 

However, we can clearly see that these are not probabilities. Instead, quality scores are encoded by a character that is associated with an *ASCII* code (equal to the *Phred-score +33*). The reason for doing it this way is so that quality scores only take up 1 byte per value in the FASTQ file. 

For example, the first base call in our sequence example above, the `C` has a quality score encoded by an `H`, which corresponds to a Q-score of 39, meaning this is a good quality base call. 

Generally, you can see this would be a good quality read if not for the strech of `#`s indicating a Q-score of 2. Looking at the FASTQ record, you can see these correspond to a string of `N` calls, which are bases that the sequencer was not able to make a base call for. Streches of Ns' are generally not useful for your analysis. 

You can read more about quality score encoding [here](https://support.illumina.com/help/BaseSpace_OLH_009008/Content/Source/Informatics/BS/QualityScoreEncoding_swBS.htm), and view the full table of symbols and *ASCII* codes used to represent Q-scores. 

**Paired-end reads:**  

If you sequenced paired-end reads, you will have two FASTQ files:  
*..._R1.fastq* - contains the forward reads  
*..._R2.fastq*- contains the reverse reads  

Most downstream analysis tools will recognize that such files are paired-end, and the reads in the forward file correspond to the reads in the reverse file, although you often have to specify the names of both files to these tools. 

It is critical that the R1 and R2 files have **the same number of records in both files**. If one has more records than the other, which can sometimes happen if there was an issue in the demultiplexing process, you will experience problems using these files as paired-end reads in downstream analyses. 

### Working with FASTQ files 

While you don't normally need to go looking within an individual FASTQ file, it is very important to be able to manipulate FASTQ files in you are going to be doing any more involved bioinformatics. There are a lot of operations we can do with a FASTQ file to gain more information about our experiment, and being able to interact with FASTQ files can be useful for troubleshooting problems that might come up in your analyses. 

Due to their large size, we often perform gzip copmpression of FASTQ files so that they take up less space, however this means we have to unzip them if we want to look inside them and perform operations on them. We can do this with the `zcat` command. 

Lets use `zcat` and `head` to have a look at the first few records in our FASTQ file. 
```bash
zcat sample.fastq.gz | head
```

How many records do we have in total? (don't forget to divide by 4..) 
```bash
zcat sample.fastq.gz | wc -l
```

How do we capture specific lines of each record recursively?
```bash
zcat sample.fastq.gz | wc -l
```

What if we want to count how many unique barcodes exist in the FASTQ file. 
```bash
zcat sample.fastq.gz | sed -n '2~4p' | head -10000 | grep -o .
```
Using sed with the -n option and '2~4p' will return the 2nd line, skip to 4 lines, and print again, recursively. We can use head to view this for the first 10,000 lines. 

Next, we can use grep to find all of the instances of each nucleotide across these lines. 
```bash
zcat sample.fastq.gz | sed -n '2~4p' | head -10000 | grep -o .
```

Finally we can sort and count up the unique nucleotides that were found..
```bash
zcat sample.fastq.gz | sed -n '2~4p' | head -10000 | grep -o . | sort | uniq -c
```
Now we have the number of each nuleotide across the reads from the first 1000 records. A quick and easy program to get GC content! 

### Searching for regular expressions

What if we wanted to find matches for a specific sequence, maybe after a start codon, in the FASTQ file 
```bash
zcat sample.fastq.gz | grep -o "ATGGGATCA" | sort | uniq -c
```

Perhaps this sequence represents some a contaminating sequence from the run that we want to quickly screen all of our samples for (e.g. from bacteria). We can do this by searching for matches and counting how many times it was found, and repeating this process for each sample using a for loop. 
```bash
ls *.fastq.gz | 
for read x; do 
    echo $x is being processed...; zcat $x | grep -o "ATGGGATCA" | sort | uniq -c; 
done
```
TO DO, need to get a real sequence, or change the exercise a little
Could we use an adapter sequence to check if the adapters need to be trimmed? I can't imagine a small string that would represent bacterial contamination, but an adapter sequence would be pretty appropriate and lead right into trimming.

- Running this using an executable shell script

- 1 Example of same stuff with a FASTA file 

- Running in the background with nohup 

The value of these sorts of tasks may not be immediately clear, but as we start piping together these operations we can calculate useful metrics and gain some basic insight into the reads in the file as a whole. Such operations are used by common bioinformatics tools under the hood to perform routine tasks.

<br>

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


<br>

### Read pre-processing & trimming  
#### Principles of read trimming: Downstream steps are more efficient
Read trimming is used to clean up the library of raw reads. Trimming can be used to trim adapter sequences, polyA tails, low quality bases, or reads that are too short. Mostly trimming is used to remove low quality bases, generally 20 is the lowest filter we would recommend for quality trimming. 


![Read alignment](../figures/read_processing.png)


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
-a trims adapters from the 3' end

-g trims adapters from the 5' end

-b trims adapters from both ends

-m trims reads that are samller than the minimum threshold (length is measured after all quality and adapter trimming)

-q qulaity threshold for trimming bases

--nextseq-trim works like q=20 except that quality scores on Gs are ignored, this accomodates the effects of dark cycles from some illumina instruments that leave a string of high quality but in correct Gs at the 3' end of reads

-o output file

<br>


![Read alignment](../figures/read_alignment.png)




### Read alignment  

Aligning millions of reads to very large reference genomes (such as the human genome) is generally done by splitting the reads and reference into a catalog of shorter reads with unique sequnece structures (kmers). It is improtant when selecting an alignement program to ensure that it is appropriate for the dataset you are working with, for example STAR (Spliced Transcripts Alignment to a Reference) is used to align reads that have come from spliced transcripts. A single read from a spliced transcriptome might map across a splice junction, such that the left side of the read and the right side of the read map hundreds of base pairs apart. If your dataset is prokaryotic (non-splicosomal) this would not be the appropriate program for you to align your reads, we would suggest looking into bwa-mem or bowtie2. If you are in a hurry and not interested in obtaining read alignments and only need count data quasi-mapping with a tool like Salmon might be a good option. 

#### STAR read aligner
STAR uses a method of seed searching, clustering, stitching, and scoring to find the most probable match in the reference sequence for each read. A seed is the longest possible match between a read and the reference sequence. By using multiple seeds on a single read, reads are able to span hundreds of base pairs across splice junctions. Once a read is mapped into multiple seeds STAR attempts to map the remaining unmapped portions of the read by extending the seed match allowing for indels and mismatches. Any portion of the read that cannot be mapped is assumed to be contamination, leftover adapter sequences, or an incorrect base call and these bases are clipped (called soft-clipping).


```bash
STAR --genomeDir myind --sjdbGTFfile mygene --runThreadN 4 --outSAMunmapped Within --outFilterType BySJout --outSAMattributes NH HI AS NM MD --outSAMtype BAM SortedByCoordinate --outFilterMultimapNmax 10 --outFilterMismatchNmax 999 --outFilterMismatchNoverReadLmax 0.04 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --readFilesIn trimed_R1_fastq --readFilesCommand zcat --outFileNamePrefix Sample_ID
```

--genomeDir the path to the directory with genome indices

--sjdbGTFfile the path to the annotated transcripts in GTF file format

--runThreadN number of threads to use in the run

--outSAMunmapped the name of the file to write unmapped reads to 

--outFilterType how mapped reads will be filtered (normal/BySJout)

--outSAMattributes the attributes to write to the outfile with mapped reads

--outSAMtype (BAM unsorted/ BAM SortedByCoordinate)

--outFilterMultimapNmax maximum number of alignments allowed for a read, reads that exceed this will be considered unmapped

--outFilterMismatchNmax maximum number of mismatches per pair, 999 turns this filter off

--outFilterMismatchNoverReadLmax maximum number of mismatches per pair relative to read length, 0.4 with 100bp allows for 8 mismatches across the pair

--alignIntronMin minimum intron length

--alignIntronMax maximum intron length

--alignMatesGapMax maximum genomic distance between mates

--alignSJoverhangMin minimum overhang for unannotated junctions

--alignSJBDoverhangMin minimum overhang for annotated junctions

--readFilesIn read files to map to reference alignment

--readFilesCommand uncompression command to apply to read files

--outFileNamePrefix prefix for outfiles generated in the run


- View & explore some reads in IGV (show difference on read distributions in full-length transcript & 3'end data)

<br>

## Viewing Alignments in IGV

The Integrative Genomics Viewer (IGV) from the Broad Institute is an extremely useful tool for visulazation of alignment files (as well as other genomic file formats). Viewing your alignments in this way can be used to explore your data, troubleshoot issues you are having downstream of alignment, and inspect coverage for and quality of reads in specific regions of interest (e.g. in variant calling). I strongly encourage you to download the IGV for your computer from their [website](http://software.broadinstitute.org/software/igv/) and play around with some BAM file to get familar with all its various features. 

Here, we will create a small subset of a BAM file, download it onto our local machines, and view it using the IGV web app (for speed). You can open the IGV web app in your browser [here](https://igv.org/app/). 

Lets go ahead and subset our BAM file for reads aligning only to chromosome 22. We also need to create an index. 
```bash
# subset for reads just on chr 22 (to make it smaller)
samtools view -b -@ 8 -o chr22.bam SRR1039508_1.Aligned.sortedByCoord.out.bam 22

# index your new bam file 
samtools index chr22.bam
``` 

Now download the files onto your local machine, so that you can load them into the IGV web app. 
```bash
# make a directory and go into it 
mkdir rnaseq_wrksp/
cd rnaseq_wrksp/

# download the file using secure copy (scp)
##### modify this for your discovery ID, AND the directory your in 
scp d41294d@discovery7.dartmouth.edu:/dartfs-hpc/rc/lab/G/GMBSR_bioinfo/workshops/rnaseq-july20/data/bam/chr22.bam* .
##### you will be promoted for your password for discovery/polaris 

# you may also need to modify the permissions 
chmod a+rwx chr22*
```

Now navigate to the IGV web app, and follow the below steps:  
1. Under *'Genomes'* in the top right, select *'GRCh38/hg38'*
2. Under *'Tracks'* select *'Local file'*, then highlight both the .bam and .bai 
3. Navigate to chromosome 22, and zoom in on ...
4. Use the setting at the side of the track to set colors by *read strand*

![IGV](../figures/igv_webapp.png)

- Can you find an example of a highly expressed gene?
- What do you notice about the orientation of the aligning reads? 
- Why might this cause ambiguity in read quantification for genes that overlap?  


### Post-alignment quality control 
Once you have your reads aligned you need to assess the quality of your alignment. Before running FastQC or MultiQC it is prudent to get more information about the aligned files. Picard has some useful tools for assesing the quality of an alignment. 

#### CollectRNASeqMetrics 
This tool will generate a file with several metrics relating to the quality of the alignment. These metrics can be used to diagnose problems with the library that can either be mitigated or taken into consideratoin in assessing the alignment of the data. Metrics reported include the dsitribution of bases within transcripts, the total number and fraction of nucleotides within genomic regions (UTR, introns, exons, intragenic regions, etc.), the number of bases that pass wuality filters, median coverage, the ratio of 5' to 3' biases, and the number of reads designated to the correct strand. 

```bash 
java -jar picard.jar CollectRnaSeqMetrics I=input.bam O=output.RNA_Metrics REF_FLAT=ref_flat.txt STRAND=SECOND_READ_TRANSCRIPTION_STRAND RIBOSOMAL_INTERVALS=ribosomal.interval_list
      
 ```
 I= input aligned bam file
 
 O= output RNAseq metrics
 
 REF_FLAT= Gene annotations in refFlat form. Format described here: http://genome.ucsc.edu/goldenPath/gbdDescriptionsOld.html#RefFlat
 
 STRAND= For strand-specific library prep. For unpaired reads, use FIRST_READ_TRANSCRIPTION_STRAND if the reads are expected to be on the transcription strand.
 
 RIBOSOMAL_INTERVALS= Location of rRNA sequences in genome, in interval_list format. If not specified no bases will be identified as being ribosomal. Format described here: http://samtools.github.io/htsjdk/javadoc/htsjdk/htsjdk/samtools/util/IntervalList.html
 
#### MarkDuplicates
This tool detects duplicate reads that might have been generated during library preparation or the sequencing run. Duplicates are reads that originate from a single piece of DNA. Duplicate reads are idnetified by comparing the sequences in the 5' end of both reads and read pairs. A high rate of duplicate reads left in your library will skew the count matrix you generate downstream, and ultimately affect the differential expression analysis. If one or more of your libraries have a high proportion of duplicate reads it may be worth removing them before generating a count matrix. 
```bash 
java -Xmx32G -jar picard.jar MarkDuplicates I=myInBam O=myOutBam M=myTxt OPTICAL_DUPLICATE_PIXEL_DISTANCE=100 CREATE_INDEX=false
```
I= input aligned bam file

O= output with duplicate reads marked

M= file to write duplication metrics to

OPTICAL_DUPLICATE_PIXEL_DISTANCE= The maximum offset between two duplicate clusters in order to consider them optical duplicates

CREATE_INDEX= (TRUE/FALSE) Whether to create a BAM index when writing a coordinate-sorted BAM file.

- RNA-seq QC metrics 
- Discuss a bad QC report 
- Using MultiQC to synthesize a QC report 

<br>

### Duplicate removal 



## Read count quantification  

For most downstream analyses in RNA-seq, especially differential expression, we care about how many reads aligned to a specific gene, as this tells us about its expression level, which we can then compare to other samples. Inherently, this means that we essentially want to make these data *count-based*, so that we can use statistical models to compare these counts between our experimental conditions of interest. 

![](../figures/quant_principle.png)

There are a number of methods that one can use to quantify reads overlapping a specific features ( in this case, exons). These methods require a .bam file as input, in addition to a list of the genomic features that we wish to count reads over. The most simplistic methods (e.g. [htseq-count](https://htseq.readthedocs.io/en/release_0.11.1/count.html), [featureCounts](http://subread.sourceforge.net/)) use a specific set of rules to count the number of reads overlapping specific features. These are a good choice if your data is less complex, e.g. 3'-end data. 

Other methods leverage probablistic modeling in order to quantify the alignments (e.g. [RSEM](https://deweylab.github.io/RSEM/)), which ascibes reads to features with a probability of this being the correct location for a given read. Generally, these methods are used on more complex data (e.g. full-length transcript and/or paired-end data) where transcript estimates are of interest. 

For our analysis, we will used htseq-count. htseq-count has three distinct modes for handling overlapping features; union, intersection_strict, and intersection_nonempty. You change change these using the `mode` option. The behavious of each setting is described in the figure below, taken from the htseq-count documnetation. Another important (and related) behaviour of htseq-count is how it handles multi-mapping reads (reads that map to > 1 place in the genome). Generally, multimapping reads are discarded during counting to prevent introduction of bias into the counts. 

![](htseq-count-mode.png)

One of the most important options in htseq-count is `strandedness`. It is critical to select the correct option for `strandedness` (`-s`) for your dataset, otherwise you may incorrectly use, or throw away, a lot of information. The default setting in htseq-count for `strandedness` is `yes`. This means reads will only be counted as overlapping a feature provided they map to the same strand as the feature. If you data was generated using an unstranded library preparation protocol, as in this experiment, we must set this option to `no`. Failiure to do so would mean you would throw away ~50% of all your reads, as they will be distributed equally across both strands for each feature in an unstranded library.  

Another important important option in htseq-count is `t` or `type` which specifies which feature type (3rd column of a GTF file) you want to count features over. The default is `exon` which works for GTF files from Ensembl, such as the file we will be using. However, this can be changed to any feature in your GTF file, so theoretically can be used to count any feature you have annotated. 

When counting paired-end data (such as in this experiemnt) your .bam files should be sorted before running htseq-count, and you can specify how your .bam is sorted using the `-r` option. `name` indicates they are sorted by read name, `pos` indicates they are sorted by genomic position. 

# run htseq-count on your .bam file 
```bash
htseq-count \
	-f bam \
	-s no \
	-r pos \
	./pipeline.out/alignment/${i}Aligned.sortedByCoord.out.bam \
	/dartfs-hpc/rc/lab/G/GMBSR_bioinfo/genomic_references/human/ensembl-annotations/Homo_sapiens.GRCh38.97.gtf > ${i}.htseq-counts
```

There are numerous settings that can be tweaked and turned on/off in htseq-count. I strongly recommend you **read the manual** before running htseq-count so that you understand all the default options and available settings. 

.... Let it run...

Lets have a look at the resulting file. 
```bash
# first few rows 
head out.htseq.counts

# importantly, lets check the last few rows as these contain some important info 
tail -n 12 out.htseq.counts
```

Exercise: 
- Can you visulaly confirm the read count returned in htseq-count by looking at the .bam file in IGV? 

# Generate the gene expression matrix of raw read counts

The final step in the pre-processing of RNA-seq data for differential expression analysis is to concatenate your read counts into a gene expression matrix that contains the counts from all your samples. We will do this at the command line, however there are also ways to directly read the output of programs like htseq-count and RSEM directly into R without concatenating them into a matrix before hand. 

Have a look at the htseq-count output files 
```bash
ls -1 *.htseq-counts | sort
```

Loop over htseq-count output files and extract the read count column 
```bash
# set up an array that we will fill with shorthand sample names
myarray=()

# loop over count files using a while loop 
while read x;  do 
	# split up sample names to remove everything after "_trim"
	sname=`echo bash misc/split.sh "$x" "_trim"`
	# extract second column of file to get read counts only 
	echo counts for `$sname` being extracted
	cut -f2 $x > `$sname`.tmp.counts
	# save shorthand sample names into an array  
	sname2=`$sname`
	myarray+=($sname2) 
done < <(ls -1 *.htseq-counts | sort) 
```

Paste all gene IDs into a file with each to make the gene expression matrix
```bash 
paste gene_IDs.txt *.tmp.counts > tmp_all_counts.txt
head tmp_all_counts.txt 
```

Save sample names in the array into text file 
```bash 
# look at the contents of the array we made with shorthand sample names 
echo ${myarray[@]}

# print contents of array into text file with each element on a new line 
printf "%s\n" "${myarray[@]}" > names.txt
cat names.txt
```

Put sample names in the file with counts to form row headers and complete the gene expression matrix
```bash 
cat <(cat names.txt | sort | paste -s) tmp_all_counts.txt > all_counts.txt
head all_counts
``` 

Remove all the tmp files 
```bash 
rm -f *tmp*
```

<br>


### Homework: 
- Catch-up if needed
- Exercises (if time)
- Bring questions of discussion points for end of day 2






### Exercises? (if time, could be take home) 
Just ideas...

1. All of your samples show good quality alignment %s except for one, which has an alignment rate of 55%?  
You noted when doing QC of the raw FASTQ files that this sample has many duplicated sequences.  
What might these sequences represent and how could you quickly test this?





