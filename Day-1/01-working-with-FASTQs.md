# Part 1 - Working with FASTQ files

### Learning objectives: 
- Understand the FASTQ file format and the formatting sequence information it stores
- Learn how to perform basic operations on FASTQ files in the command-line 

## FASTQ file format

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

## Working with FASTQ files 

### Basic operations 

While you don't normally need to go looking within an individual FASTQ file, it is very important to be able to manipulate FASTQ files in you are going to be doing any more involved bioinformatics. There are a lot of operations we can do with a FASTQ file to gain more information about our experiment, and being able to interact with FASTQ files can be useful for troubleshooting problems that might come up in your analyses. 

Due to their large size, we often perform gzip copmpression of FASTQ files so that they take up less space, however this means we have to unzip them if we want to look inside them and perform operations on them. We can do this with the `zcat` command. 

Lets use `zcat` and `head` to have a look at the first few records in our FASTQ file. 
```bash
zcat SRR1039508_1.trim.chr20.fastq.gz | head
zcat SRR1039508_2.trim.chr20.fastq.gz | head
```

How many records do we have in total? (don't forget to divide by 4..) 
```bash
zcat SRR1039508_1.trim.chr20.fastq.gz | wc -l
zcat SRR1039508_2.trim.chr20.fastq.gz | wc -l
```

What if we want to count how many unique barcodes exist in the FASTQ file. 
```bash
zcat SRR1039508_1.trim.chr20.fastq.gz | sed -n '2~4p' | head -10000 | grep -o .
```
Using sed with the -n option and '2~4p' will return the 2nd line, skip to 4 lines, and print again, recursively. We can use head to view this for the first 10,000 lines. 

Next, we can use grep to find all of the instances of each nucleotide across these lines. 
```bash
zcat SRR1039508_1.trim.chr20.fastq.gz | sed -n '2~4p' | head -10000 | grep -o .
```

Finally we can sort and count up the unique nucleotides that were found..
```bash
zcat SRR1039508_1.trim.chr20.fastq.gz | sed -n '2~4p' | head -10000 | grep -o . | sort | uniq -c
```
Now we have the number of each nuleotide across the reads from the first 1000 records. A quick and easy program to get GC content! 

### Searching for regular expressions

What if we wanted to find matches for a specific sequence, maybe after a start codon, in the FASTQ file 
```bash
zcat SRR1039508_1.trim.chr20.fastq.gz | grep -o "ATGGGATCA" | sort | uniq -c
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
