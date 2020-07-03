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
- **Row 1:** Header line that stores information about the read (always starts with an `@`), such as the *instrument ID*, *flowcell ID*, *lane on flowcell*, *file number*, *cluster coordinates*, *sample barcode*, etc.
- **Row 2:** The sequence of bases called
- **Row 3:** Usually just a `+` and sometimes followed by the read information in line 1
- **Row 4:** Individual base qualities (must be same length as line 2)

Quality scores, also known as **Phred scores**, in row 4 represent the probability that the associated base call is incorrect, which are defined by the below formula for current Illumina machines:
```
Q = -10 x log10(P), where Q = base quality, P = probability of incorrect base call
```
or 
```
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

While you don't normally need to go looking within an individual FASTQ file, it is very important to be able to manipulate FASTQ files if you are going to be doing more involved bioinformatics. There are a lot of operations we can do with a FASTQ file to gain more information about our experiment, and being able to interact with FASTQ files can be useful for troubleshooting problems that might come up in your analyses. 

Due to their large size, we often perform gzip copmpression of FASTQ files so that they take up less space, however this means we have to unzip them if we want to look inside them and perform operations on them. We can do this with the `zcat` command and a pipe (|). The pipe command is a way of linking commands, the pipe sends the output from the first command to the second command. `zcat` lists the contents of a zipped file to your screen, and head limits the output to the first ten lines. 

Lets use `zcat` and `head` to have a look at the first few records in our FASTQ file. 
```bash
zcat SRR1039508_1.chr20.fastq.gz | head
zcat SRR1039508_2.chr20.fastq.gz | head
```

How many records do we have in total? (don't forget to divide by 4..) 
```bash
zcat SRR1039508_1.chr20.fastq.gz | wc -l
zcat SRR1039508_2.chr20.fastq.gz | wc -l
```
Paired-end reads should have the same number of records! 

What if we want to count how many unique barcodes exist in the FASTQ file. To do this, we would need to print all the sequence lines of each FASTQ entry, then search those for the barcode by specifying a regular expression. To print all the sequence lines (2nd line) of each FASTQ entry, we can use a command called ***sed***, short for ***stream editor*** which allows you to streamline edits to text that are redirected to the command. You can find a tutorial on using **sed** [here](https://www.digitalocean.com/community/tutorials/the-basics-of-using-the-sed-stream-editor-to-manipulate-text-in-linux). 

First we can use sed with with the `'p'` argument to tell it that we want the output to be printed, and the `-n` option to tell sed we want to suppress automatic printing (so we don't get the results printed 2x). Piping this to `head` we can get the first line of the first 10 options in the FASTQ file (the header line). We specify `'1-4p'` as we want sed tp *print 1 line, then skip forward 4*. 
```bash
zcat SRR1039508_1.chr20.fastq.gz | sed -n '1~4p' | head -10
```

Using this same approach, we can print the second line for the first 10000 entires of the FASTQ file, and use the ***grep*** command to search for regular expressions in the output. Using the `-o` option for grep, we tell the command that we want it to print lines that match the character string. 
```bash
# print the first 10 lines to confirm we are getting bthe sequence lines 
zcat SRR1039508_1.chr20.fastq.gz | sed -n '2~4p' | head -10

# pipe the sequence line from the first 10000 FASTQ records to grep to search for our (pretend) adapter sequence
zcat SRR1039508_1.chr20.fastq.gz | sed -n '2~4p' | head -10000 | grep -o "ATGGGA"
```

This is a bit much to count by each, so lets count the how many lines were printed by grep using the ***wc*** (word count) command with the `-l` option specified for lines.
```bash
zcat SRR1039508_1.chr20.fastq.gz | sed -n '2~4p' | head -10000 | grep -o "ATGGGA" | wc -l
```

Using a similar approach, we could count up all of the instances of individual DNA bases (C,G) called by the sequencer in this sample. Here we use the ***sort*** command to sort the bases printed by grep, grep again to just get the bases we are interested in, then using the ***uniq*** command with the `-c` option to count up the unique elements. 
```bash
zcat SRR1039508_1.chr20.fastq.gz | sed -n '2~4p' | head -10000 | grep -o . | sort | grep 'C\|G' | uniq -c 
```
Now we have the number of each nuleotide across the reads from the first 10000 records. A quick and easy program to get GC content. GC content is used in basic quality control of sequence from FASTQs to check for potential contamination of the sequencing library. We just used this code to check 1 sample, but what if we want to know for our 4 samples?

## For & while loops 

Loops allow us repeat operations over a defined variable or set of files. Essentially, you need to tell Bash what you want to loop over, and what operation you want it to do to each item. 

Notice that the variable ***i*** set in the conditions for our loop is used to reference all the elements to be looped over in the operation using the term ***$i*** in this **for*** loop example: 
```bash 
# loop over numbers 1:10, printing them as we go
for i in {1..10}
do
   echo "$i"
done
``` 


Alternatively, if you do not know how many times you might need to run a loop, using a ***while*** loop may be useful, as it will continue the loop until the boolean (logical) specified in the first line evaluates to `false`. An example would be looping over all of the files in your directory to perform a specific task. e.g. 
```bash
ls *.fastq.gz | while read x; do 
echo $x is being processed...; zcat $x | head -n 4 
done
```

Perhaps this sequence represents some a contaminating sequence from the run that we want to quickly screen all of our samples for (e.g. from bacteria). We can do this by searching for matches and counting how many times it was found, and repeating this process for each sample using a for loop. 
```bash
ls *.fastq.gz | while read x; do 
echo $x;zcat $x | sed -n '2~4p' | head -4 | grep -o "ATGGGA" | wc -l;
done
```

We could use one of these loops to perform the nucleotide counting task that we performed on a single sample above. 
```bash
ls *.fastq.gz | while read x; do 
echo processing sample $x; zcat $x | sed -n '2~4p' | sed -n '1,10000p' | grep -o . | sort | uniq -c;
done
```

## Scripting in bash 

So loops are pretty useful, but what if we wanted to make it even simpler to run. Maybe we even want to share the program we just wrote with other lab members so that they can execute it on their own FASTQ files. One way to do this would be to write this series of commands into a Bash script, that can be executed at the command line, passing the files you would like to be operated on to the script. 

To generate the script (suffix `.sh`) we could use the `nano` editor: 
```bash 
nano count_GC_content.sh
```

Add our program to the script, using a shebang `#!/bin/bash` at the top of our script to let the shell know this is a bash script. As in the loops we use the `$` to specify the input variable to the script. `$1` represents the variable that we want to be used in the first argument of the script. Here, we only need to provide the file name, so we only have 1 `$`, but if we wanted to create more variables to expand the functionality of our script, we would do this using `$2`, `$3`, etc. 
```bash 
#!/bin/bash
echo processing sample "$1"; zcat $1 | sed -n '2~4p' | sed -n '1,10000p' | grep -o . | sort | grep 'C\|G' | uniq -c;
```

Now run the script, specifying the a FASTQ file as variable 1 (`$1`)
```bash
# have a quick look at our script 
cat count_GC_content.sh

# now run it with bash 
bash count_GC_content.sh SRR1039508_1.chr20.fastq.gz
```

Now we can use our while loop again to do this for all the FASTQs in our directory 
```bash
ls *.fastq.gz | while read x; do 
bash count_GC_content.sh $x
done
```

What if we wanted to write the output into a file instead of printing to the screen? We could save the output to a *Standard output* (stout) file that we can look at, save to review later, and document our findings. The `1>>` redirects the output that would print to the screen to a file.
```bash
# create the text file you want to write to
touch stout.txt

# run the loop 
ls *.fastq.gz | while read x; do 
bash count_GC_content.sh $x 1>> stout.txt
done

# view the file 
cat stout.txt
```

These example programs run fairly quickly, but stringing together mutiple commands in a bash script is common and these programs take much longer to run. In these cases we might want to close our computer and go and do some other stuff while our program is running. We can do this using `nohup` which allows us to run a series of commands in the background, but disconnects the process from the shell you initally submit it through, so you are free to close this shell and the process will continue to run until completion. e.g. 
```bash
nohup bash count_GC_content.sh SRR1039508_1.chr20.fastq.gz > result.txt &

# show the result 
cat nohup.out 
```

## Quality control of FASTQ files 

While the value of these sorts of tasks may not be immediately clear, you can imagine that if we wrote some nice programs like we did above, and grouped them together with other programs doing complimentary tasks, we would make a nice bioinformatics software package. Fortunately, people have already started doing this, and there are various collections of tools that perform specific tasks on FASTQ files. 

One excellent tool that is specifically designed assess quality of FASTQ file is [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/). FastQC is composed of a number of analysis modules that calculate various QC metrics from FASTQ files (such as GC content, distribution of base quality, etc.) and summarizes this all into a nice QC report in HTML format, that can be opened in a web browser. 

Lets have a look at some example QC reports from the FastQC documentation: 

[Good Illumina Data FastQC Report](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/good_sequence_short_fastqc.html)
[Bad Illumina Data FastQC Report](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/bad_sequence_fastqc.html)

Lets run FASTQC on our data and move the results to a new directory. 
```bash
# specify the -t option for 4 threads to make it run faster
fastqc -t 4 *.fastq.gz

# move results to a new folder 
mkdir fastqc_results
mv *fastqc* fastqc_results

# move into it and ls 
cd fastqc_results
ls -lah 
```

**Note**: FastQC does not use the entire dataset, just the first few thousand reads in the FASTQ file, therefore there could be some bias introduced by this, although we assume there isn't since entires are placed into FASTQ files randomly. 

Opening and evaluating an individual .html file for each FASTQ file is obviously going to be tedious and slow. Luckily, someone built a tool to speed this up. [MultiQC](https://multiqc.info/) *MultiQC* searches a specified directory (and subdirectories) for log files that it recognizes and synthesizes these into its own browsable, sharable, interactive .html report that can be opened in a web-browser. *MultiQC* recognizes files from a very wide range of bioinformatics tools (includeing FastQC), and allows us to compare QC metrics generated by various tools across samples and analyze our experiment as a whole. 

Lets run MultiQC on our FastQC files:
```bash 
multiqc .
```

Copy to your home directory and open in a web-broswer: 
```
cp *multiqc* $HOME
```

You can find the MultiQC report run on the complete dataset across all samples in the dataset in the github repository, under `QC-reports`. Lets open it and explore our QC data. 

**What do we think about the quality of our dataset?**

## Read pre-processing & trimming

An additional QC step we can perform before proceeding with the read alignments is to *pre-process* or *trim* the sequences in the FASTQ files to remove sequences that we are not interested in, or were not called confidently by the sequenecer. This step is **optional** in most analysis, although should be based on an empirical decision that leverages the QC assessment of raw FASTQs. For example, if we see we have a large number of adapter seqeunces in our data, or a high proportion of low-quality bases near our read ends, we may wish to trim our raw reads. In the absence of such issues, we could skip this step. 

Another reason this step is optional is that many aligners can account for mismatches or low quality bases at the end of reads during the alignment process. The aligner will *soft-clip* these bases from the alignment when it maps that read (more on this later) so we do not necessarily need to have these bases trimmed off all of the reads. 

### Principles of read trimming: Downstream steps are more efficient

Several algorithms exist for trimming reads in FASTQ format. Generally, these algorithms work by looking for matches to the sequence you specify at the 5' and 3' end of a read. You can specify the minimum number of bases you would like to be considered a match, as the algorithm will trim partial matches to the sequence you specify. Examples of sequences you might want to remove include: adapter sequences, polyA tails, or low quality bases. 

![Read alignment](../figures/read_processing.png)


### Read trimming with cutadapt 

[Cutadapt](https://cutadapt.readthedocs.io/en/stable/) is a useful tool for cleaning up sequencing reads, and allows for multiple adapters to be specified simulatenously, and has numerous options that can be tweaked to control its behaviour. 

Basic usage of cutadapt: 
```bash 
cutadapt -a ADAPTER -g ADAPT2 [options] -o output.fastq input.fastq.gz
``` 
- `-a` specifies an adapter to trim from the 3' end of read 1
- `g` specifies an adapter to trim from the 5' end of read 1
- `o` specifies name of out file

For paired-end reads:
```bash 
cutadapt -a ADAPT1 -g ADAPT2 [options] -o out1.fastq.gz -p out2.fastq input1.fastq.gz input2.fastq.gz
``` 
Capital letters are used to specify adapters for read 2. 

If we wanted to trim polyA sequences and save the output to a report called cutadapt.logout, we could use:  
```bash
cutadapt -a 'A{76}' -o out.trimmed.fastq.gz input.fastq.gz > cutadapt.logout;
```
`-a A{75}` tells cutadapt to search for streches of A bases at the end of reads, with a maximum length of the read length (75bp). 

Since the polyA and adapter sequence contamination is relatively low for this dataset, we won't trim any specific sequences, although we will perform basic quality and length processing of the raw reads. Lets make a new directory and do this for do this for one sample. 
```bash
mkdir ../results/trim 
cd ../results/trim

cutadapt \
-o SRR1039508_1.trim.chr20.fastq.gz \
-p SRR1039508_2.trim.chr20.fastq.gz \
../../raw_data/SRR1039508_1.chr20.fastq.gz ../../raw_data/SRR1039508_1.chr20.fastq.gz \
-m 1 -q 20 -j 4 > SRR1039508.cutadapt.report
```

- `-m` removes reads that are samller than the minimum threshold
- `-q` qulaity threshold for trimming bases
- `-j` number of cores/threads to use

Now lets run this on all of our samples:
```bash 
ls ../../raw_data/*_1.trim.chr20.fastq.gz | while read x; do

# save the file name
sample=`echo "$x"`
# get everything in file name after "/" and before "_" e.g. "SRR1039508"
sample=`echo "$sample" | cut -d"/" -f4|cut -d"_" -f1`
echo processing "$sample"

# run cutadapt for each sample 
cutadapt \
   -o ${sample}_1.trim.chr20.fastq.gz \
   -p ${sample}_2.trim.chr20.fastq.gz \
   ../../raw_data/${sample}_1.chr20.fastq.gz ../../raw_data/${sample}_2.chr20.fastq.gz \
   -m 1 -q 20 -j 4 > $sample.cutadapt.report;
done
```

You should now have trimmed FASTQ files in this directory that we will use for the alignment. You should also be able to see and print each of your reports from cutadapt. 
```bash
cat *cutadapt.report
```

**Additional note:** For data generated at Dartmouth, since much of the data in the Genomics core is generated using an **Illumina NextSeq 500**, we also often use the `--nextseq-trim` option in cutadapt. This works the qulaity threshold BUT ignores Q-scores for streches of G bases, as some Illumina instruments, such as the NextSeq, generate strings of Gs when when the sequencer 'falls off' the end of a fragment and dark cycles occur, and therefore provides more appropriate quality trimming for data generated on these instrucments. 
