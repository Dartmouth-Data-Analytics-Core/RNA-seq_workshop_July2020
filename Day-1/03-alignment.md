# Part 3 - Read alignment  

### Learning objectives: 
- Understand the major principles behind read alignment for RNA-seq data 
- Learn how alignment data is stored in SAM/BAM format
- Learn how to perform basic operations on BAM files using `Samtools`
- Perform an alignment with `STAR`
- How to view alignments using the `Integrative genomics Viewer (IGV)` 

Make a new directory to work in: 
```bash 
mkdir alignment
cd alignment
```

## Principles of read alignment for RNA-seq
Aligning millions of reads to very large reference genomes (such as the human genome) is generally done by splitting the reads and reference into a catalog of shorter reads with unique sequnece structures (kmers). It is improtant when selecting an alignement program to ensure that it is appropriate for the dataset you are working with, for example STAR (Spliced Transcripts Alignment to a Reference) is used to align reads that have come from spliced transcripts. A single read from a spliced transcriptome might map across a splice junction, such that the left side of the read and the right side of the read map hundreds of base pairs apart. If your dataset is prokaryotic (non-splicosomal) this would not be the appropriate program for you to align your reads, we would suggest looking into bwa-mem or bowtie2. If you are in a hurry and not interested in obtaining read alignments and only need count data quasi-mapping with a tool like Salmon might be a good option. 

![Read alignment](../figures/read_alignment.png)

#### Concepts for read alignment 

**Read clipping**
Aligners are capable of 'clipping' reads from sequence ends if they do not improve the quality of an alignment that exists for the rest of the sequence.  

There are two type of clipping:  
- *Soft-clipping*: bases at 5' and 3' ends of the read will be kept in the read sequence in the BAM file, but are NOT part of the alignment
- *Hard-clipping*: bases at 5' and 3' ends of the read will be removed from the BAM file altogether and are NOT part of the alignment 

Such clipping is commonly used by aligners to get rid of sequence contamination, e.g. adapter sequences or polyA tails from mRNAs, so that it does not affect the alignment. At least for RNA-seq, this is why you do not necessairily need to be very aggressive in read trimming and pre-processing steps. 

Clipping can be very advantageous, but also can potentially cause some issues, read more [here](https://sequencing.qcfail.com/articles/soft-clipping-of-reads-may-add-potentially-unwanted-alignments-to-repetitive-regions/). 

**Splicing**
Several different aligners exist, and one particularly important feature of an aligner is whether or not it is *splie-aware*. Splice-aware aligners, such as `STAR` and `HISAT2` are able to map reads over splice junctions by spanning the intronic region. This is obviously an important characteristic for RNA-seq data. Furthermore, if you provide coordinates of splice-junctions during alignment to aligners like STAR, it can improve the mapping over spliced regions and improve detection of novel splice-functions. 

**Genome vs transcriptome mapping?**
While there are times when one may want to map to a transcriptome, there are issues with this approach.  
- If your annotated transcriptome is not complete, you may fail to map some reads simply because the sequences aren't in your reference, which would not be as much of an issue if you mapped to the genome. 
- It is difficult to disambiguate which splice isoform a read should be aligned to in transcriptome mapping 
- Can't really idengtify novel transcripts this way

**What input do I need for an alignmnet?**
At miniumum:  
- `FASTQ` file(s)
- A reference genome (`.fasta`)

Optional: 
- `.gtf` file for the reference genome that species the geneomic feature annotation. As mentioned below, if you know where the splice-junctions in your genome are, you can give this to aligners such as STAR and they will use this information to improve the quality of mapping in these regions. 

**Alignment file formats**

Read alignments are stored in the SAM (.sam) and BAM (.bam) file format. SAM stands for *Sequence Alignment/Map* format and is in tab-delimited text format, making it a human readable file (should you dare to look inside). Bam files are the compressed, indexed, binary version of SAM files and are NOT human readable, but are much fatsre to parse and do complex operations on. You can read all about the SAM/BAM file format specification in the documentation [here](https://samtools.github.io/hts-specs/SAMv1.pdf). While you may never need to actually look inside of a SAM/BAM file, its important to have an understanding of what information is stored in one. 

Both formats contain a number of slots for each read alignment that describe key information about the alignment. 11 slots are mandatory, while others are optional and depend on the aligner used, and the settings used in that alignment.

![SAM file](../figures/sam-file.png)
The image for the example BAM file is take from the [SAM/BAM file format documentation](https://samtools.github.io/hts-specs/SAMv1.pdf)

Notes on select fields:

**FLAG** encodes important information about the read, for example, is it a primary, secondary, or supplementary alignment. Since one read will likely have a number of properties that we want to 'flag', SAM files use a special way of encoding the FLAG field to pack as much information as possible into one number. While we won't go into detail on this here, it works by decomposing large numbers into their constituents. I encourage you to go read more about FLAGs and how they are specified.  

Try running this to get basic information on FLAGs from samtools. 
```bash 
samtools flags
```

**MAPQ** corresponds to the quality of the mapping. These are calculated in the same way as the Phred scores `Q = -10 x log10(P)`, although are generally considers best guesses form the aligner. A MAPQ of 255 is used where mapping quality is not available. Some aligners also use specific values to represent certain types of alignments, which may affect use of downstream tools, so it is worth understanding those that are specific to your aligner. 

**CIGAR** is an alphanumerical string that tells you information about the alignment. For relatively short reads, these are nice, but for long reads, they are a headache. Numbers correspond to number of bases, and letters correspond to features of those bases.  

Letter key for CIGAR strings: 
M = match or mismatch  
S = soft clip  
H = hard clip  
I = insertion  
D = deletion  
N = skipping  

So for example, alignment in row 3 of our SAM file example above (`5S6M`) would describe an alignment where 5 bases are soft-clipped, followed by 6 matching bases. 

#### STAR (Spliced Transcripts Alignment to a Reference)
STAR is a very flexible, efficient, and quick read aligner. It uses a method of seed searching, clustering, stitching, and scoring to find the most probable match in the reference sequence for each read. A seed is the longest possible match between a read and the reference sequence. By using multiple seeds on a single read, reads are able to span hundreds of base pairs across splice junctions. Once a read is mapped into multiple seeds STAR attempts to map the remaining unmapped portions of the read by extending the seed match allowing for indels and mismatches. Any portion of the read that cannot be mapped is assumed to be contamination, leftover adapter sequences, or an incorrect base call and these bases are clipped (called soft-clipping). I encourage you to go look through the [user manual](https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf) if you plan to use STAR. 

#### Constructing a genome index 

Before running an alignment with STAR, you need to create an index of your reference genome, and specify the location of this index when you run the aligner. The index is in principle similar to how one might index a book, so that specific items or information can be found more quickly. For the genome index, we are indexing the genome so that the aligner can narrow down where a read may map to, help us to quickly index the genome and speed up mapping. 

For the purposes of this workshop and saving time, we have pre-built created a small genome index consisting of only chromosome 20, and will be using a subset of the total reads sequenced for sample `SRR1039508` which are known to align to chromosome 20. In practice you would use the entire genome to generate the index. This step is time consuming so we won't run it now, but the command used to create the index we will use is: 
```bash 
STAR --runThreadN 16 \
--runMode genomeGenerate \
--genomeDir hg38_chr20_index \
--genomeFastaFiles Homo_sapiens.GRCh38.dna.primary_assembly.chr20.fa \
--sjdbGTFfile Homo_sapiens.GRCh38.97.chr20.gtf \
--genomeSAindexNbases 11
```

Option details: 
- `--runThreadN`: no. of core/threads you want to use
- `--runMode`: the mode you want to run STAR in (for index generation, this should be genomeGenerate)
- `--genomeDir`: directory you want your genome to go to
- `--genomeFastaFiles`: path to genome .fasta 
- `--sjdbGTFfile`: path to genome annotation in .gtf format
- `--sjdbOverhang`: default is 100, usually set to the readlength -1

You can find the pre-built index at `/dartfs-hpc/rc/lab/`. 

Once you have generated an index, it is best not to do anything with it, except tell STAR where it is when you want to align reads. 

#### Aligning the reads

We are ready to align our reads to the genome, using the subset of reads sequenced for sample `SRR1039508` that are known to align to chromosome 20 (files: `SRR1039508_1.trim.chr20.fastq.gz` and `SRR1039508_2.trim.chr20.fastq.gz`). 

```bash
STAR --genomeDir hg38_chr20_index \
--readFilesIn SRR1039508_1.trim.chr20.fastq.gz SRR1039508_2.trim.chr20.fastq.gz \
--readFilesCommand zcat \
--sjdbGTFfile Homo_sapiens.GRCh38.97.chr20.gtf \
--runThreadN 10 \
--outSAMtype SAM \
--outFilterType BySJout \
--outFileNamePrefix SRR1039508.
```

> *NOTE:* I usually set `outSAMtype` to `BAM SortedByCoordinate`, so that I do not need to convert the defaulkt SAM file output by STAR to BAM, then sort it. However, since we want to look inside the file at the alignments, we are creatuing a SAM first, and will convert to a BAM afterwards. 

Option details: 
- `--genomeDir`: the path to the directory with genome indices
- `--sjdbGTFfile`: the path to the annotation file that includes cooordinates of splice-junctions
- `sjdbOverhang`: length of the sequence around the splice junction to be used in constructing the splice junctions database
- `--runThreadN`: number of threads to use in the run
- `--outFilterType`: how mapped reads will be filtered (normal/BySJout)
- `--outSAMtype`: (BAM unsorted/ BAM SortedByCoordinate)
- `--readFilesIn`: read files to map to reference alignment
- `--readFilesCommand`: uncompression command to apply to read files
- `--outFileNamePrefix`: prefix for outfiles generated in the run

Now, wait...

There are a number of other options you may wish to specify, dependning on your application and downstream analysis. These are the barebones options suggested for RNA-seq and optimized for mammalian genomes. Again, I encourage you to go look through the [user manual](https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf) if you plan to use STAR. 

**Aligment output**

Once the alignment has finished, you should have a number of new files in your directory. These are composed of:  
- `.bam` - your alignment file
- `Log.out` - the log of the STAR run 
- `Log.final.out` - the summary filke of the mammping statsistics
- `Log.progress.out` - a summary file that is updated with key mapping statistics as the run progresses
- `SJ.out.tab` - high-confidence splice-functions 

It is very important to do a detailed quality control review of your alignments across all your samples to identify any potential issues. We are going to (later) build a detailed multi-sample QC report using an independent tool, however we can have a quick look at the `Log.final.out` file from STAR as this contains the major mapping stats that we want to look at to assess how well reads were mapped to the reference. 

```bash
cat SRR1039508.Log.final.out
``` 

**Working with SAM/BAM files**  
We can also have a look around our SAM file to get to know it a bit better. Several tools exist that enable you to perform operations on SAM/BAM files. `Samtools` is perhaps the most widely used of these, and is used widely. You can find the documentation [here](http://www.htslib.org/doc/samtools.html). 

Lets use `Samtools` to have a look at the header for our SAM file
```bash 
samtools view -H SRR1039508.Aligned.out.sam  | head
```

Lets also have look at the first few alignment records in our BAM file
```bash 
samtools view SRR1039508.Aligned.out.sam | head
```

It is common to sort SAM/BAM files as this is required by many downstream tools that take alignment files as input. 
```bash
samtools sort SRR1039508.Aligned.out.sam -o SRR1039508.Aligned.out.sorted.sam
``` 

In practice, we can ask programs like STAR to give us indexed and sorted BAM files as output from the alignment, however this is not the case with all aligners. Now that we've looked at the alignments, we should convert our SAM to BAM for indexing and downstream analysis. 
```bash
samtools view -S -b SRR1039508.Aligned.out.sorted.sam > SRR1039508.Aligned.out.sorted.bam
``` 

We should also index this BAM, which will create a file with the same name, but the suffix `.bai`. 
```bash
samtools index SRR1039508.Aligned.out.sorted.bam
``` 

Another useful thing we might want to do with our BAM file is to count how many alignments have specific FLAG types (unique alignments, secondary, unmapped, properly paired). 
```bash
samtools flagstat SRR1039508.Aligned.out.sorted.bam
``` 

We can even use the specific FLAGs in the BAM file to extract specific alignments. For example, you might want to produce BAM files where all of the reads mapping to the forward and reverse strands are in separate files: 
```bash
# use -F option in view to filter out reads with REV alignment flag (16), leaving only FWD alignments 
samtools view -F 16 SRR1039508.Aligned.out.sorted.bam -o SRR1039508.Aligned.out.sorted.FWD.bam
samtools view -c SRR1039508.Aligned.out.sorted.FWD.bam

# use -f option in view to kepp reads with REV alignment flag (16), leaving only REV reads 
samtools view -f 16 SRR1039508.Aligned.out.sorted.bam -o SRR1039508.Aligned.out.sorted.REV.bam
samtools view -c SRR1039508.Aligned.out.sorted.REV.bam
``` 

You might just want to go straight to counting how many reads have a particular SAM flag
```bash
# count how many reads are NOT a primary alignment (FLAG=256)
samtools view -c -F 256 SRR1039508.Aligned.out.sorted.REV.bam
```

## Viewing Alignments in IGV

The Integrative Genomics Viewer (IGV) from the Broad Institute is an extremely useful tool for visulazation of alignment files (as well as other genomic file formats). Viewing your alignments in this way can be used to explore your data, troubleshoot issues you are having downstream of alignment, and inspect coverage for and quality of reads in specific regions of interest (e.g. in variant calling). I strongly encourage you to download the IGV for your computer from their [website](http://software.broadinstitute.org/software/igv/) and play around with some BAM file to get familar with all its various features. 

Here, we will create a small subset of a BAM file, download it onto our local machines, and view it using the IGV web app (for speed). You can open the IGV web app in your browser [here](https://igv.org/app/). 

Lets go ahead and subset our BAM file for reads aligning only to chromosome 22. We also need to create an index. 
```bash
# subset for reads just on chr 22 (to make it smaller)
samtools view -b -@ 8 -o chr20.bam SRR1039508.1.Aligned.sortedByCoord.out.bam 20

# index your new bam file 
samtools index chr20.bam
``` 

Now download the files onto your local machine, so that you can load them into the IGV web app. 
```bash
# make a directory and go into it 
mkdir rnaseq_wrksp/
cd rnaseq_wrksp/

# download the file using secure copy (scp)
##### modify this for your discovery ID, AND the directory your in 
scp d41294d@discovery7.dartmouth.edu:/dartfs-hpc/rc/lab/G/GMBSR_bioinfo/workshops/rnaseq-july20/data/bam/chr20.bam* .
##### you will be promoted for your password for discovery/polaris 

# you may also need to modify the permissions 
chmod a+rwx chr20*
```

Now navigate to the IGV web app, and follow the below steps:  
1. Under *'Genomes'* in the top right, select *'GRCh38/hg38'*
2. Under *'Tracks'* select *'Local file'*, then highlight both the .bam and .bai 
3. Navigate to chromosome 20, and zoom in to see genes and individual alignments
4. Use the setting at the side of the track to set colors by *read strand*
5. In the search bar, type in gene `SAMHD1`

![IGV](../figures/igv_webapp.png)

- What do you notice about the orientation of the aligning reads? 
- Do you think this gene is expressed in this sample? What about relative to nearby genes?
- Is there potentially going to be any ambiguity in read quantification for `SAMHD1`, given that our library was generared using a **stranded** protocol? 
