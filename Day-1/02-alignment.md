# Part 2 - Read alignment  

### Learning objectives: 
- Understand the major principles behind read alignment for RNA-seq data 
- Learn how alignment data is stored in SAM/BAM format
- Learn how to perform basic operations on BAM files using `Samtools`
- Perform an alignment with `STAR`
- How to view alignments using the `Integrative genomics Viewer (IGV)` 

Make a new directory to work in: 
```bash 
# go to our home dir for the wrkshop
rnaw

# make the directroy and cd into it 
mkdir results/alignment | cd
cd results/alignment
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
Several different aligners exist, and one particularly important feature of an aligner is whether or not it is *splie-aware*. Splice-aware aligners, such as `STAR` and `HISAT2` are able to map reads over splice junctions by spanning the intronic region. This is obviously an important characteristic for RNA-seq data. Furthermore, if you provide coordinates of splice-junctions to aligners like STAR, it can improve the mapping over spliced regions and improve detection of novel splice-functions. 

**Genome vs transcriptome mapping?**
While there are times when one may want to map to a transcriptome, there are issues with this approach.  
- If your annotated transcriptome is not complete, you may fail to map some reads simply because the sequences aren't in your reference, which would not an issue if you mapped to the genome. 
- With multiple splice isoforms it is difficult to disambiguate which splice isoform a read should be aligned to in transcriptome mapping. 
- You cannot identify novel transcripts this way.

**What input do I need for an alignmnet?**
At miniumum:  
- `FASTQ` file(s)
- A reference genome (`.fasta`)

Optional: 
- `.gtf` file for the reference genome that species the genomic feature annotation. As mentioned above, if you know where the splice-junctions in your genome are, you can give this to aligners such as STAR and they will use this information to improve the quality of mapping in these regions. 

**Alignment file formats**

Read alignments are stored in the SAM (.sam) and BAM (.bam) file format. SAM stands for *Sequence Alignment/Map* format and is in tab-delimited text format, making it a human readable file (should you dare to look inside, these file are huge). Bam files are the compressed, indexed, binary version of SAM files and are NOT human readable, but are much faster to parse and do complex downstream operations on. You can read all about the SAM/BAM file format specification in the documentation [here](https://samtools.github.io/hts-specs/SAMv1.pdf). While you may never need to actually look inside of a SAM/BAM file, it is important to have an understanding of what information is stored in one. 

Both formats contain a number of slots for each read alignment that describe key information about the alignment. 11 slots are mandatory, while others are optional and depend on the aligner used, and the settings used in that alignment.

![SAM file](../figures/sam-file.png)
The image for the example BAM file is take from the [SAM/BAM file format documentation](https://samtools.github.io/hts-specs/SAMv1.pdf)

Notes on select fields:

**FLAG** encodes important information about the read, for example, is it a primary, secondary, or supplementary alignment. Since a single read will likely have a number of properties that we want to 'flag', SAM files use a special way of encoding the FLAG field to pack as much information as possible into a single number. While we won't go into detail on this here, it works by decomposing large numbers into their constituents. I encourage you to go read more about FLAGs and how they are specified.  

This command will provide basic information on FLAGs from samtools. 
```bash 
samtools flags
```

**MAPQ** corresponds to the quality of the mapping. These are calculated in the same way as the Phred scores `Q = -10 x log10(P)`, although are generally considered to be a best guess form the aligner. A MAPQ of 255 is used where mapping quality is not available. Some aligners also use specific values to represent certain types of alignments, which may affect use of downstream tools, so it is worth understanding those that are specific to your aligner. 

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
STAR is a very flexible, efficient, and quick read aligner. It uses a method of seed searching, clustering, stitching, and scoring to find the most probable match in the reference sequence for each read. A seed is the longest possible match between a read and the reference sequence. By using multiple seeds on a single read, reads are able to span hundreds of base pairs across splice junctions. Once a read is mapped into multiple seeds STAR attempts to map the remaining unmapped portions of the read by extending the seed match allowing for indels and mismatches. Any portion of the read that cannot be mapped is assumed to be contamination, leftover adapter sequences, or an incorrect base call and these bases are clipped (soft-clipping). I encourage you to go look through the [user manual](https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf) if you plan to use STAR. 

#### Constructing a genome index 

Before running an alignment with STAR, you need to create an index of your reference genome, and specify the location of this index when you run the aligner. The index is in principle similar to how one might index a book, so that specific items or information can be found more quickly. For the genome index, we are indexing the genome so that the aligner can narrow down where a read may map to and speed up mapping. 

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

You can find the pre-built index at `/dartfs-hpc/scratch/rnaseq1/refs/hg38_chr20_index/`. 

Once you have generated an index, it is best not to do anything with it, except tell STAR where it is when you want to align reads. 

#### Aligning the reads

We are ready to align our reads to the genome, using the subset of reads sequenced for sample `SRR1039508` that are known to align to chromosome 20 (files: `SRR1039508_1.trim.chr20.fastq.gz` and `SRR1039508_2.trim.chr20.fastq.gz`). 

```bash
STAR --genomeDir /dartfs-hpc/scratch/rnaseq1/refs/hg38_chr20_index \
  --readFilesIn ../trim/SRR1039508_1.trim.chr20.fastq.gz ../trim/SRR1039508_2.trim.chr20.fastq.gz \
  --readFilesCommand zcat \
  --sjdbGTFfile /dartfs-hpc/scratch/rnaseq1/refs/Homo_sapiens.GRCh38.97.chr20.gtf \
  --runThreadN 4 \
  --outSAMtype SAM \
  --outFilterType BySJout \
  --outFileNamePrefix SRR1039508.
```

> *NOTE:* I usually set `outSAMtype` to `BAM SortedByCoordinate`, so that I do not need to convert the default SAM file output by STAR to BAM, then sort it. However, since we want to look inside the file at the alignments, we are creating a SAM first, and will convert to a BAM afterwards. 

Option details: 
- `--genomeDir`: the path to the directory with genome indices
- `--sjdbGTFfile`: the path to the annotation file that includes cooordinates of splice-junctions
- `sjdbOverhang`: length of the sequence around the splice junction to be used in constructing the splice junctions database
- `--runThreadN`: number of threads to use in the run
- `--outFilterType`: how mapped reads will be filtered (normal/BySJout)
- `--outSAMtype`: (SAM/BAM unsorted/ BAM SortedByCoordinate)
- `--readFilesIn`: read files to map to reference alignment
- `--readFilesCommand`: uncompression command to apply to read files
- `--outFileNamePrefix`: prefix for outfiles generated in the run

Now, wait...

There are a number of other options you may wish to specify, depending on your application and downstream analysis. These are the barebones options suggested for RNA-seq and optimized for mammalian genomes. Again, I encourage you to go look through the [user manual](https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf) if you plan to use STAR. 

**Aligment output**

Once the alignment has finished, you should have a number of new files in your directory. These are composed of:  
- `.sam` - your alignment file
- `Log.out` - the log of the STAR run 
- `Log.final.out` - the summary file of the mapping statistics
- `Log.progress.out` - a summary file that is updated with key mapping statistics as the run progresses
- `SJ.out.tab` - high-confidence splice-functions 

It is very important to do a detailed quality control review of your alignments across all your samples to identify any potential issues. We are going to (later) build a detailed multi-sample QC report using an independent tool, however we can have a quick look at the `Log.final.out` file from STAR as this contains the major mapping stats that we want to look at to assess how well reads were mapped to the reference. 

```bash
cat SRR1039508.Log.final.out
``` 

**Working with SAM/BAM files**  
We can also have a look around our SAM file to get to know it a bit better. Several tools exist that enable you to perform operations on SAM/BAM files. `Samtools` is perhaps the most widely used of these. You can find the documentation [here](http://www.htslib.org/doc/samtools.html). 

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

In practice, we can ask programs like STAR to give us indexed and sorted BAM files as output from the alignment, however this is not the case with all aligners and in these cases you will have to sort and index files after the alignment is complete. Now that we've looked at the alignments, we should convert our SAM to BAM for indexing and downstream analysis. 
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

The Integrative Genomics Viewer (IGV) from the Broad Institute is an extremely useful tool for visualization of alignment files (as well as other genomic file formats). Viewing your alignments in this way can be used to explore your data, troubleshoot issues you are having downstream of alignment, and inspect coverage and quality of reads in specific regions of interest (e.g. in variant calling). I strongly encourage you to download the IGV for your computer from their [website](http://software.broadinstitute.org/software/igv/) and play around with a BAM file to get familar with all its various features. 

Here, we will create a small subset of a BAM file, download it onto our local machines, and view it using the IGV web app (for speed). You can open the IGV web app in your browser [here](https://igv.org/app/). 

Lets go ahead and subset our BAM file for reads aligning only to chromosome 20. We also need to create an index. 
```bash
# subset for reads just on chr 20 (to make it smaller)
samtools view -b -@ 8 -o chr20.bam SRR1039508.Aligned.out.sorted.bam 20

# index your new bam file 
samtools index chr20.bam
``` 

Now download the files onto your local machine, so that you can load them into the IGV web app. 
```bash
# make a directory and go into it (ON YOUR LOCAL MACHINE)
mkdir rnaseq_wrksp/
cd rnaseq_wrksp/

# download the file using secure copy (scp)
##### modify this for your discovery ID, AND the directory your in 
scp d41294d@discovery7.dartmouth.edu:/dartfs-hpc/scratch/omw/results/alignment/chr20.bam* .
##### you will be promoted for your password for discovery 

# you may also need to modify the permissions 
chmod a+rwx chr20*
```

Now navigate to the IGV web app, and follow the below steps:  
1. Under *'Genomes'* in the top right, select *'GRCh38/hg38'*
2. Under *'Tracks'* select *'Local file'*, then highlight both the .bam and .bai 
3. Navigate to chromosome 20, and zoom in to see genes and individual alignments
4. Use the setting at the right side of the track to set colors by *read strand*
5. In the search bar, type in gene `SAMHD1`

![IGV](../figures/igv_webapp.png)

- What do you notice about the orientation of the aligning reads? 
- Do you think this gene is expressed in this sample? What about relative to nearby genes?
- Is there potentially going to be any ambiguity in read quantification for `SAMHD1`, given that our library was generared using a **stranded** protocol? 

## Run the alignment for all samples
```bash
ls ../trim/*_1.trim.chr20.fastq.gz | while read x; do

  # save the file name
  sample=`echo "$x"`
  # get everything in file name before "/" (to remove '../trim/')
  sample=`echo "$sample" | cut -d"/" -f3`
  # get everything in file name before "_" e.g. "SRR1039508"
  sample=`echo "$sample" | cut -d"_" -f1`
  echo processing "$sample"
  
  # run STAR for each sample 
  STAR --genomeDir /dartfs-hpc/scratch/rnaseq1/refs/hg38_chr20_index \
    --readFilesIn ../trim/${sample}_1.trim.chr20.fastq.gz ../trim/${sample}_2.trim.chr20.fastq.gz \
    --readFilesCommand zcat \
    --sjdbGTFfile /dartfs-hpc/scratch/rnaseq1/refs/Homo_sapiens.GRCh38.97.chr20.gtf \
    --runThreadN 4 \
    --outSAMtype BAM SortedByCoordinate \
    --outFilterType BySJout \
    --outFileNamePrefix ${sample}.
done

# index the BAMs 
samtools index *sorted.bam
```
Note that I change `--outSAMtype` to `BAM sortedByCoord` so that we dont have to convert SAM to BAM and run `sort`. 

View the reports quickly: 
```bash 
cat *Log.final.out
```

Now we can move on to do a comprehensive QC analysis of our alignments. 

