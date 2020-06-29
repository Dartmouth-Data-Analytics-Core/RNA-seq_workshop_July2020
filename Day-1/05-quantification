
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

# Run htseq-count on your .bam file 
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

![](../figures/ge-matrix.png)

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
