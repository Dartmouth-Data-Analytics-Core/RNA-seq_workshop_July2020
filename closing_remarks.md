# Closing remarks

### Workshop goals: 
- Develop a working understanding of the analytical workflow for a typical RNA-seq experiment
- Pre-process raw data in FASTQ format to generate a gene expression matrix
- Learn how to perform a detailed quality control analysis
- Develop a working understanding of the fundamental statistics behind a typical differential expression analysis using R/Bioconductor packages 
- Perform a differential expression analysis using R/Bioconductor packages 
- Learn how to explore the results and make robust insights from your data

### Day 1 analysis overview
![](figures/analysis_overview.png)

### Day 2 overview
![](figures/day2_summary.png)

### Some final take-aways from the workshop:
- Spend the time to plan, consult, practice, (and money) to generate high quality data that will provide robust inferences 
- If you are going to do a lot of Bioinformatics, you should get **really** good at the command-line (Bash), otherwise, pre-processing will be slow & painful
- Downstream differential expression analysis of a raw count matrix is best done in R, and requires some basic knowledge of statistics 
- Identify, understand, and check key QC metrics in the pre-processing and DE analysis portions, to ensure the quality of your results
- **PLEASE** correct for multiple testing!

### How to consolidate your learning: 
- Re-run the code a week or two after the workshop, as this is a great way to consolidate what you have learned at the command-line
- Edit the code, run sub-sections, read the `man` pages for commands, etc. to build a solid understanding of how everything works
- Practice pre-processing (Day-1 material) with the complete dataset (all chromosomes), that is available to you for approx. 1 month on discovery in `/dartfs-hpc/scratch/rnaseq1/data/`. This will give you experience running data from the entire genome, and an appreciation for the computational resources and required time to complete these tasks. 
- Read the methods sections of published papers that perform RNA-seq, to gain an appreciation for the range of approaches used in practice and how they are implemented 
- Read reviews like (this one)[https://pubmed.ncbi.nlm.nih.gov/31341269/] from Stark *et al*, 2019, *Nat. Rev. Genetics*, `RNA Sequencing: The Teenage Years`. 
- Ask us questions! (Bioinformatics office hours: https://dartmouth.zoom.us/s/96998379866, fridays at 1-2 pm. )

### I have my differentially expressed genes, what do I do now? 

**It depends!** What you do with your DEG results is dependent on your scientific question and hypothesis being tested. Two common downstream applications to use your DEGs for include: 

**1. Integrative genomics**  
You may have collected other types of genomics data (e.g. ChIP-seq, ATAC-seq) for your samples, or a public dataset exists from samples that are appropriate to integrate your RNA-seq data with. Once you have a confident set of DEGs, you can take an integrative genomics approach to answer questions that may not be able to be addressed with either type if data alone. 

For example, if you collected ChIP-seq for a specific transcription factor (TF) with paired RNA-seq data, you may wish to use your significant DEGs to identify genes whose expression is turned on/off by this TF under some treatment condition. 

**2. Gene ontology (GO) & pathway analyses**  
Unless you have very few significant DEGs, it may be difficult to identify complex patterns of functional gene regulation from looking at your list of DEGs. GO and pathway analysis methods are a very diverse collection of approaches that use statstical methodology to identify sets of genes, grouped using some sort of functional connection to each other (e.g. same biological pathway) that are *enriched* under one of the experimental conditions. 

Such methods can provide a lot of insight into the biological and molecular process controlling a phenptype, however there is an enormous range of tools that were designed for a specific applications (predominantly microarray data, **NOT** RNA-seq) and are not all appropriate for all data types. Selecting the appropriate tool then using it correctly are non-trivial and are commonly applied incorrectly in the literature. We encourage you to read more about these methods if you plan to use one in your own analysis, and we plan to cover this topic in future workshops. 

Here is some suggested reading regarding gene ontology and pathway analysis approaches:  
- [Gene set analysis approaches for RNA-seq data: performance evaluation and application guideline. *Briefings in Bioinformatics.* 2016.](https://doi.org/10.1093/bib/bbv069)
- [Ten Years of Pathway Analysis: Current Approaches and Outstanding Challenges. *PLoS Computational Biology.* 2012.](https://doi.org/10.1371/journal.pcbi.1002375)
- [Gene set enrichment analysis: a knowledge-based approach for interpreting genome-wide expression profiles. *PNAS* 2005.](https://doi.org/10.1073/pnas.0506580102) (the original GSEA paper)
- [Gene Set Enrichment Analysis Made Simple. *Stat Methods Med Red.* 2009.](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3134237/)

### Feedback: 

We ask that you all complete the survey that has been sent out over email so that we can gauge what worked well and what we need to improve for our next workshop. If you have additional thoughts that were not addressed in the survey, please feel free to contact any one of us, or reach out to the DAC meail directly (*DataAnalyticsCore@groups.dartmouth.edu*). 

<img src="figures/logo.jpg" width="250" height="140" >

We hope to offer this workshop again, as well as workshops covering other types of genomic data analysis and bioinformatics. If you have suggestions for workshops you would like to see, please let us know! 

Please feel free to reach out to us with questions about concepts discussed in the workshop, or for a analysis consultations. Our **bioinformatics office hours** on **Fridays 1-2pm** are a great place to do this! (currently on zoom: https://dartmouth.zoom.us/s/96998379866, pword: *bioinfo*)

### Now.... Discussion/question time! 
