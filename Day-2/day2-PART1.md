---
title: "Day 2, PART1"
author: "Prepared by the Data Analytics Core (Center for Quantitative Biology at Dartmouth)"
output:
    html_document:
      keep_md: TRUE
      theme: default
      number_sections: TRUE

---
*Bulk RNA-seq data analysis data analysis workshop, July 2020*

## Exploratory data analysis in R

### Introduction and set-up

Several popular R packges designed for exploration and statistical analysis of bulk RNA-seq data exist, including [*EdgeR*](https://www.bioconductor.org/packages/release/bioc/html/edgeR.html), [*limma-voom*](http://bioconductor.org/packages/release/bioc/html/limma.html), [*DESeq2*](https://bioconductor.org/packages/release/bioc/html/DESeq2.html). For the purposes of this workshop, we will use DESeq2 to perform the parts of the analysis, including reading in the data, normalization of read counts, and fitting statistical models to test differential expression. [Detailed tutorials](https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html) for using DESeq2 can be found on its Bioconductor page. 

DESeq2 is a very well organized package that applies robust algorithms to perform several aspects of RNA-seq data analysis. If you plan to use DESeq2 for your work, you should read both the tutorials made available on their Bioconductor page, and the original manuscript for DESeq2, in [Love *et al*, 2014](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-014-0550-8) to develop an understanding of the theory behind DESeq2 and the proccesses implemented by its functions. Despite DESeq2's extensive functionality, a different package may be appropriate for the analysis of RNA-seq data with unique experimental designs, for example using linear-mixed effects models to perform differential expression analysis of clustered datasets. 

NOTE: You must change the below line, and all other lines loading images, to the directory on your computer!!

<center>
![Overview](/Users/OwenW/Downloads/RNA-seq_workshop_July2020-master/figures/overview.png)
</center>

Set the root directory for the whole markdown. THIS MUST BE SET TO THE LOCATION OF THE FOLDER YOU DOWNLOADED!

```r
knitr::opts_knit$set(root.dir = '/Users/OwenW/Downloads/RNA-seq_workshop_July2020-master/')
```

Lets start by loading the packages we will need:

```r
library(dplyr)
library(ggplot2)
library(tximport)
library(DESeq2)
library(biomaRt)
library(vsn)
library(pheatmap)
library(gplots)
library(RColorBrewer)
library(ComplexHeatmap)
library(circlize)
library(xtable)
library(kableExtra)
```

### Read in the raw counts

The dataset that we are using comes from [this paper](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0099625). This data was collected from human airway smooth muscle cells to test gene pathways effected by exposure to Glucocorticoid, which have been historically used for their anti-inflammatory effects to treat asthmatics. Four cell lines were treated with either a control vehicle (untreated), dexamethasone (dex), albuterol (alb), or both dexamethasone and albuterol (co-treated) for 18 hours before transcriptomes were extracted. 

### Read in raw count data 

Now we can read in our data. How you read your data into DESeq2 depends on what format your raw reads counts are in (individual files for each sample, or a gene expression matrix) and how your read counts were quantified (e.g. at the gene or transcript level). DESeq2 provides a specific function (`DESeqDataSetFromHTSeqCount`) to read in gene-level read count abundances from *htseq-count*. 

```r
# read in the matrix we generated using htseq-count 
cts <- as.matrix(read.table("Day-2/all_counts.txt", 
                            sep="\t", header = TRUE, row.names=1, 
                            stringsAsFactors = F))

# quick look at the matrix 
head(cts)
```

```
##                 SRR1039508 SRR1039509 SRR1039510 SRR1039511 SRR1039512
## ENSG00000000003        698        468        758        471        893
## ENSG00000000005          0          0          0          0          0
## ENSG00000000419        464        507        455        520        608
## ENSG00000000457        258        206        237        227        262
## ENSG00000000460         58         55         88         57         39
## ENSG00000000938          0          0          0          0          2
##                 SRR1039513 SRR1039514 SRR1039515 SRR1039516 SRR1039517
## ENSG00000000003        420        985        878       1180       1074
## ENSG00000000005          0          1          0          0          0
## ENSG00000000419        357        869        447        580        781
## ENSG00000000457        163        328        217        242        330
## ENSG00000000460         35         87         47         77         64
## ENSG00000000938          0          0          1          1          0
##                 SRR1039518 SRR1039519 SRR1039520 SRR1039521 SRR1039522
## ENSG00000000003       1114       1141        802        595       1125
## ENSG00000000005          1          0          0          0          0
## ENSG00000000419        717        510        414        500        535
## ENSG00000000457        294        199        230        229        248
## ENSG00000000460         66         67         76         60         84
## ENSG00000000938          1          0          0          0          0
##                 SRR1039523
## ENSG00000000003        872
## ENSG00000000005          0
## ENSG00000000419        673
## ENSG00000000457        277
## ENSG00000000460         91
## ENSG00000000938          2
```

```r
tail(cts)
```

```
##                        SRR1039508 SRR1039509 SRR1039510 SRR1039511 SRR1039512
## ENSG00000288111                 0          0          0          0          0
## __no_feature               936867     842433    1027258     964295    1203807
## __ambiguous               1729002    1587216    1724395    1639595    2089004
## __too_low_aQual             31467      34647      38352      42305      51152
## __not_aligned              547100     741136     670279     420256     476785
## __alignment_not_unique    1033615     901230    1011957     936332    1156255
##                        SRR1039513 SRR1039514 SRR1039515 SRR1039516 SRR1039517
## ENSG00000288111                 0          0          0          0          0
## __no_feature               618312    1625007    1108750    1218372    1467646
## __ambiguous               1230442    2787173    1708909    2062034    2535139
## __too_low_aQual             33739      81543      43103      52636      68199
## __not_aligned              344865    1012304     474310     557462     576573
## __alignment_not_unique     660376    1591279     949990    1168480    1370005
##                        SRR1039518 SRR1039519 SRR1039520 SRR1039521 SRR1039522
## ENSG00000288111                 1          0          0          0          0
## __no_feature              1257366     954313     936202     888824    1284938
## __ambiguous               2262879    1621650    1540090    1673566    2037477
## __too_low_aQual             61530      42754      44372      50221      68854
## __not_aligned              557606     429632     413797     505024     589176
## __alignment_not_unique    1311827     951332     866864     954938    1128194
##                        SRR1039523
## ENSG00000288111                 0
## __no_feature              1393701
## __ambiguous               2245960
## __too_low_aQual             57949
## __not_aligned              537783
## __alignment_not_unique    1279829
```

```r
# filter out these last 5 rows 
cts <- cts[1:(nrow(cts)-5),]
tail(cts)
```

```
##                 SRR1039508 SRR1039509 SRR1039510 SRR1039511 SRR1039512
## ENSG00000288106          2          2          2          7         12
## ENSG00000288107          1          1          0          0          0
## ENSG00000288108          0          0          0          0          0
## ENSG00000288109          0          0          0          0          0
## ENSG00000288110          0          0          0          0          0
## ENSG00000288111          0          0          0          0          0
##                 SRR1039513 SRR1039514 SRR1039515 SRR1039516 SRR1039517
## ENSG00000288106          9         13         14         12          7
## ENSG00000288107          0          2          1          0          2
## ENSG00000288108          0          0          0          0          0
## ENSG00000288109          1          0          1          0          0
## ENSG00000288110          0          0          0          0          0
## ENSG00000288111          0          0          0          0          0
##                 SRR1039518 SRR1039519 SRR1039520 SRR1039521 SRR1039522
## ENSG00000288106          9          5          7          3          1
## ENSG00000288107          0          0          0          1          2
## ENSG00000288108          0          0          0          0          0
## ENSG00000288109          0          0          0          0          0
## ENSG00000288110          0          0          0          0          0
## ENSG00000288111          1          0          0          0          0
##                 SRR1039523
## ENSG00000288106          5
## ENSG00000288107          0
## ENSG00000288108          0
## ENSG00000288109          0
## ENSG00000288110          0
## ENSG00000288111          0
```

Lets drop these last 5 rows as we don't need them. 

```r
cts <- cts[-((60622-5):60622),]
dim(cts)
```

```
## [1] 60616    16
```

If you estimated transcript counts (rather than gene-level counts produced by htseq-count) using a method like *RSEM*, *Salmon*, or *kallisto*, using the tximport() function from the [tximport package](https://f1000research.com/articles/4-1521/v1). You may have estimated transcript-level counts if you used a library-preparation protocol that captures full length transcript information. Even if you only plan to do a differential expression analysis at the gene-level, it has been shown that [transcript-level estimates can improve gene-level inferences](https://f1000research.com/articles/4-1521/v1), therefore if you are able to estimate counts at the transcript-level for your data, it is beneficial to do so. Briefly, this method works by collapsing transcript-level estimates into gene-level estimates, while an offset matrix is calculated based on the average transcript length, that is used in the differential expression analysis to correct for biases that may be introduced by transcript-length differences between samples. You can read more about how to do this in the [documnetation for tximport](https://bioconductor.org/packages/release/bioc/vignettes/tximport/inst/doc/tximport.html).

If you collected 3'-end data, e.g. with the Lexogen QuantSeq assay, you should not do such a correction for length, as there is no length bias in your data. Doing this correction would introduce bias into your data and likely distort your differential expression results. For 3'-end data, it is best to read in the raw count matrix directly using (`DESeqDataSetFromHTSeqCount`) or simply (`read.table()`). 

### Read in sample metadata

We also need to read in the sample annotation (metadata) that we downloaded from the SRA, which contains sample labels, experimental labels, and sequencing run information, etc. 

```r
# read in the file from the SRA that has sample/experimental labels 
sra_res <- read.csv("Day-2/sra_result.csv", row.names=1)
head(sra_res)
```

```
##                                               Experiment.Title Organism.Name
## SRX384360   GSM1275877: N061011_Alb_Dex; Homo sapiens; RNA-Seq  Homo sapiens
## SRX384359       GSM1275876: N061011_Alb; Homo sapiens; RNA-Seq  Homo sapiens
## SRX384358       GSM1275875: N061011_Dex; Homo sapiens; RNA-Seq  Homo sapiens
## SRX384357 GSM1275874: N061011_untreated; Homo sapiens; RNA-Seq  Homo sapiens
## SRX384356   GSM1275873: N080611_Alb_Dex; Homo sapiens; RNA-Seq  Homo sapiens
## SRX384355       GSM1275872: N080611_Alb; Homo sapiens; RNA-Seq  Homo sapiens
##                    Instrument Submitter Study.Accession
## SRX384360 Illumina HiSeq 2000       GEO       SRP033351
## SRX384359 Illumina HiSeq 2000       GEO       SRP033351
## SRX384358 Illumina HiSeq 2000       GEO       SRP033351
## SRX384357 Illumina HiSeq 2000       GEO       SRP033351
## SRX384356 Illumina HiSeq 2000       GEO       SRP033351
## SRX384355 Illumina HiSeq 2000       GEO       SRP033351
##                                                                                  Study.Title
## SRX384360 Human Airway Smooth Muscle Transcriptome Changes in Response to Asthma Medications
## SRX384359 Human Airway Smooth Muscle Transcriptome Changes in Response to Asthma Medications
## SRX384358 Human Airway Smooth Muscle Transcriptome Changes in Response to Asthma Medications
## SRX384357 Human Airway Smooth Muscle Transcriptome Changes in Response to Asthma Medications
## SRX384356 Human Airway Smooth Muscle Transcriptome Changes in Response to Asthma Medications
## SRX384355 Human Airway Smooth Muscle Transcriptome Changes in Response to Asthma Medications
##           Sample.Accession Sample.Title Total.Size..Mb Total.RUNs Total.Spots
## SRX384360        SRS508581           NA        2269.81          1    31099538
## SRX384359        SRS508582           NA        2006.00          1    28406812
## SRX384358        SRS508580           NA        2449.87          1    41152075
## SRX384357        SRS508579           NA        2109.13          1    34575286
## SRX384356        SRS508577           NA        2006.23          1    31533740
## SRX384355        SRS508578           NA        2354.22          1    30441200
##           Total.Bases Library.Name Library.Strategy Library.Source
## SRX384360  3918541788           NA          RNA-Seq TRANSCRIPTOMIC
## SRX384359  3562039593           NA          RNA-Seq TRANSCRIPTOMIC
## SRX384358  4072315905           NA          RNA-Seq TRANSCRIPTOMIC
## SRX384357  3518623962           NA          RNA-Seq TRANSCRIPTOMIC
## SRX384356  3380181525           NA          RNA-Seq TRANSCRIPTOMIC
## SRX384355  3835591200           NA          RNA-Seq TRANSCRIPTOMIC
##           Library.Selection
## SRX384360              cDNA
## SRX384359              cDNA
## SRX384358              cDNA
## SRX384357              cDNA
## SRX384356              cDNA
## SRX384355              cDNA
```

```r
# rename the Sample.Accession column as Sample so that we can merge with the run info, which we will load next 
sra_res$Sample <- sra_res$Sample.Accession
head(sra_res)
```

```
##                                               Experiment.Title Organism.Name
## SRX384360   GSM1275877: N061011_Alb_Dex; Homo sapiens; RNA-Seq  Homo sapiens
## SRX384359       GSM1275876: N061011_Alb; Homo sapiens; RNA-Seq  Homo sapiens
## SRX384358       GSM1275875: N061011_Dex; Homo sapiens; RNA-Seq  Homo sapiens
## SRX384357 GSM1275874: N061011_untreated; Homo sapiens; RNA-Seq  Homo sapiens
## SRX384356   GSM1275873: N080611_Alb_Dex; Homo sapiens; RNA-Seq  Homo sapiens
## SRX384355       GSM1275872: N080611_Alb; Homo sapiens; RNA-Seq  Homo sapiens
##                    Instrument Submitter Study.Accession
## SRX384360 Illumina HiSeq 2000       GEO       SRP033351
## SRX384359 Illumina HiSeq 2000       GEO       SRP033351
## SRX384358 Illumina HiSeq 2000       GEO       SRP033351
## SRX384357 Illumina HiSeq 2000       GEO       SRP033351
## SRX384356 Illumina HiSeq 2000       GEO       SRP033351
## SRX384355 Illumina HiSeq 2000       GEO       SRP033351
##                                                                                  Study.Title
## SRX384360 Human Airway Smooth Muscle Transcriptome Changes in Response to Asthma Medications
## SRX384359 Human Airway Smooth Muscle Transcriptome Changes in Response to Asthma Medications
## SRX384358 Human Airway Smooth Muscle Transcriptome Changes in Response to Asthma Medications
## SRX384357 Human Airway Smooth Muscle Transcriptome Changes in Response to Asthma Medications
## SRX384356 Human Airway Smooth Muscle Transcriptome Changes in Response to Asthma Medications
## SRX384355 Human Airway Smooth Muscle Transcriptome Changes in Response to Asthma Medications
##           Sample.Accession Sample.Title Total.Size..Mb Total.RUNs Total.Spots
## SRX384360        SRS508581           NA        2269.81          1    31099538
## SRX384359        SRS508582           NA        2006.00          1    28406812
## SRX384358        SRS508580           NA        2449.87          1    41152075
## SRX384357        SRS508579           NA        2109.13          1    34575286
## SRX384356        SRS508577           NA        2006.23          1    31533740
## SRX384355        SRS508578           NA        2354.22          1    30441200
##           Total.Bases Library.Name Library.Strategy Library.Source
## SRX384360  3918541788           NA          RNA-Seq TRANSCRIPTOMIC
## SRX384359  3562039593           NA          RNA-Seq TRANSCRIPTOMIC
## SRX384358  4072315905           NA          RNA-Seq TRANSCRIPTOMIC
## SRX384357  3518623962           NA          RNA-Seq TRANSCRIPTOMIC
## SRX384356  3380181525           NA          RNA-Seq TRANSCRIPTOMIC
## SRX384355  3835591200           NA          RNA-Seq TRANSCRIPTOMIC
##           Library.Selection    Sample
## SRX384360              cDNA SRS508581
## SRX384359              cDNA SRS508582
## SRX384358              cDNA SRS508580
## SRX384357              cDNA SRS508579
## SRX384356              cDNA SRS508577
## SRX384355              cDNA SRS508578
```

```r
# read in the file that contains specific info for each sequencing run 
sra_run <- read.csv("Day-2/SraRunInfo.csv", row.names=1)
head(sra_run)
```

```
##            ReleaseDate       LoadDate    spots      bases spots_with_mates
## SRR1039508 1/2/14 9:16 11/26/13 16:39 22935521 2889875646         22935521
## SRR1039509 1/2/14 9:16 11/26/13 16:38 21155707 2665619082         21155707
## SRR1039510 1/2/14 9:16 11/26/13 16:41 22852619 2879429994         22852619
## SRR1039511 1/2/14 9:16 11/26/13 16:42 21938637 2764268262         21938637
## SRR1039512 1/2/14 9:16 11/26/13 16:47 28136282 3545171532         28136282
## SRR1039513 1/2/14 9:16 11/26/13 17:02 43356464 3791311776         16823088
##            avgLength size_MB AssemblyName
## SRR1039508       126    1588           NA
## SRR1039509       126    1480           NA
## SRR1039510       126    1593           NA
## SRR1039511       126    1519           NA
## SRR1039512       126    2055           NA
## SRR1039513        87    2271           NA
##                                                                                      download_path
## SRR1039508 https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos1/sra-pub-run-5/SRR1039508/SRR1039508.1
## SRR1039509 https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos1/sra-pub-run-5/SRR1039509/SRR1039509.1
## SRR1039510 https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos1/sra-pub-run-5/SRR1039510/SRR1039510.1
## SRR1039511 https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos1/sra-pub-run-5/SRR1039511/SRR1039511.1
## SRR1039512 https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos1/sra-pub-run-5/SRR1039512/SRR1039512.1
## SRR1039513 https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos1/sra-pub-run-5/SRR1039513/SRR1039513.1
##            Experiment LibraryName LibraryStrategy LibrarySelection
## SRR1039508  SRX384345          NA         RNA-Seq             cDNA
## SRR1039509  SRX384346          NA         RNA-Seq             cDNA
## SRR1039510  SRX384347          NA         RNA-Seq             cDNA
## SRR1039511  SRX384348          NA         RNA-Seq             cDNA
## SRR1039512  SRX384349          NA         RNA-Seq             cDNA
## SRR1039513  SRX384350          NA         RNA-Seq             cDNA
##             LibrarySource LibraryLayout InsertSize InsertDev Platform
## SRR1039508 TRANSCRIPTOMIC        PAIRED          0         0 ILLUMINA
## SRR1039509 TRANSCRIPTOMIC        PAIRED          0         0 ILLUMINA
## SRR1039510 TRANSCRIPTOMIC        PAIRED          0         0 ILLUMINA
## SRR1039511 TRANSCRIPTOMIC        PAIRED          0         0 ILLUMINA
## SRR1039512 TRANSCRIPTOMIC        PAIRED          0         0 ILLUMINA
## SRR1039513 TRANSCRIPTOMIC        PAIRED          0         0 ILLUMINA
##                          Model  SRAStudy  BioProject Study_Pubmed_id ProjectID
## SRR1039508 Illumina HiSeq 2000 SRP033351 PRJNA229998               2    229998
## SRR1039509 Illumina HiSeq 2000 SRP033351 PRJNA229998               2    229998
## SRR1039510 Illumina HiSeq 2000 SRP033351 PRJNA229998               2    229998
## SRR1039511 Illumina HiSeq 2000 SRP033351 PRJNA229998               2    229998
## SRR1039512 Illumina HiSeq 2000 SRP033351 PRJNA229998               2    229998
## SRR1039513 Illumina HiSeq 2000 SRP033351 PRJNA229998               2    229998
##               Sample    BioSample SampleType TaxID ScientificName SampleName
## SRR1039508 SRS508568 SAMN02422669     simple  9606   Homo sapiens GSM1275862
## SRR1039509 SRS508567 SAMN02422675     simple  9606   Homo sapiens GSM1275863
## SRR1039510 SRS508570 SAMN02422668     simple  9606   Homo sapiens GSM1275864
## SRR1039511 SRS508569 SAMN02422667     simple  9606   Homo sapiens GSM1275865
## SRR1039512 SRS508571 SAMN02422678     simple  9606   Homo sapiens GSM1275866
## SRR1039513 SRS508572 SAMN02422670     simple  9606   Homo sapiens GSM1275867
##            g1k_pop_code source g1k_analysis_group Subject_ID Sex Disease Tumor
## SRR1039508           NA     NA                 NA         NA  NA      NA    no
## SRR1039509           NA     NA                 NA         NA  NA      NA    no
## SRR1039510           NA     NA                 NA         NA  NA      NA    no
## SRR1039511           NA     NA                 NA         NA  NA      NA    no
## SRR1039512           NA     NA                 NA         NA  NA      NA    no
## SRR1039513           NA     NA                 NA         NA  NA      NA    no
##            Affection_Status Analyte_Type Histological_Type Body_Site CenterName
## SRR1039508               NA           NA                NA        NA        GEO
## SRR1039509               NA           NA                NA        NA        GEO
## SRR1039510               NA           NA                NA        NA        GEO
## SRR1039511               NA           NA                NA        NA        GEO
## SRR1039512               NA           NA                NA        NA        GEO
## SRR1039513               NA           NA                NA        NA        GEO
##            Submission dbgap_study_accession Consent
## SRR1039508  SRA114259                    NA  public
## SRR1039509  SRA114259                    NA  public
## SRR1039510  SRA114259                    NA  public
## SRR1039511  SRA114259                    NA  public
## SRR1039512  SRA114259                    NA  public
## SRR1039513  SRA114259                    NA  public
##                                     RunHash                         ReadHash
## SRR1039508 65345B22A4745B616DD5840CF2AD90EE 54B084AA7AD2A02B39D37CBBD5B1C3DC
## SRR1039509 E2F5FBD3278F29EC6EA7B00FC95C28B2 294D13FE5364BA88B5FEAD68A53CB11B
## SRR1039510 2019240C219172C77646ED066033D588 16F7F5A040CECD58D73F344710F13D7D
## SRR1039511 D7FCAEEB12F4F5EE4956EDED736E4E4F 0D8679513EAC88ED62C14B488E2134E5
## SRR1039512 828A997474D0BBD6E5E14DEDB0D9E31D 5834A34CA3CE7EF962A396098543C9E6
## SRR1039513 0EB840FF1F4AA22A5017CB754C039225 5A5DD90EDCD5998A26F7F7BC099D4010
```

```r
# lets set the row names (SRR numbers from the SRA, stands for SRA run accession) as a new variable in the dataframe, as we will need it later 
sra_run$SRR <- rownames(sra_run)

# merge sra_res and sra_run together so that we get all of the identifiers and experiment info that we will need in the same dataframe
colData <- merge(sra_res, sra_run, merge="Sample")

# order by SRA run accession 
colData <- colData[order(colData$SRR),]

# quick look 
head(colData)
```

```
##      Sample                                     Experiment.Title Organism.Name
## 2 SRS508568  GSM1275862: N61311_untreated; Homo sapiens; RNA-Seq  Homo sapiens
## 1 SRS508567        GSM1275863: N61311_Dex; Homo sapiens; RNA-Seq  Homo sapiens
## 4 SRS508570        GSM1275864: N61311_Alb; Homo sapiens; RNA-Seq  Homo sapiens
## 3 SRS508569    GSM1275865: N61311_Alb_Dex; Homo sapiens; RNA-Seq  Homo sapiens
## 5 SRS508571 GSM1275866: N052611_untreated; Homo sapiens; RNA-Seq  Homo sapiens
## 6 SRS508572       GSM1275867: N052611_Dex; Homo sapiens; RNA-Seq  Homo sapiens
##            Instrument Submitter Study.Accession
## 2 Illumina HiSeq 2000       GEO       SRP033351
## 1 Illumina HiSeq 2000       GEO       SRP033351
## 4 Illumina HiSeq 2000       GEO       SRP033351
## 3 Illumina HiSeq 2000       GEO       SRP033351
## 5 Illumina HiSeq 2000       GEO       SRP033351
## 6 Illumina HiSeq 2000       GEO       SRP033351
##                                                                          Study.Title
## 2 Human Airway Smooth Muscle Transcriptome Changes in Response to Asthma Medications
## 1 Human Airway Smooth Muscle Transcriptome Changes in Response to Asthma Medications
## 4 Human Airway Smooth Muscle Transcriptome Changes in Response to Asthma Medications
## 3 Human Airway Smooth Muscle Transcriptome Changes in Response to Asthma Medications
## 5 Human Airway Smooth Muscle Transcriptome Changes in Response to Asthma Medications
## 6 Human Airway Smooth Muscle Transcriptome Changes in Response to Asthma Medications
##   Sample.Accession Sample.Title Total.Size..Mb Total.RUNs Total.Spots
## 2        SRS508568           NA        1588.22          1    22935521
## 1        SRS508567           NA        1480.19          1    21155707
## 4        SRS508570           NA        1593.15          1    22852619
## 3        SRS508569           NA        1519.57          1    21938637
## 5        SRS508571           NA        2055.72          1    28136282
## 6        SRS508572           NA        2271.23          1    43356464
##   Total.Bases Library.Name Library.Strategy Library.Source Library.Selection
## 2  2889875646           NA          RNA-Seq TRANSCRIPTOMIC              cDNA
## 1  2665619082           NA          RNA-Seq TRANSCRIPTOMIC              cDNA
## 4  2879429994           NA          RNA-Seq TRANSCRIPTOMIC              cDNA
## 3  2764268262           NA          RNA-Seq TRANSCRIPTOMIC              cDNA
## 5  3545171532           NA          RNA-Seq TRANSCRIPTOMIC              cDNA
## 6  3791311776           NA          RNA-Seq TRANSCRIPTOMIC              cDNA
##   ReleaseDate       LoadDate    spots      bases spots_with_mates avgLength
## 2 1/2/14 9:16 11/26/13 16:39 22935521 2889875646         22935521       126
## 1 1/2/14 9:16 11/26/13 16:38 21155707 2665619082         21155707       126
## 4 1/2/14 9:16 11/26/13 16:41 22852619 2879429994         22852619       126
## 3 1/2/14 9:16 11/26/13 16:42 21938637 2764268262         21938637       126
## 5 1/2/14 9:16 11/26/13 16:47 28136282 3545171532         28136282       126
## 6 1/2/14 9:16 11/26/13 17:02 43356464 3791311776         16823088        87
##   size_MB AssemblyName
## 2    1588           NA
## 1    1480           NA
## 4    1593           NA
## 3    1519           NA
## 5    2055           NA
## 6    2271           NA
##                                                                             download_path
## 2 https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos1/sra-pub-run-5/SRR1039508/SRR1039508.1
## 1 https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos1/sra-pub-run-5/SRR1039509/SRR1039509.1
## 4 https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos1/sra-pub-run-5/SRR1039510/SRR1039510.1
## 3 https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos1/sra-pub-run-5/SRR1039511/SRR1039511.1
## 5 https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos1/sra-pub-run-5/SRR1039512/SRR1039512.1
## 6 https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos1/sra-pub-run-5/SRR1039513/SRR1039513.1
##   Experiment LibraryName LibraryStrategy LibrarySelection  LibrarySource
## 2  SRX384345          NA         RNA-Seq             cDNA TRANSCRIPTOMIC
## 1  SRX384346          NA         RNA-Seq             cDNA TRANSCRIPTOMIC
## 4  SRX384347          NA         RNA-Seq             cDNA TRANSCRIPTOMIC
## 3  SRX384348          NA         RNA-Seq             cDNA TRANSCRIPTOMIC
## 5  SRX384349          NA         RNA-Seq             cDNA TRANSCRIPTOMIC
## 6  SRX384350          NA         RNA-Seq             cDNA TRANSCRIPTOMIC
##   LibraryLayout InsertSize InsertDev Platform               Model  SRAStudy
## 2        PAIRED          0         0 ILLUMINA Illumina HiSeq 2000 SRP033351
## 1        PAIRED          0         0 ILLUMINA Illumina HiSeq 2000 SRP033351
## 4        PAIRED          0         0 ILLUMINA Illumina HiSeq 2000 SRP033351
## 3        PAIRED          0         0 ILLUMINA Illumina HiSeq 2000 SRP033351
## 5        PAIRED          0         0 ILLUMINA Illumina HiSeq 2000 SRP033351
## 6        PAIRED          0         0 ILLUMINA Illumina HiSeq 2000 SRP033351
##    BioProject Study_Pubmed_id ProjectID    BioSample SampleType TaxID
## 2 PRJNA229998               2    229998 SAMN02422669     simple  9606
## 1 PRJNA229998               2    229998 SAMN02422675     simple  9606
## 4 PRJNA229998               2    229998 SAMN02422668     simple  9606
## 3 PRJNA229998               2    229998 SAMN02422667     simple  9606
## 5 PRJNA229998               2    229998 SAMN02422678     simple  9606
## 6 PRJNA229998               2    229998 SAMN02422670     simple  9606
##   ScientificName SampleName g1k_pop_code source g1k_analysis_group Subject_ID
## 2   Homo sapiens GSM1275862           NA     NA                 NA         NA
## 1   Homo sapiens GSM1275863           NA     NA                 NA         NA
## 4   Homo sapiens GSM1275864           NA     NA                 NA         NA
## 3   Homo sapiens GSM1275865           NA     NA                 NA         NA
## 5   Homo sapiens GSM1275866           NA     NA                 NA         NA
## 6   Homo sapiens GSM1275867           NA     NA                 NA         NA
##   Sex Disease Tumor Affection_Status Analyte_Type Histological_Type Body_Site
## 2  NA      NA    no               NA           NA                NA        NA
## 1  NA      NA    no               NA           NA                NA        NA
## 4  NA      NA    no               NA           NA                NA        NA
## 3  NA      NA    no               NA           NA                NA        NA
## 5  NA      NA    no               NA           NA                NA        NA
## 6  NA      NA    no               NA           NA                NA        NA
##   CenterName Submission dbgap_study_accession Consent
## 2        GEO  SRA114259                    NA  public
## 1        GEO  SRA114259                    NA  public
## 4        GEO  SRA114259                    NA  public
## 3        GEO  SRA114259                    NA  public
## 5        GEO  SRA114259                    NA  public
## 6        GEO  SRA114259                    NA  public
##                            RunHash                         ReadHash        SRR
## 2 65345B22A4745B616DD5840CF2AD90EE 54B084AA7AD2A02B39D37CBBD5B1C3DC SRR1039508
## 1 E2F5FBD3278F29EC6EA7B00FC95C28B2 294D13FE5364BA88B5FEAD68A53CB11B SRR1039509
## 4 2019240C219172C77646ED066033D588 16F7F5A040CECD58D73F344710F13D7D SRR1039510
## 3 D7FCAEEB12F4F5EE4956EDED736E4E4F 0D8679513EAC88ED62C14B488E2134E5 SRR1039511
## 5 828A997474D0BBD6E5E14DEDB0D9E31D 5834A34CA3CE7EF962A396098543C9E6 SRR1039512
## 6 0EB840FF1F4AA22A5017CB754C039225 5A5DD90EDCD5998A26F7F7BC099D4010 SRR1039513
```

Now that we have the sample metadata loaded, we need to clean up some variables to get the information we need for downstream analysis. In this case the experiment title tells us which group each sample belongs, this is important so that we can make the correct comparisons later in our analysis. This information is used to create our experimental design.

The design formula lets DESeq2 know your experimental design and what groups you want to compare in the differential expression analysis. It requires a variable that is present in colData to be specified, and will be used for the statistical analysis in DESeq2. This variable should be a (`factor`) class variable, with the reference group for your comparison set as the first level of the variable.


```r
# lets chop up the 'Experiment.Title' variable to extrat the experimental group labels for the experiment 
colData$group <- sapply(as.character(colData$Experiment.Title), function(x) strsplit(x, "1_")[[1]][2])
colData$group
```

```
##  [1] "untreated; Homo sapiens; RNA-Seq" "Dex; Homo sapiens; RNA-Seq"      
##  [3] "Alb; Homo sapiens; RNA-Seq"       "Alb_Dex; Homo sapiens; RNA-Seq"  
##  [5] "untreated; Homo sapiens; RNA-Seq" "Dex; Homo sapiens; RNA-Seq"      
##  [7] "Alb; Homo sapiens; RNA-Seq"       "Alb_Dex; Homo sapiens; RNA-Seq"  
##  [9] "untreated; Homo sapiens; RNA-Seq" "Dex; Homo sapiens; RNA-Seq"      
## [11] "Alb; Homo sapiens; RNA-Seq"       "Alb_Dex; Homo sapiens; RNA-Seq"  
## [13] "untreated; Homo sapiens; RNA-Seq" "Dex; Homo sapiens; RNA-Seq"      
## [15] "Alb; Homo sapiens; RNA-Seq"       "Alb_Dex; Homo sapiens; RNA-Seq"
```

```r
colData$group <- sapply(as.character(colData$group), function(x) strsplit(x, ";")[[1]][1])
colData$group
```

```
##  [1] "untreated" "Dex"       "Alb"       "Alb_Dex"   "untreated" "Dex"      
##  [7] "Alb"       "Alb_Dex"   "untreated" "Dex"       "Alb"       "Alb_Dex"  
## [13] "untreated" "Dex"       "Alb"       "Alb_Dex"
```

```r
# now make this a factor as it will be the variable we will use define groups for the differential expression analysis 
colData$group <- factor(colData$group, levels=c("untreated", "Dex", "Alb", "Alb_Dex"))
colData$group
```

```
##  [1] untreated Dex       Alb       Alb_Dex   untreated Dex       Alb      
##  [8] Alb_Dex   untreated Dex       Alb       Alb_Dex   untreated Dex      
## [15] Alb       Alb_Dex  
## Levels: untreated Dex Alb Alb_Dex
```

### Construct the DESEq2 data set & explore the characteristics of the data 

DESeq2 uses an object class called the (`DESeqDataSet`) that stores the read counts, metadata, experimental design, and all the intermediate values calculated during the analysis. (`DESeqDataSet`) extends the (`SummarizedExperiment`) class object from the (`SummarizedExperiment`) R/Bioconductor package that is commonly used to store data from expression studies and other genomics assays in R. 

Three elements are required to generate the (`DESeqDataSet`):  
- matrix of raw counts
- sample metadata (colData)
- a design formula 
 

Lets create the (`DESeqDataSet`) object. 

```r
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = colData,
                              design = ~ group)
```

We could have also done this using the (`DESeqDataSetFromHTSeqCount()`) function by specifying a (`SampleTable`) that includes the path to the htseq-count files, however since we compiled the read counts into one file, we can just load the dataset directly. 

Before moving on, lets explore our DESeq2 class object a bit to get to familar with its contents. 

```r
# have a quick look at the object 
dds

# print structure 
str(dds)

# several accessor functions exist to access specific data 'slots'
head(counts(dds))
head(colData(dds))

# specific slots can also be accessed using the '@'
dds@colData
```

Lets drop genes that have less than 10 reads across all samples, as there just isn't enough 
information for these genes to fit robust statistical models to. 

```r
# drop genes with low counts 
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
```

Lets also save the DESeq object at this point (so that we don't have to do the above everytime we want to work with our data). 

```r
save(dds, file = "DESeq2.rdata")
```

### Normalization of raw counts 

Before comparing expression levels of specific genes between samples, and performing any differential expression analyses, it is important that we normalize the data to account for variation in expression that is not related to true differential expression. There are two major sources variation that we need to adjust for in the normalization process for RNA-seq data when we wish to compare expression levels **between** samples:  

#### Library size/sequencing depth  

Although we generally try to pool samples together (each sample is tagged with a barcode and samples are combined) at similar concentrations in a sequencing run, some samples will end up being sequenced more than others, leading to slight differences in how many reads are produced for that sample, and therefore sequencing depth and size. Furthermore, if samples are sequenced on separate runs, their sequencing depths may be very different. If we don't account for this variation in sequencing depth, we might conclude some genes are expressed at greater levels in a sample that has simply been sequenced to a higher depth.  

![](/Users/OwenW/Downloads/RNA-seq_workshop_July2020-master/figures/library_size.png)

#### Library composition 

The presence of truly differentially expressed genes (in particular, DEGs with very large fold changes) between samples will cause the number of reads for other genes in those samples to be skewed. For example, in the below example, gene C is differentially expressed between the two samples, with much higher expression in sample 1. This high number of reads causes fewer reads to be detected for other genes in this sample, making it appear that these other genes are expressed at lower levels than in sample 2, however this is simply an artifact of library composition differences between the samples. 

![](/Users/OwenW/Downloads/RNA-seq_workshop_July2020-master/figures/library_composition.png)

To correct for **library size** AND **library composition**, DESeq2 uses a algorithm referred to as the median-of-ratios method. Although we won't go over how the algorithm works in detail, a brief summary of the steps is:  

1. Take the log of all values in raw count matrix  
2. Average each row (genes)
3. Filter out genes with Infinity values
4. Subtract average log count value from log of count for each cell (due to the laws of working with logarithms, this is essentially calculating the ratio of the counts for gene X in 1 sample to the average counts for gene X across all samples)
5. Calculate the median of the ratios in each sample (column)
6. Take exponents of medians to get the **size factors** for each sample/library. 
7. Divide the count for each gene in each sample by the size factor calculated for that sample. 

This procedure will generate a matrix of read counts that are corrected for both **library size** and **library composition**, and are stored in our (`DESeqDataset`) object. DESeq2 uses the function (`estimateSizeFactors()`) to perform this algorithm and calculate size factors for each sample. Lets do this for our (`DESeqDataset`). 

```r
dds <- estimateSizeFactors(dds)
```

Note: [This video](https://www.youtube.com/watch?v=UFB993xufUU) from StatQuest provides an excellent summary of the steps performed by (`estimateSizeFactors()`) in order to calculate these size factors.  

Once we have calculated the size factors, it can be helpful to look at their distribution to get a feel for how they vary and how much normalization between the samples is required. 

```r
sizeFactors(dds)
```

```
## SRR1039508 SRR1039509 SRR1039510 SRR1039511 SRR1039512 SRR1039513 SRR1039514 
##  0.9546639  0.8379264  0.9559347  0.8841828  1.1063105  0.6305985  1.5002929 
## SRR1039515 SRR1039516 SRR1039517 SRR1039518 SRR1039519 SRR1039520 SRR1039521 
##  0.9404767  1.1086740  1.3092383  1.1900257  0.8972838  0.8647837  0.8868375 
## SRR1039522 SRR1039523 
##  1.1550681  1.2038272
```

```r
hist(sizeFactors(dds), 
     breaks=6, col = "cornflowerblue",
     xlab="Size factors", ylab="No. of samples", 
     main= "Size factor distribution over samples")
```

<img src="day2-PART1_files/figure-html/unnamed-chunk-11-1.png" style="display: block; margin: auto;" />

After we have calculated the size factors, we can use the `counts()` function, with `normalized` set to `TRUE`), to return the matrix of counts where each column (each library/sample) have been divided by the size factors calculated by the `estimateSizeFactors()` function. 

```r
counts_norm <- counts(dds, normalized=TRUE)
head(counts_norm)
```

```
##                 SRR1039508 SRR1039509 SRR1039510 SRR1039511 SRR1039512
## ENSG00000000003  731.14734  558.52161   792.9412   532.6952  807.18751
## ENSG00000000419  486.03491  605.06508   475.9739   588.1136  549.57447
## ENSG00000000457  270.25217  245.84498   247.9249   256.7342  236.82321
## ENSG00000000460   60.75436   65.63822    92.0565    64.4663   35.25231
## ENSG00000000971 3048.19305 3935.90655  3595.4340  4330.5522 5038.36860
## ENSG00000001036 1465.43714 1243.54597  1450.9360  1104.9751 1528.50400
##                 SRR1039513 SRR1039514 SRR1039515 SRR1039516 SRR1039517
## ENSG00000000003  666.03391  656.53845  933.56907 1064.33448  820.32433
## ENSG00000000419  566.12882  579.22022  475.29086  523.14746  596.53008
## ENSG00000000457  258.48459  218.62397  230.73404  218.27877  252.05496
## ENSG00000000460   55.50283   57.98868   49.97465   69.45233   48.88339
## ENSG00000000971 6041.87900 7435.21467 6154.32548 5468.69486 7524.22254
## ENSG00000001036 1368.54110 1240.42443 1601.31551 1259.16181 1070.85169
##                 SRR1039518 SRR1039519 SRR1039520 SRR1039521 SRR1039522
## ENSG00000000003  936.11422  1271.6156  927.39954  670.92334  973.96853
## ENSG00000000419  602.50799   568.3821  478.73243  563.80113  463.17614
## ENSG00000000457  247.05348   221.7805  265.96246  258.22092  214.70595
## ENSG00000000460   55.46099    74.6698   87.88325   67.65614   72.72298
## ENSG00000000971 8363.68479  6652.2991 5330.81279 8039.80410 6509.57275
## ENSG00000001036 1052.91842  1292.7906 1535.64412 1223.44845 1523.71965
##                 SRR1039523
## ENSG00000000003  724.35645
## ENSG00000000419  559.05033
## ENSG00000000457  230.09947
## ENSG00000000460   75.59224
## ENSG00000000971 7908.94239
## ENSG00000001036 1145.51324
```

Comparing the normalized to the raw counts, we can clearly see that they are different. 

```r
head(counts(dds, normalized=FALSE))
```

```
##                 SRR1039508 SRR1039509 SRR1039510 SRR1039511 SRR1039512
## ENSG00000000003        698        468        758        471        893
## ENSG00000000419        464        507        455        520        608
## ENSG00000000457        258        206        237        227        262
## ENSG00000000460         58         55         88         57         39
## ENSG00000000971       2910       3298       3437       3829       5574
## ENSG00000001036       1399       1042       1387        977       1691
##                 SRR1039513 SRR1039514 SRR1039515 SRR1039516 SRR1039517
## ENSG00000000003        420        985        878       1180       1074
## ENSG00000000419        357        869        447        580        781
## ENSG00000000457        163        328        217        242        330
## ENSG00000000460         35         87         47         77         64
## ENSG00000000971       3810      11155       5788       6063       9851
## ENSG00000001036        863       1861       1506       1396       1402
##                 SRR1039518 SRR1039519 SRR1039520 SRR1039521 SRR1039522
## ENSG00000000003       1114       1141        802        595       1125
## ENSG00000000419        717        510        414        500        535
## ENSG00000000457        294        199        230        229        248
## ENSG00000000460         66         67         76         60         84
## ENSG00000000971       9953       5969       4610       7130       7519
## ENSG00000001036       1253       1160       1328       1085       1760
##                 SRR1039523
## ENSG00000000003        872
## ENSG00000000419        673
## ENSG00000000457        277
## ENSG00000000460         91
## ENSG00000000971       9521
## ENSG00000001036       1379
```

We can use this table of normalized read counts to compare values for individual genes across samples. We might want to use this to (sanity) check the expression of a few genes of interest, before we actually do any statistical modelling. The abstract of the paper describes *DUSP1*, a phosphatase with dual specificity for tyrosine and threonine, as a well-known glucocorticoid-responsive gene. 

```r
# lets make a function to generate a quick plot of the normalized counts 
gene_plot <- function(ENSG, gene_symbol){
  # save the normalized counts in a dataframe 
  cnts <- counts(dds, normalized=TRUE)
  colnames(cnts) <- colData(dds)$SRR
  # extract the counts for specified ENSG ID and add sample group data 
  df1 <- data.frame(log2(cnts[ENSG,]), colData(dds)$group)
  colnames(df1) <- c(paste0("log2_gene"), "sample_group")
  # use ggplot2 to make a plot of counts vs sample group 
  p1<- ggplot(df1, aes(sample_group, log2_gene)) + 
    geom_jitter(aes(color = sample_group)) + 
    ggtitle(paste0(gene_symbol), " - Log2 Normalized counts")
  # print the plot 
  print(p1)
}

# now apply the function to print a plot for a specified gene 
gene_plot(ENSG = "ENSG00000120129", gene_symbol = "DUSP1")
```

<img src="day2-PART1_files/figure-html/unnamed-chunk-14-1.png" style="display: block; margin: auto;" />

DUSP1 expression is consistently higher in the DEX samples than the untreated, suggesting this gene is differentially expressed after DEX treatment, validating prior knowledge and giving us confidence that our experiment worked, sample labels are all correct, and we are well positioned to make new discoveries with these data. 

**Important note:** the normalized count matrix is normalized for **library size and composition**, which means we can compare expression levels of individual genes across samples. The read counts are NOT normalized for gene length, so we cannot use this matrix to compare expression levels between genes within the same sample. This is important because some genes may simply pick up more reads than others because they are larger, making them appear more highly expressed than a smaller gene, which may not be the case. For such comparisons between genes, we need to use measures such as counts per million (CPM), transcripts per million (TPM), or fragments per kilobase million (FPKM)/reads per kilobase million (RPKM). [This video](https://www.rna-seqblog.com/rpkm-fpkm-and-tpm-clearly-explained/) provides an excellent explanation of RPKM, FRPM, & TPM, and explains why it is better to use TPM if you need to correct for library size AND gene length. 

<center>
![](/Users/OwenW/Downloads/RNA-seq_workshop_July2020-master/figures/gene_length.png)
</center>


### Exploratory data analysis & quality control

Before we run the differential expression analysis, we should explore our dataset to learn a little more about it. Importantly, we should determine how samples are related to each other based on their gene expression profiles, ensure replicates cluster together, assess the potential for any batch effects (technical variation between samples run on different machines or by different technicians) that may exist, and more generally build expectations for the DE analysis. For these exploratory analyses, we use several unsupervised methods to assess how our samples cluster relative to one another. These unsupervised approached include **principal components analysis (PCA)**, and **unsupervised Hierarchical clustering** and are commonly referred to as *data reduction* appraoches. 

#### Principal components analysis (PCA)

PCA is a mathematical procedure that calculates vectors that explain varition in the dataset (in this case, variation in gene expression), and orders samples along these vectors. We would expect samples that are more similar to each other, e.g. replicates, to be very close to each other along these axes of varaition, while we might expect samples in different treatment groups to be further away. Each vector that is calculated is called a principal component (PC), and each principal component explains less varaition in our dataset than the last. e.g, PC1 explains more variation in gene expression differences between the samples than PC2. If we plot PC1 against PC2, samples will 'cluster' with other samples that have similar gene expression profiles, and be further away from samples with more distant expression profiles. 

<center>
![](/Users/OwenW/Downloads/RNA-seq_workshop_July2020-master/figures/pca_example.png)
</center>


Again, StatQuest has an excellent [video](https://www.youtube.com/watch?v=_UVHneBUBW0) that explains the fundamental concepts of PCA, and provoides more details how to the PCs themselves are calculated. 

To perform mathematical procedures such as PCA, it is better to work with a transformed version of the counts, rather than the original untransformed values. DESeq2 actually utilizes its own transformation procedure, called the regularized logatrithm (rlog) implemented in the `rlog()` function. The rlog is similar in principle to a standard log transformation of the data, but is able to more appropriately transform the counts for genes with low expression values. 

DESeq2 is also capable of implementing a variance stabilising transformation (VST) for count data, this is generally recommended for larger datasets due to increased speed. 

For this analysis, we will use the rlog, which produces values on the log2 scale that are also normalized for library size during the rlog procedure. Lets perform the rlog transformation on our data. 

```r
rld <- rlog(dds, blind = FALSE)
head(assay(rld))
```

```
##                 SRR1039508 SRR1039509 SRR1039510 SRR1039511 SRR1039512
## ENSG00000000003   9.553915   9.293960   9.634062   9.248913   9.651762
## ENSG00000000419   8.975079   9.189020   8.954938   9.161026   9.094500
## ENSG00000000457   8.020479   7.932336   7.940192   7.972483   7.897862
## ENSG00000000460   5.957300   6.012633   6.272144   5.999742   5.602442
## ENSG00000000971  11.859681  12.108263  12.019403  12.202862  12.354590
## ENSG00000001036  10.468316  10.302750  10.458210  10.185346  10.511347
##                 SRR1039513 SRR1039514 SRR1039515 SRR1039516 SRR1039517
## ENSG00000000003   9.463634   9.448013   9.797468   9.930933   9.667877
## ENSG00000000419   9.123299   9.146374   8.953586   9.046305   9.175449
## ENSG00000000457   7.978110   7.824048   7.874069   7.823220   7.955786
## ENSG00000000460   5.898212   5.922434   5.823606   6.055080   5.804083
## ENSG00000000971  12.539505  12.754418  12.558474  12.437680  12.766821
## ENSG00000001036  10.398871  10.300084  10.558848  10.315172  10.153987
##                 SRR1039518 SRR1039519 SRR1039520 SRR1039521 SRR1039522
## ENSG00000000003   9.800475  10.113804   9.790676   9.469947   9.840594
## ENSG00000000419   9.185251   9.127448   8.960671   9.119496   8.928404
## ENSG00000000457   7.937042   7.838148   8.005290   7.977851   7.808157
## ENSG00000000460   5.891935   6.108145   6.233475   6.034777   6.090024
## ENSG00000000971  12.877758  12.638696  12.411750  12.836135  12.616315
## ENSG00000001036  10.137456  10.341643  10.515925  10.286453  10.508171
##                 SRR1039523
## ENSG00000000003   9.544554
## ENSG00000000419   9.111313
## ENSG00000000457   7.871274
## ENSG00000000460   6.119933
## ENSG00000000971  12.819004
## ENSG00000001036  10.220749
```

We can illustrate the benefit of using the rlog over standard log transformation (+ a pseudo-count for genes with 0 counts where the log of 0 is infinity) by comparing the transformed values for two samples against each other. 

```r
par(mfrow=c(1,2))
plot(log2(cts[,1]+1), log2(cts[,2]+1), col = "cornflowerblue", xlab = "Sample 1", ylab = "Sample 2", main = "Log2 + 1")
plot(assay(rld)[,1], assay(rld)[,2], col = "indianred", xlab = "Sample 1", ylab = "Sample 2", main = "rlog")
```

<img src="day2-PART1_files/figure-html/unnamed-chunk-16-1.png" style="display: block; margin: auto;" />

We can use these transformed values to investigate how many features (genes) in our dataset exhibit variability across samples. This is useful to know as we only want to use variable features for PCA. Genes that don't explain any variation in the dataset aren't useful for helping us explore differences between the samples. 

```r
# calculate gene expression level variance between samples 
var <- rev(rowVars(assay(rld))[order(rowVars(assay(rld)))])

# plot variance for genes accross samples
plot(var, las = 1, main="Sample gene expression variance", xlab = "Gene", ylab = "Variance")
abline(v=1000, col="red") ; abline(v=500, col="green") ; abline(v=250, col="blue")
```

<img src="day2-PART1_files/figure-html/unnamed-chunk-17-1.png" style="display: block; margin: auto;" />

At around 500 the variance starts to spike upwards, so this is the number of variable features (genes) to use. Lets restrict the dataset to 500 genes for purposes of the PCA. Now lets extract and visualize the variance explained by each PC to determine which are most informative. 

```r
# modify variable feature number to be used in PCA and hierachical clutering based on no. of most variable features 
var_feature_n <- 500 

# perform PCA and order by variance 
rv <- rowVars(assay(rld))
select <- order(rv, decreasing = TRUE)[seq_len(min(var_feature_n, length(rv)))]
pca <- prcomp(t(assay(rld)[select, ]))

# extract the varioance explained by each PC 
percentVar <- pca$sdev^2/sum(pca$sdev^2)
names(percentVar)[1:5] <- c("PC1", "PC2", "PC3", "PC4", "PC5")
percentVar <- percentVar[1:5]

# plot variance for top 10 PCs 
barplot(percentVar[1:5], col = "indianred", las = 1, ylab = "% Variance", cex.lab = 1.2)
```

<img src="day2-PART1_files/figure-html/unnamed-chunk-18-1.png" style="display: block; margin: auto;" />
We can see that the majority of variance is explained by the first few PCs, therefore visualizing where samples fall along these PCs will be the most informative way to identify major differences between them, based on their gene expression profiles. 
Lets generate a PCA plot for PC1 vs PC2. 

```r
# construct data frame w/ PC loadings and add sample labels 
pca_df <- as.data.frame(pca$x)
pca_df$group <- dds@colData$group
pca_df$sample_ids <- colnames(dds)

# add colors for plotting to df 
pca_df$col <- NA
for(i in 1:length(levels(pca_df$group))){
  ind1 <- which(pca_df$group == levels(pca_df$group)[i])
  pca_df$col[ind1] <- i
}

# plot PC1 vs PC2
plot(pca_df[, 1], pca_df[, 2], 
     xlab = paste0("PC1 (", (round(percentVar[1], digits=3)*100), "% variance)"), 
     ylab = paste0("PC2 (", (round(percentVar[2], digits=3)*100), "% variance)"),
     main=paste0("PC1 vs PC2 for ", var_feature_n, " most variable genes"),
     pch=16, cex=1.35, cex.lab=1.3, cex.axis = 1.15, las=1, 
     panel.first = grid(),
     col=pca_df$col)
text((pca_df[, 2])~(pca_df[, 1]), labels = pca_df$group, cex=0.6, font=2, pos=4)
```

<img src="day2-PART1_files/figure-html/unnamed-chunk-19-1.png" style="display: block; margin: auto;" />

Looking at the plot, we can see some clear clustering by treatment group for some of the samples. For example, **untreated** samples appear to have much lower PC1 values than all the **Dex** samples, suggesting that some of the largest variability in gene expression differences between samples in this dataset explains differences between **untreated** and **Dex** treated samples. Therefore we expect the most substantial differential expression to be found for the **untreated** vs **Dex** analysis. 

The **Alb** treated samples, and the co-treated samples **(Alb_Dex)** do not seem to consistently cluster along PC1 and PC2. One explanation for this could be that treatment with **Alb** or co-treatment with **Alb + Dex** have incocnsistent effects on gene expression across study subjects. This is often the case with *in vivo* data sets, rather than those collected in *in vitro* systems such as cell culture. This could indicate that a greater number of replicates is required to accurately assess the impact of such treatment on gene expression, as there are usually a greater number of uncontrollable variables in *in vivo* datasets (e.g. genotype). 

Another explanation could be there is another variable, unrelated to sample groupings by treatment, that explains a lot of the variability in the dataset, such as a batch effect, that is driving the clustering observed in the PCA plot. For example, if the two Alb treated samples with high values for PC1 and PC2 were known to have been processed in a different batch to the Alb treated sample with lower values for these PCs, we could infer that there is a batch effect that is driving the differences between these samples. Unfortunately, no batch information was provided for this dataset, so we cannot confidnetly attribute the lack of clustering to batch effects. 

It is interesting that 4 of the samples, 1 from each treatment group, have much higher PC2 values than all other samples. This is a pattern we might expect to see if these 4 samples were processed in a separate batch, and the gene expression variation captured by PC2 is related to the batch effect, rather than our biological factor of interest (treatment group). For purposes of this example, lets create a fake variable for sample batch, and include this on the PCA plot. 


```r
pca_df$batch <- NA
pca_df$batch[pca_df$PC2 < 10] <- "Batch 1"
pca_df$batch[pca_df$PC2 > 10] <- "Batch 2"
pca_df$batch <- factor(pca_df$batch, levels = c("Batch 1", "Batch 2"))

# replot PC1 vs PC2 with the shape set to batch 
plot(pca_df[pca_df$batch=="Batch 1", 1], pca_df[pca_df$batch=="Batch 1", 2], 
     xlab = paste0("PC1 (", (round(percentVar[1], digits=3)*100), "% variance)"), 
     ylab = paste0("PC2 (", (round(percentVar[2], digits=3)*100), "% variance)"),
     main=paste0("PC1 vs PC2 for ", var_feature_n, " most variable genes"),
     pch=16, cex=1.35, cex.lab=1.3, cex.axis = 1.15, las=1, 
     panel.first = grid(),
     col=pca_df$col, 
     ylim=c(-14,20))
points(pca_df[pca_df$batch=="Batch 2", 1], pca_df[pca_df$batch=="Batch 2", 2], 
       col=pca_df$col, 
       pch=2, cex=1.35)
legend(9.5, 10.5, levels(pca_df$batch), pch = c(16, 2))
legend(1.5, 11.5, levels(pca_df$group), pch = 16, col = pca_df$col)
text((pca_df[, 2])~(pca_df[, 1]), labels = pca_df$group, cex=0.6, font=2, pos=4)
```

<img src="day2-PART1_files/figure-html/unnamed-chunk-20-1.png" style="display: block; margin: auto;" />

If this was a real batch variable, this would clearly indicate a batch effect. If batch effects are not accounted for in the differential expression analysis, these samples will increase the per gene variance (and therefore dispersion paracmeter in the negative-binomial model) and prevent the detection of genes that are truly differentially expressed between our experimental conditions of interest, effectively reducing statistical power. This could lead to both false-negatives and false-positives.

If you detect a batch effect in your data, you can try to: 
* use a statistical procedure to remove the batch effect from your data
* adjust for batch as a term in your statistical model when running differential expression 
* remove the samples driving the batch effect, if you have a particular reason to suspect these specific samples 

How you handle a batch effect is a complicated issue, and is largely dependent on the extent of the batch effect. If the batch effect is very large, it may be too diffcult to effectively remove it statistically, or regress out variation attributable to it in the DE analysis. This is where practicing your protcol and confirming you can get consistent results across replciates and batches comes in. If you do this work ahead of time, you reduce the risk of having to deal with this complicated batch effect in the analysis. 

If your experiment includes multiple batches, you should always include them in your unsupervised analyses to check for a batch effect. 

#### Hierarchical clustering

Hierarchical clustering is another complimentary approach to explore the relationships between your samples. While supervised clustering approaches exist, we will perform an unsupervised analysis so that we do not impose any restrictions on the clustering of the samples. 

Hierachical clustering is often associated with heatmaps, as it is a useful way to explore the results of hierachical clustering. Here we represent genes are rows, and individual samples as columns. The denrograms on the rows and the columns represent the 'distances' calculated between each of the genes/samples. Presenting the data in this way is useful as it allows us to identify samples whose patterns of gene expression are similar to each other, but also modules of genes that change in a similar way across our samples, and may share some common function of interest. 

The first step in a hierachical clustering analaysis is to *scale your data*. This means that expression levels are all transformed onto the same scale before clustering. This is important to do as we can only visualize so many colors at once, and a very highly expressed gene would mean that all the other genes would essentially invisible on this scale. Scaling for clustering in this way is typically performed by calculating Z-scores, where the mean for each gene across all the samples is subtracted from each of the individual expression values for each gene, this centers the expression values around 0. We then divide these values by the standard deviation to make sure the data is more tightly grouped, and we can represent lots of genes in the same scale. 

<center>
![](/Users/OwenW/Downloads/RNA-seq_workshop_July2020-master/figures/heatmaps.png)
</center>

Although we will not go into full detail here on how the actual clustering algorithm works to group samples and genes, once more StatQuest has an [excellent video](https://www.youtube.com/watch?v=oMtDyOn2TCc) on this topic. 

Similarly to the PCA, we perform the clustering using the rlog transformed data and the 500 most variable features, as features that do not vary across samples are not informative for dimension reduction appraoches. 


```r
# select top X no. of variable genes 
topVarGenes <- head(order(rowVars(assay(rld)), decreasing=TRUE), var_feature_n)

# set up gene expression matrix 
mat1 <- assay(rld)[topVarGenes,]

# scale matrix by each col. values 
mat_scaled = t(apply(mat1, 1, scale))

# set up colors for heatmap 
col = colorRamp2(c(-3, 0, 3), c("blue", "white", "red"))
cols1 <- brewer.pal(11, "Paired")
cols2 <- brewer.pal(9, "Greens")

# set up annotation bar for samples 
ha1 = HeatmapAnnotation(Group = colData(dds)$group, 
                        col = list(Group = c("untreated" = cols1[1], "Dex" = cols1[2], 
                                             "Alb" = cols1[5], "Alb_Dex" = cols1[6])), 
                                   show_legend = TRUE)

# se up column annotation labels (samples)
ha = columnAnnotation(x = anno_text(colData(dds)$SRR, 
                                    which="column", rot = 45, 
                                    gp = gpar(fontsize = 10)))

# generate heatmap object 
ht1 = Heatmap(mat_scaled, name = "Expression", col = col, 
              top_annotation = c(ha1), 
              bottom_annotation = c(ha),
              show_row_names = FALSE)

# plot the heatmap 
draw(ht1, row_title = "Genes", column_title = "Top 500 most variable genes")
```

<img src="day2-PART1_files/figure-html/unnamed-chunk-21-1.png" style="display: block; margin: auto;" />

As we saw in the PCA, the Alb and co-treated samples do not form any clear clusters. We may want to remove them and perform the clustering again so that we can compare the untreated and Dex samples more easily. 


```r
ind_to_keep <- c(which(colData(rld)$group=="untreated"), which(colData(rld)$group=="Dex"))
topVarGenes <- head(order(rowVars(assay(rld)[,ind_to_keep]), decreasing=TRUE), var_feature_n)

# set up gene expression matrix 
mat1 <- assay(rld)[topVarGenes, ind_to_keep]

# scale matrix by each col. values 
mat_scaled = t(apply(mat1, 1, scale))

# set up colors for heatmap 
col = colorRamp2(c(-3, 0, 3), c("blue", "white", "red"))
cols1 <- brewer.pal(11, "Paired")
cols2 <- brewer.pal(9, "Greens")

# subset coldata for samples in untx and ex groups
colData_sub <- colData(dds)[ind_to_keep, ]

# set up annotation bar for samples 
ha1 = HeatmapAnnotation(Group = colData_sub$group, 
                        col = list(Group = c("untreated" = cols1[1], "Dex" = cols1[2])), 
                                   show_legend = TRUE)

# se up column annotation labels (samples)
ha = columnAnnotation(x = anno_text(colData_sub$SRR, 
                                    which="column", rot = 45, 
                                    gp = gpar(fontsize = 10)))

# generate heatmap object 
ht1 = Heatmap(mat_scaled, name = "Expression", col = col, 
              top_annotation = c(ha1), 
              bottom_annotation = c(ha),
              show_row_names = FALSE)

# plot the heatmap 
draw(ht1, row_title = "Genes", column_title = "Top 500 most variable genes")
```

<img src="day2-PART1_files/figure-html/unnamed-chunk-22-1.png" style="display: block; margin: auto;" />

There does indeed seem to be relatively good clustering between untreated and Dex samples, suggesting there are unique gene expression programs defining the Dex samples from the untreated. However, one of these samples in the Dex group seems to be clustered further away from the other Dex samples. This could be the Dex treated sample that clustered away from the other Dex treated samples on the PCA. We can add an annotation bar for the fake batch effect we created earlier to this plot to confirm this. 


```r
colData_sub$batch <- "Batch 1"
colData_sub$batch[colData_sub$SRR=="SRR1039516"] <- "Batch 2"
colData_sub$batch[colData_sub$SRR=="SRR1039517"] <- "Batch 2"

# set up annotation bar for samples 
ha1 = HeatmapAnnotation(group = c(as.character(colData_sub$group)), 
                        batch = c(as.character(colData_sub$batch)),
                        col = list(group = c("untreated" = cols1[1], "Dex" = cols1[2]),
                                   batch = c("Batch 1" = cols1[5], "Batch 2" = cols1[6])), 
                                   show_legend = TRUE)

# generate heatmap object 
ht2 = Heatmap(mat_scaled, name = "Expression", col = col, 
              top_annotation = c(ha1), 
              bottom_annotation = c(ha),
              show_row_names = FALSE)

# plot the heatmap 
draw(ht2, row_title = "Genes", column_title = "Top 500 most variable genes")
```

<img src="day2-PART1_files/figure-html/unnamed-chunk-23-1.png" style="display: block; margin: auto;" />

Based on our newly labeled plot it does seem that these 2 samples are outliers based on the hierachical clustering, which would support the presence of a batch effect if these data were infact collected in multiple batches. Again without the metadata to indicate this is a batch effect using a method to correct for a batch effect is inappropriate. In this case I would reach out to the seqeuncing center or authors of the paper where the data was published and ask for additional metadata to confirm our suspicion about batch effects in this data.

We may also wish to cluster based on the pairwise correlations between gene expression samples, another complimenatry aproach to help us understand how patterns of gene expression in the samples are rleated to each other. 


```r
rld_cor <- cor(assay(rld)[topVarGenes[1:50], ind_to_keep])   

cols = colorRamp2(c(0, 1), c("yellow3", "royalblue2"), transparency = 0.1)
cols1 <- brewer.pal(11, "Paired")

# make heatmap annotations 
ha1 = HeatmapAnnotation(group = as.factor(colData_sub$group), 
                        show_legend = FALSE, 
                        show_annotation_name = FALSE,
                        #text = anno_text(colData_sub$SRR, location = unit(1, "npc"), rot = 50, just = "right"),
                        col = list(group = c("untreated" = cols1[5], "Dex" = cols1[6])))

# generate heatmap object 
ht_sub1 = Heatmap(rld_cor, name = "Perason correlation", col = cols, 
                  rect_gp = gpar(col= "white"),
                  top_annotation = ha1, row_labels = rep("", 8))

# plot it 
draw(ht_sub1, column_title = "Correlation for top variable genes")
```

<img src="day2-PART1_files/figure-html/unnamed-chunk-24-1.png" style="display: block; margin: auto;" />

That wraps up exploratory analysis. In the next R markdown, we will perform the differential expression analysis. 

## Session Information

```r
sessionInfo()
```

```
## R version 3.6.1 (2019-07-05)
## Platform: x86_64-apple-darwin15.6.0 (64-bit)
## Running under: macOS Catalina 10.15.5
## 
## Matrix products: default
## BLAS:   /Library/Frameworks/R.framework/Versions/3.6/Resources/lib/libRblas.0.dylib
## LAPACK: /Library/Frameworks/R.framework/Versions/3.6/Resources/lib/libRlapack.dylib
## 
## locale:
## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
## 
## attached base packages:
##  [1] grid      parallel  stats4    stats     graphics  grDevices utils    
##  [8] datasets  methods   base     
## 
## other attached packages:
##  [1] kableExtra_1.1.0            xtable_1.8-4               
##  [3] circlize_0.4.8              ComplexHeatmap_2.2.0       
##  [5] RColorBrewer_1.1-2          gplots_3.0.3               
##  [7] pheatmap_1.0.12             vsn_3.54.0                 
##  [9] biomaRt_2.42.0              DESeq2_1.26.0              
## [11] SummarizedExperiment_1.16.1 DelayedArray_0.12.2        
## [13] BiocParallel_1.20.1         matrixStats_0.55.0         
## [15] Biobase_2.46.0              GenomicRanges_1.38.0       
## [17] GenomeInfoDb_1.22.0         IRanges_2.20.2             
## [19] S4Vectors_0.24.3            BiocGenerics_0.32.0        
## [21] tximport_1.14.0             ggplot2_3.3.0              
## [23] dplyr_0.8.5                
## 
## loaded via a namespace (and not attached):
##   [1] colorspace_1.4-1       rjson_0.2.20           htmlTable_1.13.3      
##   [4] XVector_0.26.0         GlobalOptions_0.1.1    base64enc_0.1-3       
##   [7] clue_0.3-57            rstudioapi_0.11        farver_2.0.3          
##  [10] affyio_1.56.0          bit64_0.9-7            AnnotationDbi_1.48.0  
##  [13] xml2_1.2.2             splines_3.6.1          geneplotter_1.64.0    
##  [16] knitr_1.28             Formula_1.2-3          annotate_1.64.0       
##  [19] cluster_2.1.0          dbplyr_1.4.2           png_0.1-7             
##  [22] readr_1.3.1            BiocManager_1.30.10    compiler_3.6.1        
##  [25] httr_1.4.1             backports_1.1.5        assertthat_0.2.1      
##  [28] Matrix_1.2-18          limma_3.42.2           acepack_1.4.1         
##  [31] htmltools_0.4.0        prettyunits_1.1.1      tools_3.6.1           
##  [34] gtable_0.3.0           glue_1.3.1             GenomeInfoDbData_1.2.2
##  [37] affy_1.64.0            rappdirs_0.3.1         Rcpp_1.0.3            
##  [40] vctrs_0.2.3            gdata_2.18.0           preprocessCore_1.48.0 
##  [43] xfun_0.12              stringr_1.4.0          rvest_0.3.5           
##  [46] lifecycle_0.2.0        gtools_3.8.1           XML_3.99-0.3          
##  [49] zlibbioc_1.32.0        scales_1.1.0           hms_0.5.3             
##  [52] yaml_2.2.1             curl_4.3               memoise_1.1.0         
##  [55] gridExtra_2.3          rpart_4.1-15           latticeExtra_0.6-29   
##  [58] stringi_1.4.6          RSQLite_2.2.0          genefilter_1.68.0     
##  [61] checkmate_2.0.0        caTools_1.18.0         shape_1.4.4           
##  [64] rlang_0.4.5            pkgconfig_2.0.3        bitops_1.0-6          
##  [67] evaluate_0.14          lattice_0.20-40        purrr_0.3.3           
##  [70] labeling_0.3           htmlwidgets_1.5.1      bit_1.1-15.2          
##  [73] tidyselect_1.0.0       magrittr_1.5           R6_2.4.1              
##  [76] Hmisc_4.3-1            DBI_1.1.0              pillar_1.4.3          
##  [79] foreign_0.8-76         withr_2.1.2            survival_3.1-11       
##  [82] RCurl_1.98-1.1         nnet_7.3-13            tibble_2.1.3          
##  [85] crayon_1.3.4           KernSmooth_2.23-16     BiocFileCache_1.10.2  
##  [88] rmarkdown_2.1          jpeg_0.1-8.1           GetoptLong_0.1.8      
##  [91] progress_1.2.2         locfit_1.5-9.1         data.table_1.12.8     
##  [94] blob_1.2.1             webshot_0.5.2          digest_0.6.25         
##  [97] openssl_1.4.1          munsell_0.5.0          viridisLite_0.3.0     
## [100] askpass_1.1
```



