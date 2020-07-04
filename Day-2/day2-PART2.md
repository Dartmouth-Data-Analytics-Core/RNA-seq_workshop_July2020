---
title: "Day 2, PART-2"
author: "Prepared by the Data Analytics Core (Center for Quantitative Biology at Dartmouth)"
output:
    html_document:
      keep_md: TRUE
      theme: default
      number_sections: TRUE
---
*Bulk RNA-seq data analysis data analysis workshop, July 2020*

## Differential expression analysis in R 

### Introduction to DEG in R

After exploring our data, we are ready to run the differential expression analysis. As our exploratory analysis showed the Alb treated and co-treated samples did not cluster together, so going forward we will focus on the comparison between untreated samples and Dex treated samples. We will continue using *DESeq2* to perform the analysis. 

*NOTE:* You must change the below line, and all other lines loading images, to the directory on your computer!!

</center>
![Overview](/Users/OwenW/Downloads/RNA-seq_workshop_July2020-master/figures/overview.png)
</center>

**In this script, we will:**
1. Run the differential expression analysis comparing experimental groups of interest  
2. Annotate the results with information on each gene  
3. Explore our results:  
  * How many significantly DE genes did you detect?  
  * Which genes were the most confident DEGs identified?  
  * What is the overall distribution of P-values, fold-change, and expression levels among our results?  
4. Correct for multiple hypothesis testing (and discuss why we bother)  
5. Hierachical clustering of samples based on significant DEGs  

Set the root directory for the whole markdown  

```r
knitr::opts_knit$set(root.dir = '/Users/OwenW/Downloads/RNA-seq_workshop_July2020-master/')
```

Lets start by loading the required libraries again. 

```r
library(tximport)
library(DESeq2)
library(biomaRt)
library(vsn)
library(dplyr)
library(pheatmap)
library(gplots)
library(RColorBrewer)
library(ComplexHeatmap)
library(readr)
library(circlize)
library(EnhancedVolcano)
library(apeglm)
library(xtable)
library(kableExtra)
```

Read in the DESeq2 dataset we created in PART-1, which contains the raw counts, normalization factors, and sample metadata. 

```r
load("Day-2/DESeq2.rdata")
```

***

### Apply the DESeq2 procedure to the data 

Now we apply the DEseq() function to the dataset to perform the analysis. This function is the main powerhouse of the DESeq2 package and does **a lot** under the hood. It is important you understand the general principles of this analysis of before running your own analysis. 

**At a high level, the major function performed by `DESeq2` are:**
* estimation of size factors (`estimateSizeFactors()`)
* estimation of dispersion (`estimateDispersions`)
* fitting of the negative bionmial generalized linear model (GLM) and wald statistics for differential expression testing (`nbinomWaldTest`)

Lets run `DESeq2` on our dataset:

```r
# run the DEseq2 analysis 
dds <- DESeq(dds)
```

```
## estimating size factors
```

```
## estimating dispersions
```

```
## gene-wise dispersion estimates
```

```
## mean-dispersion relationship
```

```
## final dispersion estimates
```

```
## fitting model and testing
```

Before running the differential expression analysis, lets have a look at some of the standard characteristics of RNA-seq data. The first and most obvious thing to do is look at how the distribution of the raw counts. 

```r
hist(counts(dds, normalized=FALSE)[,5], breaks = 500, col="blue",
     xlab="Raw expression counts", ylab="Number of genes",
     main = "Count distribution for sample X")
```

<img src="day2-PART2_files/figure-html/unnamed-chunk-4-1.png" style="display: block; margin: auto;" />

Perhaps the most obvious feature of this distribution is the large number of genes with very low count values. This occurs as there are many genes expressed at low levels relative to the highly expressed genes, which are fewer in number. This causes the distribution to have a long right tail, ultimately caising the dynamic range of RNA-seq data to be very large. 

These features of how RNA-seq data is distributed are important in selecting the statistical model used to test differential expression. Importantly, we can see from the histogram that the data is **not** normally distributed, therefore any statistical model based on the normal distribution is not appropriate for this dataset. By looking again at the matrix of raw counts, it is actually clear that RNA-seq is integer count data, therefore we should use a statistical model for count-based data. 


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
At this point it might be useful to define a few terms that are really important to know in order to understand how these statistical models are applied. Each of these terms are calculated from the count matrix by DESeq2 and are used to fit the appropriate statistical model to the data, and compare counts within and between samples.

*mean - the average count of a gene across samples
*variance - the spread of count values across samples for a gene
*dispersion - the amount that the variance deviates from the mean

Usually, the most appropriate statistical model for DE analysis of RNA-seq data is a variation of the binomial distribution (e.g. probability of getting 10 heads when you toss a coin 20 times), called the negative binomial (NB) distribution. One feature of RNA-seq data in particular makes the NB distribution a good fit for such data: overdispertion, i.e. the variance is usually greater than the mean for RNA-seq data. 

We can visualize overdispersion in RNA-seq data by plotting the mean-variance relationship for a group of replicates in our data. 

```r
# calculate mean and varaince for group of replicates
mean_counts <- apply(counts(dds, normalized=FALSE)[,1:3], 1, mean)
variance_counts <- apply(counts(dds, normalized=FALSE)[,1:3], 1, var)

# plot the mean variance trend 
plot(log10(mean_counts), log10(variance_counts), 
     ylim=c(0,9), xlim=c(0,9), 
     ylab = "log10 (mean counts)", xlab = "log10 (varaince)", 
     main = "Mean-variance trend", las = 1)

# add line for x=y
abline(0,1,lwd=2,col="red")
```

![](day2-PART2_files/figure-html/unnamed-chunk-6-1.png)<!-- -->

We can clearly see a few features of the mean variance trend from this plot:  
1. The data does not fall along the x = y line, as it would if the mean = varaince. Instead, the varaince is generally greater than the mean, making the varaince overdispersed. 
2. There is more difference in the varaince between low count genes than there is amongst higher count genes, therefore the varaince is unequal across the range of count values. This property of the data is referred to as heteroscedasticity. 

The negative binomial distribution is a good fit for count data that is overdispersed in this way, as it includes a **dispersion parameter** that is able to account for the relationship between the mean and the varaince, which is clearly important in these data as the varaince changes dramatically depending on the expression level of the gene you are observing.  

We can plot a few different NB distributions to examine how the dispersion parameter affects the spread of the data.

```r
# generate a random varaible using the negative binomial distribution
### dispersion = 10
par(mfrow=c(3,1))
hist(rnbinom(n = 10000, mu = 100, size = 1/0.001), 
     xlim = c(0, 300), xlab = "", breaks = 500, 
     main = " Dispersion 0.001")
### dispersion = 10
hist(rnbinom(n = 10000, mu = 100, size = 1/0.01), 
     xlim = c(0, 300), xlab = "", breaks = 500, 
     main = " Dispersion 0.01")
### dispersion = 10
hist(rnbinom(n = 10000, mu = 100, size = 1/0.1), 
     xlim = c(0, 300), xlab = "", breaks = 500, 
     main = " Dispersion 0.1")
```

<img src="day2-PART2_files/figure-html/unnamed-chunk-7-1.png" style="display: block; margin: auto;" />

Note: The above example for plotting NB distributions at various disperions was adapted from the *Data Analysis for the Life Sciences series* available on edX and at [rafalab](https://rafalab.github.io/pages/harvardx.html), and is an excellent resource to learn more about how we model RNA-seq data for differential expression analysis. 

It is clear that as the disperion increases, the varaition around the mean also increases. The mean, variance, and dispersion are linked by the equation: 

variance = mean + dispersion x 2 mean-squared ( var = mu + disp. * mu^2 )

In order to accurately model differential expression for the genes in our dataset, DESeq2 uses this equation to obtain estimates for the dispersion of each gene within each sample group (e.g. Control and Dex separately). **However,** for the small number of replicates avaiable in RNA-seq data, these estimates of dispersion at the gene-level are often inaccurate (yet another reason to use more replicates..). 

To improve these gene-level estimates of dispersion, DESeq2 uses another statistical model called empirical bayes to 'shrink' these inital dispersion estimates toward a *'prior'* mean, which is calculated by fitting a curve to the inital dispersion estimates. This procedure is able to produce more accuarate estimates of disperion as it shares information across genes with similar expression levels to predict a more approriate dispersion for those genes. This is rational as the formula linking the mean, variance, and dispersion tells us that the variance is the only thing affecting the magnitude of the dispersion for genes with the similar mean expression. 

The major factors affecting how much a gene's dispersion is shrunk toward the prior mean are:
1. the number of samples in the group under consideration (use more replicates!)
2. how far the inital dispersion is from the prior mean

<center>
![DEseq2 dispersion estimation](/Users/OwenW/Downloads/RNA-seq_workshop_July2020-master/figures/dispersion_estimation.png)
</center>

This Figure taken frm the DESeq2 paper demonstrates the process of shrinkage, where the inital dispersion esimates for each gene (estimated by maximum-likelihood) are shrunken towards the prior mean (based on the fitted curve in red) to a final MAP estimate. For dispersion estimates further away from the line, you can see that their estimates are shrunken more than those are are originally closer to the line. 

This curve also shows the expected trend for dispersion estimates over a range of expression levels. Importantly, the dispersion tends to decrease as the mean increases. Through inspecting disperion estimates for our own data, we can determine if the model is a good fit for our data, and therefore if it can be used to accurately test DE.   

**We can plot the dispersion estimates for our own data using:**

```r
plotDispEsts(dds)
```

<img src="day2-PART2_files/figure-html/unnamed-chunk-8-1.png" style="display: block; margin: auto;" />

This is an example of a well calibrated set of dispersion estimates due to these two features: the final MAP estimates are well scattered around the fitted line, and the dispersion trend decreases with increasing mean expression. 

If the MAP estimates were more structured in these plots, we would be concerned that the model is not estimating dispersions well for our data, indicating something may be wrong with the dataset, e.g. outlier samples, a batch effect, low quality samples/data, potential contamination etc. 

**It is important to confirm your dispersion estimates are well calibrated before performing your differential expression analysis, as accurate estimation of dispersion is critical in controlling the false-positive rate in experiments with smaller sample sizes (most RNA-seq experiments)**. 

***

### Differential expression analysis - Hypothesis testing

Now that we understand how the dispersions are estimated, we are ready to fit the data and test each gene for differential expression. We fit the data using a generalized linear model (GLM) of the negative binomial family. GLM's are a family of statistical models that generalize linear regression such that response variables can be modeled by distributions other than the normal distribution (which we know isn't approrpiate for RNA-seq data). Essentially, the GLM allows us to model a *response variable* (in this case, the experimental group) as a linear combination of a set of variables, gene expression across samples. 

</center>
![](/Users/OwenW/Downloads/RNA-seq_workshop_July2020-master/figures/neg-binom.png)
</center>

In order to fit the GLM, we need the mean count of each gene across the samples in each experimental group, and the dispersion of that gene in those groups. The mean count is a combination of the expected expression level and the size factor, so that our model is corrected for library size and composition. 

The process of fitting the model to the expression and dispersion values for each gene results in final set of *model coefficients* for each sample group, which can be interpreted as the log2 fold-change in expression for that gene between the baseline group and each comparison group. 

Each of the model coefficients has an associated *standard error* associated with it, which we can use to calculate a P-value and perform a process called hypothesis testing. Through hypothesis testing we test the *null hypothesis* that the log2 fold-change between experimental groups for an individual gene is not significnatly different from 0 (no change in expression). The default test used by DESeq2 for hypothesis testing is the *Wald-test*, which is implemented as follows:  

1. The *coefficient (log 2 fold-change)* is divided by the *standard error* (measure of statistical accuracy of the measurement). 
2. The resulting *Z-statistic* is compared to a standard normal distribution (mean = 0, sd = 1) in order to compute a P-value. 
3. If the P-value is less than our pre-determined threshold for significance, we reject the null hypothesis and accept the alternative, that the gene is significantly DE. 

**Note:** DESeq2 can also implement a *likelihood ratio test* (LRT), which is used to compare expression accross more than two groups. For example, if you collected samples over a range of time points and you wanted to test if gene expression changed significantly over these time points, you could use the LRT instead of the wald-test. 

`DESeq2` already performed all of the steps for hypothesis testing using the wald-test for us when we ran the `DESeq2()` function. All we have to do is tell DESeq2 which results we want to look at, which can be done using the `results()` function, and specifying the coefficients that we want by using the `names` agument. 

```r
# quickly check the available coefficients we could extract 
resultsNames(dds)
```

```
## [1] "Intercept"                  "group_Dex_vs_untreated"    
## [3] "group_Alb_vs_untreated"     "group_Alb_Dex_vs_untreated"
```

```r
# get results for DEG analysis (and order by Pval) by specifying design 
res <- results(dds, 
  name = "group_Dex_vs_untreated", 
  alpha = 0.05, 
  lfcThreshold = 0)
```

A couple of things to note here:  
- `alpha` is set to 0.05 (5%) to correct the P-values for multiple hypothesis testing (example with more detail on this coming up below). By default, the "BH" method is used (Benjamini & Hochberg) which controls the false discovery rate (FDR). Corrected P-values are found in the `padj` column of the `results()` output, while the uncorrected P-values are found in the `pvalue` column. Other methods to control for multiple hypothesis testing can be specified using the `pAdjustMethod` argument in the `results()` function, such as the more conservative **Bonferonni** method. 

- `lfcThreshold` is set to 0, and is the default value. This tests the hypothesis that the log2 fold change values between our experimental conditions are equal to 0. Different fold change values can be specified, which can be useful if you observe a large number of significantly differentially expressed genes with small fold changes, and you want to restrict the test to the genes with the largest differences (fold changes) between your conditions (we could also achieve this by restricting the results to genes with significant P-values AND have an absolute fold change > a specific threshold, however when we do this, the P-values loose some of their meaning).  

#### Note on P-values: 
P-value thresholds do not need to be set at 0.05 for every experiment. You can be more or less stringent than this dependning on the nature of your experiment: if you want to be very conservative and restrict your results to few results that are likely to be true positives, you may wish to restrict the results to a more stringent threshold. If your experiment is very preliminary and you care less about capturing some false positives than missing true positives, you may wish to relax your threshold.  

**Additional note:** To extract the results, we could also use the `contrast` argument in a similar way to how we used the `names` agument. The first group specified to `contrast` is used as the numerator in calculating the fold change, and the second group is used as the denominator, therefore the second group is used as the baseline for the comparison. 

```r
res <- results(dds, alpha = 0.05, 
  contrast = c("group", "Dex", "untreated"), 
  lfcThreshold = 0)
```

This is useful when we have multiple levels in the experimental design variable and we wish to extract coefficients for the results from testing specific levels against one another. It is generally the same as using the `names` argument to extract coefficients, with some exceptions that are discussed in the DESeq2 documentation. 

Lets have a quick look at the results and how many genes were statistically significant at an adjusted P-value threshold of 0.05. 

```r
# order by adj Pval 
res_ord <- res[order(res$padj),] 

# quick check for how many DEGs with significance @ 5% level in either FC direction 
sum(res$padj < 0.05, na.rm=TRUE)
```

```
## [1] 1673
```

```r
sum(res$padj < 0.05 & res$log2FoldChange>2, na.rm=TRUE)
```

```
## [1] 102
```

```r
sum(res$padj < 0.05 & res$log2FoldChange < -2, na.rm=TRUE)
```

```
## [1] 71
```

You may have noticed I am using `na.rm=TRUE` in the `sum()` function above. Why might this be? 

```r
table(is.na(res$padj))
```

```
## 
## FALSE  TRUE 
## 14918  9501
```

This is not a mistake, but rather part of a deliberate filtering process conducted by DESeq2, in order to flag genes that have little or no change of being differentially expressed. This is of value as it means we can correct for fewer total tests and increase our statistical power to identify true positives.The three ways which DESeq2 filters results are:   
- Genes with counts = 0 in all samples
- Genes with extreme outliers (determined using Cook's distance)
- *Independent filtering* (identifying genes with low counts)

*Independent filtering*, DESeq2 carries out an iterative process where it maximizes the value of the number of rejections over the quantiles of the mean normalized counts. Once the maximum number of rejections is identified, DESeq2 will select the quantile of the normalized counts that is 1 standard deviation below this maximum, and filter any results with mean counts below this threshold. It is essentially a fancy (and cool) way of reducing the number of tests we need to run. 

We can plot the number of rejections of the null hypotesis against mean counts, along with a vertical line, to help us understand at which mean count value DESeq2 chose to filter results for. Any genes with a mean expression value below this line will have their `padj` values set to NA, and discarded during multiple testing correction. 

```r
plot(metadata(res_ord)$filterNumRej, 
     type="b", ylab="number of rejections",
     xlab="quantiles of filter (mean norm. counts)")
lines(metadata(res_ord)$lo.fit, col="red")
abline(v=metadata(res_ord)$filterTheta)
```

<img src="day2-PART2_files/figure-html/unnamed-chunk-13-1.png" style="display: block; margin: auto;" />

Its worth removing these results with NAs before moving forward to make our lives a little easier when handling the adjusted P-values. 

```r
res_ord <- res_ord[!is.na(res_ord$padj),]
```

### Add gene annotation to the results  

We want to also add the annotation data for each gene (symbol, genome coordinates, etc.) to the results. Since we used Ensembl version 97 to annotate these data, we need to use the Ensembl 97 annotation data to annotate these results. We can obtain this for our species of interest in a flat file format using the [BioMart on the Ensembl website](http://uswest.ensembl.org/biomart/martview/b0399bb192186dea3aedf87d82a4580c). 

```r
# read in the flat file we downloaded and have a look at it 
anno <- read.delim("Day-2/GRCh38.p12_ensembl-97.txt", stringsAsFactors = T, header = T)
anno <- anno[order(anno$Chromosome.scaffold.name),]
dim(anno)
```

```
## [1] 66832    10
```

```r
# have a look at the first few rows 
head(anno)
```

```
##     Gene.stable.ID Gene.stable.ID.version Chromosome.scaffold.name
## 47 ENSG00000201900      ENSG00000201900.1                        1
## 60 ENSG00000243173      ENSG00000243173.3                        1
## 66 ENSG00000276103      ENSG00000276103.1                        1
## 67 ENSG00000276470      ENSG00000276470.1                        1
## 75 ENSG00000212257      ENSG00000212257.1                        1
## 84 ENSG00000275823      ENSG00000275823.1                        1
##    Gene.start..bp.
## 47       114727720
## 60       162777730
## 66       150608507
## 67        11843812
## 75        65022968
## 84        63549294
##                                                              Gene.description
## 47                     RNY1 pseudogene 13 [Source:HGNC Symbol;Acc:HGNC:50876]
## 60  RNA, 7SL, cytoplasmic 861, pseudogene [Source:HGNC Symbol;Acc:HGNC:46877]
## 66                                                                           
## 67                                                                           
## 75 RNA, U6 small nuclear 1176, pseudogene [Source:HGNC Symbol;Acc:HGNC:48139]
## 84                                                                           
##    Gene.end..bp. Strand Transcript.count Gene.type  Gene.name
## 47     114727824      1                1  misc_RNA    RNY1P13
## 60     162778014      1                1  misc_RNA  RN7SL861P
## 66     150608623     -1                1     snRNA    RF00015
## 67      11843984      1                1  misc_RNA    RF02156
## 75      65023074     -1                1     snRNA RNU6-1176P
## 84      63549366     -1                1  misc_RNA    RF02107
```

Lets have a look at the Chromosome distribution of features

```r
tab1 <- table(anno$Chromosome.scaffold.name)
tab1[1:22]
```

```
## 
##    1   10   11   12   13   14   15   16   17   18   19    2   20   21   22    3 
## 5471 2332 3360 3054 1397 2282 2221 2556 3060 1242 2992 4196 1457  872 1384 3185 
##    4    5    6    7    8    9 
## 2651 2983 3059 3014 2482 2327
```

Lets also quickly check that nothing is duplicated in the ENSG ID column of our annotation, as this would cause problems when merging with our results. 

```r
any(duplicated(anno$Gene.stable.ID))
```

```
## [1] FALSE
```

Now lets add the annotation for each gene name directly to the results. 

```r
# use match() to find corresponding indicies (rows) for each ENSG ID 
mat1 <- match(rownames(res_ord), anno$Gene.stable.ID)
table(is.na(mat1))
```

```
## 
## FALSE 
## 14918
```

```r
# add gene names to results as a new column 
res_ord$gene <- as.character(anno$Gene.name[mat1])
head(res_ord, 20)
```

```
## log2 fold change (MLE): group Dex vs untreated 
## Wald test p-value: group Dex vs untreated 
## DataFrame with 20 rows and 7 columns
##                         baseMean     log2FoldChange             lfcSE
##                        <numeric>          <numeric>         <numeric>
## ENSG00000179094 879.162412923049    3.1955116536459 0.520089631550247
## ENSG00000109906 427.519588336145     7.102329380743  1.26387280388192
## ENSG00000157613 2837.46444827888 -0.977869256859263 0.188367598351614
## ENSG00000179593 74.1106067372775   9.60746204046761  1.85298907731736
## ENSG00000182010 104.205807294111  -1.89767754411958 0.361015967801541
## ...                          ...                ...               ...
## ENSG00000183160 3003.58826886215  -1.17316952492041 0.242135177240336
## ENSG00000185950 2099.90917446653   2.16736650931075 0.439155979458348
## ENSG00000198517 598.229869113694   0.70436711786436 0.145082798393017
## ENSG00000124151 1605.14338263049   1.39771472406699 0.290457449620042
## ENSG00000152583 1214.50066415466    4.5976490815425 0.961457272963187
##                              stat               pvalue                 padj
##                         <numeric>            <numeric>            <numeric>
## ENSG00000179094  6.14415566047903 8.03899782915565e-10 1.19925769615344e-05
## ENSG00000109906  5.61949696118832 1.91514247454437e-08 0.000142850477176265
## ENSG00000157613  -5.1912816504351 2.08851413522592e-07  0.00064503652440172
## ENSG00000179593  5.18484547916316 2.16194035528127e-07  0.00064503652440172
## ENSG00000182010 -5.25649199307101 1.46829166165317e-07  0.00064503652440172
## ...                           ...                  ...                  ...
## ENSG00000183160 -4.84510155976203 1.26546929480719e-06  0.00109369064101773
## ENSG00000185950  4.93530000885781 8.00276479495294e-07  0.00109369064101773
## ENSG00000198517  4.85493198136619 1.20428072665044e-06  0.00109369064101773
## ENSG00000124151  4.81211525438713  1.4934118349848e-06  0.00117256409233175
## ENSG00000152583  4.78195881484433 1.73595250807773e-06  0.00129484697577518
##                        gene
##                 <character>
## ENSG00000179094  AC129492.1
## ENSG00000109906      ZBTB16
## ENSG00000157613     CREB3L1
## ENSG00000179593     ALOX15B
## ENSG00000182010       RTKN2
## ...                     ...
## ENSG00000183160     TMEM119
## ENSG00000185950        IRS2
## ENSG00000198517        MAFK
## ENSG00000124151       NCOA3
## ENSG00000152583     SPARCL1
```

Lets also add some other columns that might be of interest to us when reviewing the results. 

```r
res_ord$chr <- as.character(anno$Chromosome.scaffold.name[mat1])
res_ord$start <- as.character(anno$Gene.start..bp.[mat1])
res_ord$end <- as.character(anno$Gene.end..bp.[mat1])
res_ord$strand <- as.character(anno$Strand[mat1])
```

***

### Visualization of Differential Expression 

#### Volcano plot

Volcano plots are a useful visualization for exploring your results, the log2 fold change (x-axis) is plotted against the -log10 P-value. Since the -log10() of a really small number is a very large value, any gene that has a very small P-value and was significantly differentially expressed, will appear higher up along the y-axis. In contrast, the -log10 of 1 (`-log10(1)`) is equal to `0`, therefore genes with low statistical significance (P-values approaching 1) will appear lower down on the y-axis. 

Similarly, genes with larger fold changes will appear further along the x-axis, in both directions. Genes with a positive fold change represent genes whose expression was greater than the group of the experimental design variable used as baseline, while genes with a negative fold change represent genes whose expression was lower than in the baseline group. The fold-change value of genes with non-significant fold changes is not meaningful, as there is not enough statistical confidence in these fold changes. 

```r
plot(res$log2FoldChange, -log10(res$pvalue), 
     main = "Volcano plot", 
     las = 1, col = "indianred",
     ylab = "- log10 P-value", xlab = "log2 Fold change")

# add horizontal lines to help guide interpretation
abline(h=-log10(0.05/nrow(res)), lty = 2, col = "black") # Bonferonni 
abline(h=-log10(0.05), lty = 2, col = "black") # nominal P-value 
```

<img src="day2-PART2_files/figure-html/unnamed-chunk-20-1.png" style="display: block; margin: auto;" />
  
Here we can clearly see that there are quite a few genes above our significance threshold in both the up and downregulation directions (+ve and -ve fold changes), that also have absolute log2 fold change values of at least 2 or more. Of particular interest, there seem to be a few genes with very large fold change values & -log10 P-values, making them especially interesting as their effect size is large AND our confidence in this fold change is good. 

It is a little hard to make specific inferences from this plot at the individual gene level, so some labels for interesting data points ( and some colors) would definitely improve this volcano plot, and make it more informative. We will use the **ggpolot2** R package to do this, and we will color each point based on a combination of fold change and P-value, as these determine which genes are of most interest to us. 

```r
# save a dataframe from the results() output
res_tmp <- as.data.frame(res_ord)

# add a column that will be used to save the colors we want to plot 
res_tmp$cols <- c()

# set the significance cut off (alpha) and fold change threshold to be used for coloring of genes 
alpha <- 0.05/nrow(res)
fc_cutoff <- 2

# loop through our dataframe and add values to the color column based on magnitude of alpha and LFCs 
res_tmp$cols <- NA
for(i in 1:nrow(res_tmp)){
    if(is.na(res_tmp$pvalue[i])){
      res_tmp$cols[i] <- NA
    }
    else if(res_tmp$pvalue[i]<=alpha & res_tmp$log2FoldChange[i] > fc_cutoff){
      res_tmp$cols[i] <- "indianred"
    } 
    else if(res_tmp$pvalue[i]<=alpha & res_tmp$log2FoldChange[i] < -fc_cutoff){
      res_tmp$cols[i] <- "indianred"
    } 
    else if(res_tmp$pvalue[i]<=alpha & res_tmp$log2FoldChange[i]>-fc_cutoff & res_tmp$log2FoldChange[i]<fc_cutoff){
      res_tmp$cols[i] <- "cornflowerblue"
    } 
    else if(res_tmp$pvalue[i]>alpha & res_tmp$log2FoldChange[i] > fc_cutoff){
      res_tmp$cols[i] <- "gray47" 
    }
    else if(res_tmp$pvalue[i]>alpha & res_tmp$log2FoldChange[i] < -fc_cutoff){
      res_tmp$cols[i] <- "gray47" 
    }
    else if(res_tmp$pvalue[i]>alpha & res_tmp$log2FoldChange[i] < fc_cutoff){
      res_tmp$cols[i] <- "gray10" 
    }
}
  

res_tmp$ENSG <- rownames(res_tmp)

# generate the splot 
p = ggplot(res_tmp, aes(log2FoldChange, -log10(pvalue))) + 
    geom_point(aes(col=col), alpha = 0.5, size =2.5, colour = res_tmp$cols, fill = res_tmp$cols)  + 
    xlab("Log2 fold change") + ylab("-log10 Q-value") +
    ylim(0, 9) + 
    xlim(-5, 11) +
    geom_hline(yintercept = -log10(alpha), color = "black", linetype = "dashed", size = 0.4) + 
    theme(legend.key = element_blank()) + 
    ggtitle("Control vs Dex") 

# print the plot 
print(p)
```

```
## Warning: Removed 1 rows containing missing values (geom_point).
```

<img src="day2-PART2_files/figure-html/unnamed-chunk-21-1.png" style="display: block; margin: auto;" />

This is nice, but some labels for potentially interesting genes would be useful. Lets add some using the **ggrepel** package. 

```r
p2 <- p + 
  # add labels to genes w/ LFC > 2 and above alpha threshold
  geom_label_repel(data = subset(res_tmp, log2FoldChange > 2 & pvalue < alpha), aes(label = gene), 
                     box.padding   = 0.35,
                     nudge_x = 0.1,
                     nudge_y = 0.1,
                     point.padding = 1,
                     label.size = 0.1,
                     segment.size = 0.3,
                     segment.color = 'grey50', size = 3) +
  # add labels to genes w/ LFC < -2 and above alpha threshold
  geom_label_repel(data = subset(res_tmp, log2FoldChange < -2 & pvalue < alpha), aes(label = gene), 
                     box.padding   = 0.35,
                     nudge_x = -0.1,
                     nudge_y = 0.1,
                     point.padding = 1,
                     label.size = 0.1,
                     segment.size = 0.3,
                     segment.color = 'grey50', size = 3) +
  # add vertical fold change lines 
  geom_vline(xintercept = fc_cutoff, colour = "black", linetype="dotted") + 
  geom_vline(xintercept = -fc_cutoff, colour = "black", linetype="dotted")

# print the plot 
print(p2)
```

```
## Warning: Removed 1 rows containing missing values (geom_point).
```

```
## Warning: Removed 1 rows containing missing values (geom_label_repel).
```

<img src="day2-PART2_files/figure-html/unnamed-chunk-22-1.png" style="display: block; margin: auto;" />

This looks a lot better, and gives us a lot more information than the first, very basic plot we generated. 

Food for thought: detecting truly differentially expressed genes is dependent on the technical variance between your replicates. If the technical variance is high, you generally need a large fold-change to achieve statistical significance. The more replicates you have, the more you are able to reduce this technical variance, which increases your statistical power, and enables you to confidently detect differential expression of smaller fold changes. For example, for an experiment where there are 300 truly differentially expressed genes between your conditions, you may detect 200 of these with 3 replicates, while you may detect 250 with 5 replicates. 

Save our results to .csv files

```r
# subset @ 5% adjusted pval sig. level 
res_order_FDR_05 <- res_ord[res_ord$padj<0.05,]
nrow(res_order_FDR_05)
```

```
## [1] 1673
```

```r
# write both to csv files
write.csv(as.data.frame(res_ord), file= "DE_results.csv")
write.csv(as.data.frame(res_order_FDR_05), file="DE_results.FDR.0.05.csv")
```

#### Why must we correct for multiple hypothesis testing? 

P-values are defined as the probability that we would observe a result as extreme as the one we observed, simply due to chance. In the case of RNA-seq, we are testing the probability that we would observe the log2 FC that we do for a given gene, if this result is due to chance. Therefore if we use 0.05 as a P-value threshold, and we test 20,000 genes for DE, this means that 5% of those genes we tested will have a log 2 FC that has a P-value < 0.05 simply due to chance. 5% of 20,000 is 1000 genes, which is obviously an unacceptable amount of false-positives. 

We address this problem through multiple testing correction. While several methods that control different aspects of the multiple testing problem, we commonly use methods that control the false-discovery rate (FDR) in RNA-seq DE experiments. Controlling the false discovery rate at 10% means that we are accepting that 1 in 10 of the genes with a significant adjusted P-value, is actually a false-positive, and not truly differentially expressed. RNA-seq DE studies are usually hypothesis generating in nature, so this is usually an acceptable compromise, however if your experiment requires more stringency, you may wish to use a method that controls the family-wise error rate (FWER), such as **Bonferonni** correction. 

Lets work through an example to demonstrate the importance of multiple tetsing. We will create a dataset with scrambeled sample labels, so that the null hypothesis (there is no differential expression) is true for all the genes. How many genes with an unadjusted P-value < 0.05 do you think we will get? 


```r
# create a new object and scramble the sample labels 
dds2 <- dds

# take random sample without replacement to get a scrambeled design variable 
colData(dds2)$group <- sample(colData(dds2)$group, length(colData(dds2)$group))

# check sample number in each group is the same 
table(colData(dds)$group)
```

```
## 
## untreated       Dex       Alb   Alb_Dex 
##         4         4         4         4
```

```r
table(colData(dds2)$group)
```

```
## 
## untreated       Dex       Alb   Alb_Dex 
##         4         4         4         4
```

```r
# re-run the DEseq2 analysis using the new group variable as the design variable 
dds2 <- DESeq(dds2)

# extract the DEG results just like before 
res2 <- results(dds2, 
  name = "group_Dex_vs_untreated", 
  alpha = 0.05, 
  lfcThreshold = 0)

# drop the NA values in P-value column 
res2 <- as.data.frame(res2)
res2 <- res2[-which(is.na(res2$padj)),]
```


```r
# how many P-values < 0.05 
sum(res2$pvalue < 0.05, na.rm=TRUE)
```

```
## [1] 325
```

```r
# how many FDR adjusted P-values < 0.05 
sum(res2$padj < 0.05, na.rm=TRUE)
```

```
## [1] 0
```

```r
# how many with Bonferonni adjusted P-values < 0.05 
sum(res2$pvalue < (0.05/nrow(res2)), na.rm=TRUE)
```

```
## [1] 0
```

```r
# plot the results 
plot(res2$log2FoldChange, -log10(res2$pvalue), 
     main = "Volcano plot - DEG w/ scrambled sample labes", 
     las = 1, col = "cornflowerblue",
     ylab = "- log10 P-value", xlab = "log2 Fold change", ylim = c(0,7))

# add significance lines 
abline(h= -log10(0.05), lty = 2, col = "red") # nominal P-value 
abline(h= -log10(0.05/nrow(res2)), lty = 2, col = "black") # Bonferonni 
```

<img src="day2-PART2_files/figure-html/unnamed-chunk-25-1.png" style="display: block; margin: auto;" />

You can see that there are minimal results with statistical signficance after correction, which is true since we scrambled the sample labels and created a fake dataset that should have no true DE. However, if we used the unadjusted P-values, we would identify A LOT of potentially interesting genes, that would infact be false-positives. 

This example highlights the short coming of hypothesis testing approaches, and demonstrates how important it is to correct for multiple hypothesis testing. 

### Other visualizations - MA plots

MA plots are also useful ways to visualize results from a DE analysis of RNA-seq data. These involve plotting the log2 fold-change (the so called M-value, representing the *M* in *MA-plot*) against the average expression level of a gene (the *A* in *MA-plot*). The MA-plot allows us to inspect the full range of expression values over which we detected significant DEGs, and what the magnitude of these fold-changes is. In a typical experiment, we expect to see DEGs across most of the range of expression values. To help identify genes that were significantly DE, any gene with an adjusted P-value of < 0.05 (or whatever threshold is set) is colored in red. 

```r
plotMA(res_ord, ylim=c(-6,6), main = "Raw Log2 Fold change")
```

![](day2-PART2_files/figure-html/unnamed-chunk-26-1.png)<!-- -->

The log2 fold-change plotted above is the raw LFC vaue estimated by the negative binomial GLM that we used in modeling. However, as we discussed above, the individual estimates of variance or dispersion for a single gene are often unreliable, and this holds true `log2 fold change` also. 

To obtain more useful LFC estimates, DESeq2 performs a statsitical procedure that involves shrinking the raw fold change estimates toward zero for genes that are less likely to contain reliable or highly important information. This is done in a very similar way to the shrinkage using empirical bayes that we discussed for the dispersion estimates. For shrinking LFC values, LFCs are penalized for properties such as:  
- low count values 
- high dispersion (& thus reduced confidence in expression levels)

DESeq2 provides a function `lfcShrink()` that must be implemented separately of the standard workflow implemented using `DESeq2()`. 

```r
# calculate shrunken fold change estimate
res_shrink <- lfcShrink(dds, 
                    coef=paste0(resultsNames(dds)[which(resultsNames(dds)=="group_Dex_vs_untreated")]), 
                    type="apeglm")
```

```
## using 'apeglm' for LFC shrinkage. If used in published research, please cite:
##     Zhu, A., Ibrahim, J.G., Love, M.I. (2018) Heavy-tailed prior distributions for
##     sequence count data: removing the noise and preserving large differences.
##     Bioinformatics. https://doi.org/10.1093/bioinformatics/bty895
```

After performing the shrinkage procedure, we compare the raw and shrunken LFCs to assess the impact of shrinkage. 

**Raw estimates of log2 FC:**

```r
plotMA(res_ord, ylim=c(-6,6), main = "Raw Log2 Fold change")
```

<img src="day2-PART2_files/figure-html/unnamed-chunk-28-1.png" style="display: block; margin: auto;" />

**Shrunken estimates of log2 FC:**

```r
plotMA(res_shrink, ylim=c(-6,6), main = "Shrunken Log2 Fold change")
```

<img src="day2-PART2_files/figure-html/unnamed-chunk-29-1.png" style="display: block; margin: auto;" />


We can see that significantly DE genes are detected across the full range of expression values (x-axis), which is a good sign that our differential expression modeling has worked well. We can also see that we have a handful of genes with larger expression values (> LFC 2) which potentially represent the most important individual genes, while the majority of our DEGs have a LFC < 1.5 (ish). 

Comparing to the raw LFCs, we can also see that the majority of genes with lower expression values have have their LFCs shrunk toward zero. This is important as genes with low counts may simply end up with a large LFC since this is easy to do at small count values, but these are unlikely to be accurate fold-changes, so we don't want to prioritize their importance by giving them a large LFC. It's always good to look at the shrunken estimates, to confirm that you don't have a lot of DEGs at very small count values. If you do, you may want to look at the expression levels for those genes to investigate these findings in more detail. 

As the mean or counts increase, it is evident that the level of shrinkage is less, although may still be high for genes with greater dispersion estimates. As we move toward the more highly expessed genes, you can see how more genes at lower fold change values are able to be identified as significant, which is due to the fact that there is more information avaiable for these genes, so we can be more confident during hypothesis tetsing of these genes. 

**Note:** This shrinkage does not really change the hypothesis testing, therefore is performed independently, as is for use in prioritizing your results further for visual inspection or some sort of functional analysis (e.g. pathway analysis). 

#### Hierachical clustering on the DEGs

A final visualization that is useful to generate is a heatmap based on unsupervised hierachical clustering of the DEGs identified. We can do this by limiting the matrix of rlog values to only those for the DEGs, and then performing the clustering specifically on these data. 

```r
rld <- rlog(dds, blind = FALSE)
ind_to_keep <- c(which(colData(rld)$group=="untreated"), which(colData(rld)$group=="Dex"))

# set up gene expression matrix 
mat1 <- assay(rld)[rownames(res_order_FDR_05), ind_to_keep]

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
draw(ht1, row_title = "Genes", column_title = "Hierachical clustering of DEGs (padj<0.05)")
```

<img src="day2-PART2_files/figure-html/unnamed-chunk-30-1.png" style="display: block; margin: auto;" />

***

## Recap of the full workflow

In the above analysis we got into the details of how the statistics behind differential expreission work and ran some extra code to demonstrate the utility of these statistics. However, the entire DESeq2 workflows boils down to just a few functions run sequentially. Lets do a quick recap of these functions to help consolidate what we have learnt. 

Read in the data:

```r
cts <- as.matrix(read.table("Day-2/all_counts.txt", 
                            sep="\t", header = TRUE, row.names=1, 
                            stringsAsFactors = F))
```

Read in the metadata:

```r
sra_res <- read.csv("Day-2/sra_result.csv", row.names=1)
sra_res$Sample <- sra_res$Sample.Accession
sra_run <- read.csv("Day-2/SraRunInfo.csv", row.names=1)
```

Construct a DESeq2 dataset from the raw counts, the metadata, and the desired design variable to be tested for differential expression. 

```r
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = colData,
                              design = ~ group)
```

Apply the DESeq2 analysis pipeline:

```r
dds <- DESeq(dds)
```

Perform regularized log transformation:

```r
rld <- rlog(dds, blind = FALSE)
```

Use rlog to perform exploratory analyses: 
- Principal components analysis (PCA)
- Unsupervised hierachical clustering

Check the disperion estimates to evalute model fit:

```r
plotDispEsts(dds)
```

Extract, order, annotate, and subset the results from the DESeq2 object

```r
res <- results(dds, 
  name = "group_Dex_vs_untreated", 
  alpha = 0.05, 
  lfcThreshold = 0)

# order by adj Pval 
res_ord <- res[order(res$padj),] 

# add gene annotation to results 
anno <- read.delim("Day-2/GRCh38.p12_ensembl-97.txt", stringsAsFactors = T, header = T)
anno <- anno[order(anno$Chromosome.scaffold.name),]
mat1 <- match(rownames(res_ord), anno$Gene.stable.ID)
res_ord$gene <- as.character(anno$Gene.name[mat1])

# subset results for only genes with adjusted P-values < 0.05
res_order_FDR_05 <- res_ord[res_ord$padj<0.05,]
```

Perform empirical bayes shrinkage of raw fold-change estimates: 

```r
res_shrink <- lfcShrink(dds, 
                    coef=paste0(resultsNames(dds)[which(resultsNames(dds)=="group_Dex_vs_untreated")]), 
                    type="apeglm")
```

Generate visualizations: 
- MA plots (raw vs shrunken fold-changes)
- Volcano plots
- Heatmaps 

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
##  [3] apeglm_1.8.0                EnhancedVolcano_1.4.0      
##  [5] ggrepel_0.8.2               ggplot2_3.3.0              
##  [7] circlize_0.4.8              readr_1.3.1                
##  [9] ComplexHeatmap_2.2.0        RColorBrewer_1.1-2         
## [11] gplots_3.0.3                pheatmap_1.0.12            
## [13] dplyr_0.8.5                 vsn_3.54.0                 
## [15] biomaRt_2.42.0              DESeq2_1.26.0              
## [17] SummarizedExperiment_1.16.1 DelayedArray_0.12.2        
## [19] BiocParallel_1.20.1         matrixStats_0.55.0         
## [21] Biobase_2.46.0              GenomicRanges_1.38.0       
## [23] GenomeInfoDb_1.22.0         IRanges_2.20.2             
## [25] S4Vectors_0.24.3            BiocGenerics_0.32.0        
## [27] tximport_1.14.0            
## 
## loaded via a namespace (and not attached):
##   [1] colorspace_1.4-1       rjson_0.2.20           htmlTable_1.13.3      
##   [4] XVector_0.26.0         GlobalOptions_0.1.1    base64enc_0.1-3       
##   [7] clue_0.3-57            rstudioapi_0.11        farver_2.0.3          
##  [10] affyio_1.56.0          bit64_0.9-7            mvtnorm_1.1-0         
##  [13] AnnotationDbi_1.48.0   xml2_1.2.2             splines_3.6.1         
##  [16] geneplotter_1.64.0     knitr_1.28             Formula_1.2-3         
##  [19] annotate_1.64.0        cluster_2.1.0          dbplyr_1.4.2          
##  [22] png_0.1-7              BiocManager_1.30.10    compiler_3.6.1        
##  [25] httr_1.4.1             backports_1.1.5        assertthat_0.2.1      
##  [28] Matrix_1.2-18          limma_3.42.2           acepack_1.4.1         
##  [31] htmltools_0.4.0        prettyunits_1.1.1      tools_3.6.1           
##  [34] coda_0.19-3            gtable_0.3.0           glue_1.3.1            
##  [37] GenomeInfoDbData_1.2.2 affy_1.64.0            rappdirs_0.3.1        
##  [40] Rcpp_1.0.3             bbmle_1.0.23.1         vctrs_0.2.3           
##  [43] gdata_2.18.0           preprocessCore_1.48.0  xfun_0.12             
##  [46] stringr_1.4.0          rvest_0.3.5            lifecycle_0.2.0       
##  [49] gtools_3.8.1           XML_3.99-0.3           MASS_7.3-51.5         
##  [52] zlibbioc_1.32.0        scales_1.1.0           hms_0.5.3             
##  [55] yaml_2.2.1             curl_4.3               memoise_1.1.0         
##  [58] gridExtra_2.3          emdbook_1.3.12         bdsmatrix_1.3-4       
##  [61] rpart_4.1-15           latticeExtra_0.6-29    stringi_1.4.6         
##  [64] RSQLite_2.2.0          genefilter_1.68.0      checkmate_2.0.0       
##  [67] caTools_1.18.0         shape_1.4.4            rlang_0.4.5           
##  [70] pkgconfig_2.0.3        bitops_1.0-6           evaluate_0.14         
##  [73] lattice_0.20-40        purrr_0.3.3            labeling_0.3          
##  [76] htmlwidgets_1.5.1      bit_1.1-15.2           tidyselect_1.0.0      
##  [79] plyr_1.8.6             magrittr_1.5           R6_2.4.1              
##  [82] Hmisc_4.3-1            DBI_1.1.0              withr_2.1.2           
##  [85] pillar_1.4.3           foreign_0.8-76         survival_3.1-11       
##  [88] RCurl_1.98-1.1         nnet_7.3-13            tibble_2.1.3          
##  [91] crayon_1.3.4           KernSmooth_2.23-16     BiocFileCache_1.10.2  
##  [94] rmarkdown_2.1          jpeg_0.1-8.1           GetoptLong_0.1.8      
##  [97] progress_1.2.2         locfit_1.5-9.1         data.table_1.12.8     
## [100] blob_1.2.1             webshot_0.5.2          digest_0.6.25         
## [103] numDeriv_2016.8-1.1    openssl_1.4.1          munsell_0.5.0         
## [106] viridisLite_0.3.0      askpass_1.1
```

