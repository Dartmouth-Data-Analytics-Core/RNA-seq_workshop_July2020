# Welcome to the DAC RNAseq workshop #

Before you attend the workshop there are a couple of things we would like you to do to get setup for using the tools that will be required during the course of the workshop.  

For those of you that indicated that you did not have an account on *discovery* you should have received an email from me explaing how to set that up, please make sure this is done and you are able to log into your account **BEFORE** the workshop begins. 

## Downloading the data ##

For this workshop we will be using a dataset downloaded from the short read archive (SRA), a public repository of genomic data. This dataset comes from [this paper](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0099625), and was collected from human airway smooth muscle cells to test gene pathways effected by exposure to Glucocorticoid drugs, which have been historically used for their anti-inflammatory effects to treat asthma. Four cell lines were treated with either a control vehicle (untreated), dexamethasone (dex), albuterol (alb), or both dexamethasone and albuterol (co-treated) for 18 hours before transcriptomes were extracted.

The commands that you will be following can be found in markdown `(.md)` files where there is a brief description of the command and how it is applied to the data and what it does followed by an example command that you can copy and paste into the terminal window. The majority of day 1 will be using the terminal window on your local machine, with an open `ssh` connection to discovery7, as we will be running `bash` code. For day 2, you will be using RStudio on your local machine to run the R-markdown files (`.rmd`) located in this GitHub repo. 

In your terminal window navigate to where you want to download the files needed for this workshop onto your local machine. Then execute the following command:

```bash
git clone https://github.com/Dartmouth-Data-Analytics-Core/RNA-seq_workshop_July2020/
```

## Setting up a Conda Environment ## 

Conda is a package management system that helps you find, install, and organize groups of packages needed for a particular task. Conda environments are really useful when working on high performance computing (HPC) environments like Dartmouth's Discovery system because you can install packages locally without needing administrator permission. Conda environments are also useful for project continuity, the versions of the packages that you install in your environment and all of their dependencies will remain the same (unless you update them). We will be using a conda environment to make sure we all have the same version of many different bioinformatics software programs available to us. 

Before you begin using conda environments on discovery you will need to ensure that you have run the source code to enable the conda commands. Log into discovery and run the following command:

```bash
source /optnfs/common/miniconda3/etc/profile.d/conda.sh
```

We recommend that you add the above line of code to your `.bashrc` file in your home directory, otherwise you will need to run this command each time you start a new session on discovery. You can follow the [tutorial](https://services.dartmouth.edu/TDClient/1806/Portal/KB/ArticleDet?ID=72888) on the *Research Computing* teams website. 

Next you will have to run the following command to create a .conda/ directory in your home drive to store all of your personal conda environments. You only have to run this command once to make this directory, so it does not need to be added to your .bashrc file.

```bash
cd ~
mkdir -p .conda/pkgs/cache .conda/envs
```

Now you will need to create the conda environment that we will be using for the course into your personal directory of accessible conda environments.  This takes about 15 minutes to execute and you will see all of the packages that are loaded into this environment. The number of packages should indicate why conda environments are so useful, imagine having to load all of these packages individually it is much easier to load them with a single command in a conda environment.

```bash
conda env create -f /dartfs-hpc/scratch/rnaseq1/environment.yml
```

When you are ready activate the conda environment, which you will need for the work we are doing for day 1 of the workshop you can use the following command. 

```bash
conda activate rnaseq1/
```

There is one more conda environment you will need to create in order to run some of the QC metrics for our alignment. We will install the picard program (yes from Star Trek) using the bioconda channel. 

```bash
conda install -c bioconda picard
```
During this installation you will be shown a list of all the dependencies that picard requires (this is why conda environments are so handy) and asked if you would like to proceed with installing these packages in your picard conda environment, type yes to install these packages in our conda environment. 

You will see that the activate command has worked when it reads (rnaseq1) rather than (base) to the left of the prompt. When you are finished using a conda environment it is good practice to deactivate your session with the following command.

```bash
conda deactivate
```


## Setting up an R project ##

We will be using R-Studio for differential expression analysis on day 2, we ask that you have R and R-Studio installed prior to attending the workshop, you will need to be running at least version 3.6 to nesure all of the packages needed will run smoothly. The latest versions for R and R-Studio can be found [here](https://cran.r-project.org) and [here](https://rstudio.com/products/rstudio/download/).

Next you will need to set up a new project in R-Studio for this workshop. Projects in R are like containers for various jobs that you will perform, the history of the project will be loaded when you open a new project. By using a project you can install all of the tools ahead of time and they will be there for you when you need them during the workshop. In R-Studio under File select New directory and then select New Project and name the project something you will remember (RNAseq_workshop).

Now that you have your project loaded run the following code to ensure the required packages (tximport, DESeq2, biomaRt, vsn, complexHeatmap, kableExtra, rlang, enhancedVolcano, apeglm) are installed:

```r
if (!any(rownames(installed.packages()) == "tximport")){
  if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
  BiocManager::install("tximport")
}
library(tximport)

if (!any(rownames(installed.packages()) == "DESeq2")){
  if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
  BiocManager::install("DESeq2")
}
library(DESeq2)

if (!any(rownames(installed.packages()) == "biomaRt")){
  if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
  BiocManager::install("biomaRt")
}
library(biomaRt)

if (!any(rownames(installed.packages()) == "vsn")){
  if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
  BiocManager::install("vsn")
}
library(vsn)

if (!any(rownames(installed.packages()) == "ComplexHeatmap")){
  if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
  BiocManager::install("ComplexHeatmap")
}
library(ComplexHeatmap)

if (!any(rownames(installed.packages()) == "kableExtra")){
      install.packages("kableExtra")
}
library(kableExtra)

if (!any(rownames(installed.packages()) == "rlang")){
  if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
  BiocManager::install("rlang")
}
library(rlang)

if (!any(rownames(installed.packages()) == "enhancedVolcano")){
  if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
  BiocManager::install("enhancedVolcano")
}
library(enhancedVolcano)

if (!any(rownames(installed.packages()) == "apeglm")){
  if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
  BiocManager::install("apeglm")
}
library(apeglm)

sessionInfo()
```
