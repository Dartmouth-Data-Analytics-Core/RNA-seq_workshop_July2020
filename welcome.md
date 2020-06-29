# Welcome to the DAC RNAseq workshop #

Before you attend the workshop there are a couple of things we would like you to do to get setup for using the tools that will be required during the course of the workshop.
For those of you that indicated that you did not have an account on discovery you should have received an email from me explaing how to set that up, please make sure this is done and you are able to log into your account BEFORE the workshop begins. 

## Setting up a Conda Environment ## 

Conda is a package management system that helps you find, install, and organize groups of packages needed for a particular task. Conda environments are really useful when working on high performance computing (HPC) environments like Dartmouth's Discovery system because you can install packages locally without needing administrator permission.Conda environments are also useful for project continuity, the versions of the packages that you install in your environment and all of their dependencies will remain the same (unless you update them).

Before you begin using conda environments on discovery you will need to ensure that you have run the source code to enable the conda commands. 

```bash
source /optnfs/common/miniconda3/etc/profile.d/conda.sh
```

We recommend that you add this line of code to your .bashrc file in your home directory, otherwise you will need to run this command each time you start a new session on discovery.

Next you will have to run this command to create a .conda/ directory in your home drive to store all of your personal conda environments. You only have to run this command once to make this directory, so it does not need to be added to your .bashrc file.

```bash
cd ~
mkdir -p .conda/pkgs/cache .conda/envs
```

Now you will need to load the conda environment that we will be using for the course into your personal directory of accessible conda environments. 

```bash
conda env create -f environment.yml
```

When you are ready to activate the conda environment, which you will need for the work we are doing for day 1 of the workshop you can use the following command. 

```bash
conad activate ranseq1/
```
You will see that the activate command has worked when it reads (rnaseq1) rather than (base) to the left of the prompt. 
When you are finished using a conda environment it is good practice to deactivate your session with the following command.

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
  if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
  BiocManager::install("kableExtra")
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

Next you will need to download the github repository with all of the files and information for day 2 by running the following command:
