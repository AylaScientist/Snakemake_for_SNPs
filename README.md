# Snakemake for SNPs: A pipeline for calling SNPs and quantify them in an unbiased manner
Snakemake for SNPs is a flexible and user-friendly SNPs analysis workflow.

[![Snakemake](https://img.shields.io/badge/snakemake-≥7.18.2.1-brightgreen.svg)](https://snakemake.github.io)
[![GitHub actions status](https://github.com/<owner>/<repo>/workflows/Tests/badge.svg?branch=main)](https://github.com/<owner>/<repo>/actions?query=branch%3Amain+workflow%3ATests)

Snakemake for SNPs can be applied to both model and non-model organisms. It supports mapping RNA-Seq raw reads to the reference genome (can be downloaded from public database or can be homemade by users) and it can do both Allele Specific Expression for SNPs and obtain Differential Expressed Genes (DEGs), which in turn can be cross between them. It requires basic python programming skill for use. If you're beginner at programming, just jump on the config file and adapt it to your experiments!

If you use our pipeline you need to cite us:

WARNING: adapt the citation to our link:

## Workflow
<img src="https://github.com/AylaScientist/Snakemake_for_SNPs/blob/master/Figure%201%20Pipeline%20white%20background.png" width="800">

The usage of this workflow is described in the [Snakemake Workflow Catalog](https://snakemake.github.io/snakemake-workflow-catalog/?usage=<owner>%2F<repo>).


## Quick start

Clone the repository:

#git clone https://github.com/AylaScientist/Snakemake_for_SNPs.git

Create the environment:

`conda create -n pipeline python=3.7`

Activate the environment:

`conda activate pipeline`

### Installation
Install the packages including the bio tools:

`pip install git+https://github.com/snakemake/snakemake`

`conda install -c bioconda trimmomatic=0.39`

`conda install -c bioconda fastqc=0.11.9`

`conda install -c bioconda star=2.7.10a`

`conda install -c bioconda htseq=0.11.3`

`conda install -c bioconda picard=2.26`

`conda install -c bioconda gatk4=4.2.5.0`

`conda install -c bioconda samtools=1.6`

`conda install -c bioconda bcftools=1.9`

`conda install -c bioconda vcftools=0.1.16`

`conda install -c anaconda perl=5.26.2`

`conda install -c anaconda pandas`

`conda install -c conda-forge matplotlib`

`conda install -c conda-forge py-bgzip`


### Run a dry run for the pipeline with the example data to build a dag of jobs
If this point doesn't work, please contact me: ayla.bcn@gmail.com
`snakemake -n `

### Run the pipeline with the desired resources. This is an example for 4 threads at 4GB
`snakemake --cores 4 --mem_mb 40000 `


### Set up configuration for your personal project
Customize the workflow based on your need in `configs/config_main.yaml`.
In this file you should also change the species and the different databases for gene/transcript/protein/GO_function/KEGG correct annotation and mining of the data

Modify the metafiles describing your data and the experiment:

`config/Experiemntal_design.csv`

`config/Experiemntal_groups.csv`

`config/Pseudogenome_codes.csv`

`config/Sample_names.csv`

`config/Samples_MAE.csv`

`config/samples.csv`

`config/tb1_colnames.csv`

`config/tb2_colnames.csv`


## Evaluation
The pipeline for SNPs has been evaluated on 4 datasets including 2 non-model organism (Nile and Mozambique tilapias). 
WARNING: Put here the link to the article
