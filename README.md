# Snakemake for SNPs: A pipeline for calling SNPs and quantify them in an unbiased manner
RASflow is a modular, flexible and user-friendly RNA-Seq analysis workflow.

[![Snakemake](https://img.shields.io/badge/snakemake-≥6.3.0-brightgreen.svg)](https://snakemake.github.io)
[![GitHub actions status](https://github.com/<owner>/<repo>/workflows/Tests/badge.svg?branch=main)](https://github.com/<owner>/<repo>/actions?query=branch%3Amain+workflow%3ATests)

Snakemake for SNPs can be applied to both model and non-model organisms. It supports mapping RNA-Seq raw reads to the reference genome (can be downloaded from public database or can be homemade by users) and it can do both Allele Specific Expression for SNPs and obtain Differential Expressed Genes (DEGs), which in turn can be cross between them. It requires basic python programming skill for use. If you're beginner at programming, just jump on the config file and adapt it to your experiments!

If you use our pipeline you need to cite us:

WARNING: adapt the citation to our link:

## Workflow
#<img src="https://github.com/zhxiaokang/RNA-Seq-analysis/blob/master/workflow/workflow_chart.jpg" width="450">

The usage of this workflow is described in the [Snakemake Workflow Catalog](https://snakemake.github.io/snakemake-workflow-catalog/?usage=<owner>%2F<repo>).


WARNING: Prepare the repository, make the step by step with conda and let aside the Doker

## Quick start
### Installation
#### Manual mode

Clone the repository:

#git clone https://github.com/zhxiaokang/RASflow.git`

Create the environment:

`conda env create -n pipeline -f env.yaml`

Activate the environment:

`conda activate pipeline`


### Set up configuration
Modify the metafiles describing your data and the experiment:
`config/Experiemntal_design.csv`
`config/Experiemntal_groups.csv`
`config/Pseudogenome_codes.csv`
`config/Sample_names.csv`
`config/Samples_MAE.csv`
`config/samples.csv`
`config/tb1_colnames.csv`
`config/tb2_colnames.csv`

Customize the workflow based on your need in `configs/config_main.yaml`.

WARNING: How to run the pipeline

### Run a dry run for the pipeline
`snakemake -n `

### Run the pipeline with the desired resources. This is an example for 4 threads at 4GB
snakemake --cores 4 --mem_mb 40000

## Tutorial
A more detailed tutorial of how to use this workflow can be found here: [Tutorial](https://github.com/zhxiaokang/RASflow/blob/master/Tutorial.pdf)

## Evaluation
The pipeline for SNPs has been evaluated on 4 datasets including one model organism (human) and 2 non-model organism (Nile and Mozambique tilapias). To keep this repository as light as possible, the evaluation of the pipeline on real datasets is deposited here:
WARNING: Put here the
