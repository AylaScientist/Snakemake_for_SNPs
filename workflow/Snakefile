#include: "rules/common.smk"

configfile: '../config/config.yaml'

print(config)

##### Target rules #####

import pandas as pd
import os

sample_names = pd.read_csv ( config['Sample_names'] )
pseudo_codes = pd.read_csv ( config['Pseudogenome_codes'])
# Create arrays of the sample names
samples = sample_names['Sample_name'].values
pseudos = pseudo_codes['PSGs'].values
print(samples)


rule all:
    input:
          o1= "results/SNPs_function.csv",
          o2= "results/SNP_dictionary.csv",
          o3= "results/Treatment_SNPs_sig_all_for_Venn.csv",
          o4= "results/Summary_of_polymorphisms.csv",
          o5= "results/Intronic_SNPs_caryotype.csv",
          bam=expand("recal/{sample}_ref.bam", sample = samples),
          bqsr=expand("recal_tables/{sample}_after_ref.table", sample = samples),
          valid_bqsr=expand("logs/picard/ValidateSamFile/Validate_bqsr_{sample}.log", sample = samples),
          valid_bqsr_PSG=expand("logs/picard/ValidateSamFile_{pseudo}/Validate_bqsr_{sample}_{pseudo}.log", sample = samples, pseudo = pseudos),
          recal_table=expand("recal_tables/{sample}_after_ref.table", sample = samples),
          degs=expand("gene_counts/{sample}.count_ref.txt", sample = samples),
          PSG_dict = expand("pseudogenomes/{pseudo}/{pseudo}_GCF_001858045.2.dict", pseudo = pseudos),
          PSG_index= expand("pseudogenomes/{pseudo}/SAindex", pseudo = pseudos),
          PSG_fai= expand("pseudogenomes/{pseudo}/{pseudo}_GCF_001858045.2.fa.fai", pseudo = pseudos),
          #error_mae ="results/SNPs_ready_mae_valid_genome.csv",
          error="results/SNPs_ready_no_mae.csv",


##### Modules #####
###################

# Preparation #
###############
include: "rules/trimmomatic.smk"#In the trimming R3 input is called 2 for the output
include: "rules/star_gi_ref.smk"

# First SNP callling on the reference genome #
##############################################
include: "rules/star_pe_ref.smk"        # map with STAR
include: "rules/replace_rg_ref.smk"     # Replace or add read Groups:
include: "rules/mark_duplicates_ref.smk"# Mark and eliminate Duplicates:
include: "rules/splitNcigar_ref.smk"    # Split Ncigars
inlcude: "rules/validate_bam_ref.smk"   # Validate bam
include: "rules/htseq.smk"              # Make count matrix from each sample for differentially expressed genes
include: "rules/haplotype_caller_ref.smk"# Haplotype caller
include: "rules/combine_gvcfs.smk"      # Combine gvcf
include: "rules/genotype_gvcfs.smk"     # Genotype gvcf
include: "rules/gatk_select.smk"        # Select biallelic sites
include: "rules/gatk_variantstotable.smk" # Make the table


## Pseudogenomes pipeline ##
############################
# From here it starts the part of the pipeline that creates the pseudogenomes
# and prepares the vcf files for the elimination of the mapping bias.

# Subset by keeping positions with GQ >= 30 and DP >=5
# Not keeping indels
include: "rules/subset_ref.smk"         # Subset the genome
include: "rules/pseudogenomes.smk"      # Create the Pseudogenomes

## Now you can start the pipeline again from the fastq file for the pseudogenomes:
include: "rules/star_pseudogenomes.smk"       # Map against pseudogenomes
include: "rules/replace_rg_pseudo.smk"        # Replace or add read groups
include: "rules/mark_duplicates_pseudo.smk"    # Mark and eliminate duplicates
include: "rules/splitNcigar_pseudo.smk"       # Split Ncigars
include: "rules/validate_bam_pseudo.smk"      # Validate bam
include: "rules/haplotype_caller_pseudo.smk"  # Haplotype caller
include: "rules/combine_gvcfs_pseudo.smk"     # Combine gvcf
include: "rules/genotype_gvcfs_pseudo.smk"    # Genotype gvcf
include: "rules/select_pseudo.smk"            # Select biallelic variants

## Annotate vcf files ##
########################
# Convert th vcf file to annovar format thus extracting the SNPs
# Make tokens for the annotation script
# Collect the annotations in the db
# Make tokens for the tables with the annotation
# Annotate the file. This file can be used for the creation of the pseudogenomes
include: "rules/annotation.smk"


## Data science ##
##################
## Create the count table for the data wrangling
include: "rules/variants_to_table.smk"  # Make the table
include: "rules/table2df.smk"           # Preapre the tables for the workflow format

# Workflow
include: "rules/workflow.smk"           # Processing the data
# WARNING! Need to chose one of the next:
#include: "rules/workflow_mae.smk"       # WARNING! Considering monoallelic expression in different tissues from teh same individual
include: "rules/workflow_one_tissue.smk"# WARNING! Only one tissue

# Stats only in non-mae
include: "rules/stats.smk"              # Apply Chi^2, Fisher and Binomial on NO MAE

# Tables and figures for publication
include: "rules/results.smk"
include: "rules/results_mae.smk"        # For MAE
