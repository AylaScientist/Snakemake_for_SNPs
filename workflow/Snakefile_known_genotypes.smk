configfile: '../config/config_known_genotypes.yaml'

print("Configfile: ")
print(config)

##### Target rules #####

import pandas as pd
import os

sample_names = pd.read_csv ( config['Sample_names'] )
pseudo_codes = pd.read_csv ( config['Pseudogenome_codes'])
# Create arrays of the sample names
samples = sample_names['Sample_name'].values
pseudos = pseudo_codes['PSGs'].values
print("Sample names: ")
print(samples)
print("pseudogenomes:")
print(pseudos[0])
print(pseudos[1])
test = expand(config["input_files"]["r1"],sample = samples)
print(test)
test2 = expand("trimmed/{sample}_R1.fastq",sample = samples)
print(test2)


rule all:
    input:
        ## Zero call
        #bam=expand("calls/{sample}_ref.g.vcf", sample = samples),
        ## First call 
        #degs = expand("gene_counts/{sample}.count_ref.txt", sample = samples), # Counts matrix for the DEGs
        #psgd = expand(config['files_pseudogenomes']['PSG_dict'], pseudo = pseudos), # Dictionary of the pseudogenomes
        #psgi = expand(config['files_pseudogenomes']['PSG_index'], pseudo = pseudos), # Index of the pseudogenomes
        #psgf = expand(config['files_pseudogenomes']['PSG_fai'], pseudo = pseudos), # fai index of the pseudogenomes
        ## Uncomment for the second call:
        o1 = config['results']['o1'],
        o2 = config['results']['o2'],
        o3 = config['results']['o3'],
        o4 = config['results']['o4'],
        o5 = config['results']['o5'],
        ## Uncomment for the second call if you have the genotypes of the DNA
        #om1 = config['results_mae']['om1'],
        #om2 = config['results_mae']['om2'],
        #om3 = config['results_mae']['om3'],
        #om4 = config['results_mae']['om4'],
        #om5 = config['results_mae']['om5'],
        ## Helpers for debug
        #bam=expand("{sample}_ref_Aligned.sortedByCoord.out.bam",sample = samples),
        #bam=expand("calls/{sample}_ref.g.vcf", sample = samples),
        #bam0 = expand ("{sample}_{pseudo}_Aligned.sortedByCoord.out.bam", sample = samples, pseudo = pseudos[0]),
        #bam = expand ("split/{sample}_{pseudo}.bam", sample = samples, pseudo = pseudos),
        #gvcf = expand ("calls/{sample}_{pseudo}.g.vcf", sample = samples, pseudo = pseudos),
        #gvcf="calls/all_ref_g.vcf",
        #vcf="calls/all_ref.vcf",
        #vcf="calls/selected_ref_RNA.vcf",
        #gvcfs=expand("calls/{sample}_{pseudo}.g.vcf", sample = samples, pseudo = pseudos),
        #gvcf0="calls/all_g_"+pseudos[0]+".vcf",
        #gvcf1="calls/all_g_"+pseudos[1]+".vcf",
        #vcf=expand("calls/all_{pseudo}.vcf", pseudo = pseudos),
        #vcf=expand("calls/selected_{pseudo}.vcf", pseudo = pseudos),
        #annotated=expand(config['params']['annotation']['convert']+'{pseudo}', pseudo = pseudos),
        #o1=expand(config['params']['annotation']['output_annotate']+"{pseudo}.variant_function", pseudo = pseudos),
        #o1=expand("annotated/annotated_all_snps_{pseudo}"+config['params']['annotation']['output'], pseudo = pseudos),
        #vcf=expand("variants/AD_GT_counts_bi_{pseudo}.table", pseudo = pseudos),
        #table1=expand("variants/AD_GT_counts_bi_{pseudo}_step1.csv", pseudo = pseudos),
        #csv=expand("workflow/AD_GT_counts_bi_{pseudo}.csv", pseudo = pseudos),
        #merged="workflow/merged_df.csv",
        #result1="results/Uniform_to_validate_SNPs_valid_genome.csv",
        #result2="results/SNP_panel_valid_genome.csv", #This SNP panel contains the SNPs found expressed in the experimental popilation
        #result3="results/Genotype_models_valid_genome.csv",
        #result_mae="results/SNPs_ready_mae.csv", # SNPs that follow monoallelic expression for one tissue
        #result_no_mae="results/SNPs_ready_no_mae.csv", # SNPs that do not follow monoallelic expression
        #vcf="variants/AD_GT_counts_bi.table",
        #vcf="pseudogenomes/subset_vcf_file_wo_indels.recode.vcf",
        #vcf="pseudogenomes/subset_vcf_file_wo_indels.recode.vcf.gz",
        #tabix="pseudogenomes/subset_vcf_file_wo_indels.recode.vcf.gz.tbi",
        #pseudo = expand("pseudogenomes/{pseudo}/{pseudo}_GCF_001858045.2.fa", pseudo = pseudos),
        #star = expand("{sample}_{pseudo}_Aligned.sortedByCoord.out.bam", sample = samples, pseudo = pseudos)
        #rg = expand("fixed-rg/{sample}_{pseudo}.bam", sample = samples, pseudo = pseudos),
        #rg_ref = expand("fixed-rg/{sample}.bam", sample = samples),
        #error = config['error'], # Estimate the error of your SNP call if you have genomic data with known SNPs
        #bam = expand(config['Other_intermediate_files']['bam'], sample = samples), 
        #bqsr = expand(config['Other_intermediate_files']['bqsr'], sample = samples),
        #val_bqsr = expand(config['Other_intermediate_files']['valid_bqsr'], sample = samples),
        #val_psg = expand(config['Other_intermediate_files']['valid_bqsr_PSG'], sample = samples, pseudo = pseudos),
        #recal = expand(config['Other_intermediate_files']['recal_table'], sample = samples)


##### MODULES #####
###################

# Preparation #
###############
ends = config['params']['star']['mode'] #Information about paired ends
print(ends)

# First SNP callling on the reference genome #
##############################################
"""
if ends == "SE":
    include: "rules/trimmomatic_se.smk"#In the trimming R3 input is called 2 for the output
    include: "rules/star_gi_ref.smk"
    # First SNP callling on the reference genome #
    ##############################################
    include: "rules/star_se_ref.smk"        # map with STAR
else:
    include: "rules/trimmomatic_pe.smk"#In the trimming R3 input is called 2 for the output
    include: "rules/star_gi_ref.smk"
    # First SNP callling on the reference genome #
    ##############################################
    include: "rules/star_pe_ref.smk"        # map with STAR



include: "rules/replace_rg_ref.smk"     # Replace or add read Groups
include: "rules/mark_duplicates_ref.smk"# Mark and eliminate Duplicates
include: "rules/splitNcigar_ref.smk"    # Split Ncigars
inlcude: "rules/validate_bam_ref.smk"   # Validate bam
include: "rules/htseq.smk"              # Make count matrix from each sample for differentially expressed genes
include: "rules/haplotype_caller_ref.smk"# Haplotype caller
include: "rules/combine_gvcfs.smk"      # Combine gvcf
include: "rules/genotype_gvcfs.smk"     # Genotype gvcf
include: "rules/gatk_select.smk"        # Select biallelic sites
include: "rules/gatk_variantstotable.smk" # Make the table


## Creation of Pseudogenomes ##
###############################
# From here it starts the part of the pipeline that creates the pseudogenomes
# and prepares the vcf files for the elimination of the mapping bias.

# Subset by keeping positions with GQ >= 30 and DP >=5
# Not keeping indels
include: "rules/subset_ref.smk"         # Subset the genome
include: "rules/pseudogenomes.smk"      # Create the Pseudogenomes

"""


## Second calling of SNPs on the pseudogenomes ##
#################################################
# Uncomment this section for the second call 
## Now you can start the pipeline again from the fastq file for the pseudogenomes:

if ends == "SE":
    include: "rules/star_se_pseudogenomes.smk"        # map with STAR
else:
    include: "rules/star_pe_pseudogenomes.smk"        # map with STARinclude: "rules/star_pseudogenomes.smk"       # Map against pseudogenomes

include: "rules/replace_rg_pseudo.smk"        # Replace or add read groups
include: "rules/mark_duplicates_pseudo.smk"   # Mark and eliminate duplicates
include: "rules/splitNcigar_pseudo.smk"       # Split Ncigars
include: "rules/validate_bam_pseudo.smk"      # Validate bam
include: "rules/haplotype_caller_pseudo.smk"  # Haplotype caller
include: "rules/combine_gvcfs_pseudo.smk"     # Combine gvcf
include: "rules/genotype_gvcfs_pseudo.smk"    # Genotype gvcf
include: "rules/select_pseudo.smk"            # Select biallelic variants


## Annotate vcf files ##
########################
# Convert the vcf file to annovar format thus extracting the SNPs
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
include: "rules/stats.smk"              # Apply Chi^2, Fisher and Binomial on NO MAE

# Tables and figures for publication
include: "rules/results.smk"
include: "rules/results_mae.smk"        # For study of the SNPs in MAE
