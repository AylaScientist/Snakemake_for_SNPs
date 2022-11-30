# Here you can include a step to separate sites that are homozygot for some sample
# and  monoallelic expression (MAE) in case you have the genotypes from the genome
# The monoallelic expression (MAE) can result from homozygots as well as imprinted genes.
# In order to distinguish each case, we need to know the genotype.
rule mae:
    input:
        result1="results/Uniform_to_validate_SNPs_valid_genome.csv",
        names="config/Sample_names.csv",
        mae="config/Samples_MAE.csv"
    output:
        result_mae="results/SNPs_ready_mae.csv", # SNPs that follow monoallelic expression for one tissue
        result_no_mae="results/SNPs_ready_no_mae.csv" # SNPs that do not follow monoallelic expression
    params:
        extra="",
        java_opts="-XX:MinRAMPercentage=80.0 -Xms200G -XX:-UseConcMarkSweepGC -XX:ParallelGCThreads=20 -XX:+UseTLAB",
    threads: 20
    resources:
        mem_mb=200000
    conda:
        "envs/python.yaml"
    log:
        "logs/workflow/mae.log"
    script:
        "scripts/ASE_workflow_MAE.py"



rule valids_no_mae:
    input:
        result1="results/SNPs_ready_no_mae.csv",
    output:
        result4="results/SNPs_ready_no_mae_valid_genome.csv",
        error = "results/Error_in_SNP_calling_no_mae_valid_genome.csv"
    log:
        "logs/python/ASE_data_wrangling/valids.log"
    conda:
        "envs/python.yaml"
    resources:
        mem_mb=200000
    threads:20
    script:
        "scripts/valid.py"



rule QC_no_mae:
    input:
        i1="results/SNPs_ready_no_mae_valid_genome.csv",
        i2="config/Experimental_groups.csv",
        sn1="config/Sample_names.csv"
    output:
        results="results/QC_no_mae_valid_genome.csv"
    conda:
        "envs/python.yaml"
    resources:
        mem_mb=200000
    threads:20
    log:
        "logs/python/QC/QC.log"
    script:
        "scripts/QC.py"
