rule merge_dataframes:
    input:
        csv = expand("workflow/AD_GT_counts_bi_{pseudo}.csv", pseudo = pseudos)
    output:
        "workflow/merged_df.csv"
    params:
        extra="",
        java_opts="-XX:MinRAMPercentage=80.0 -Xms200G -XX:-UseConcMarkSweepGC -XX:ParallelGCThreads=20 -XX:+UseTLAB",
    threads: 20
    resources:
        mem_mb=200000
    conda:
        "envs/python.yaml"
    log:
        "logs/workflow/merge_dataframes.log"
    script:
        "scripts/mergedataframes.py"


rule biallelic_sites:
    input:
        csv="workflow/merged_df.csv",
        #sn1="config/Sample_names_complete_project3_Nile.csv",
        sn1="config/Sample_names.csv",
        psc="config/Pseudogenome_codes.csv"
    output:
        "workflow/biallelic.csv"
    params:
        extra="",
        java_opts="-XX:MinRAMPercentage=80.0 -Xms200G -XX:-UseConcMarkSweepGC -XX:ParallelGCThreads=20 -XX:+UseTLAB",
    threads: 20
    resources:
        mem_mb=200000
    conda:
        "envs/python.yaml"
    log:
        "logs/workflow/biallelic_sites.log"
    script:
        "scripts/biallelicsites.py"


rule average:
    input:
        csv="workflow/biallelic.csv",
        sn1="config/Sample_names.csv",
        psc="config/Pseudogenome_codes.csv"
    output:
        "workflow/average.csv"
    params:
        extra="",
        java_opts="-XX:MinRAMPercentage=80.0 -Xms200G -XX:-UseConcMarkSweepGC -XX:ParallelGCThreads=20 -XX:+UseTLAB",
    threads: 20
    resources:
        mem_mb=200000
    conda:
        "envs/python.yaml"
    log:
        "logs/workflow/average.log"
    script:
        "scripts/average.py"



rule AD10:
    input:
        csv="workflow/average.csv",
        sn1="config/Sample_names.csv"
    output:
        "workflow/ad_valid_genome.csv"
    params:
        extra="",
        java_opts="-XX:MinRAMPercentage=80.0 -Xms200G -XX:-UseConcMarkSweepGC -XX:ParallelGCThreads=20 -XX:+UseTLAB",
    threads: 20
    resources:
        mem_mb=200000
    conda:
        "envs/python.yaml"
    log:
        "logs/workflow/ad10.log"
    script:
        "scripts/ad10.py"


rule genotype:
    input:
        csv="workflow/ad_valid_genome.csv",
        sn1="config/Sample_names.csv"
    output:
        result1="results/Uniform_to_validate_SNPs_valid_genome.csv",
        result2="results/SNP_panel_valid_genome.csv", #This SNP panel contains the SNPs found expressed in the experimental popilation
        result3="results/Genotype_models_valid_genome.csv"
    params:
        extra="",
        java_opts="-XX:MinRAMPercentage=80.0 -Xms200G -XX:-UseConcMarkSweepGC -XX:ParallelGCThreads=20 -XX:+UseTLAB",
    threads: 20
    resources:
        mem_mb=200000
    conda:
        "envs/python.yaml"
    log:
        "logs/workflow/genotype.log"
    script:
        "scripts/genotype.py"

# Here you can include a step to separate sites that are homozygot for some sample
# and  monoallelic expression (MAE) in case you have the genotypes from the genome
# The monoallelic expression (MAE) can result from homozygots as well as imprinted genes.
# In order to distinguish each case, we need to know the genotype.
