rule merge_dataframes:
    input:
        csv = expand("workflow/AD_GT_counts_bi_{pseudo}.csv", pseudo = pseudos)
    output:
        "workflow/merged_df.csv"
    params:
        extra="",
        java_opts=config['java_opts'],
    threads: config['threads']
    resources:
        mem_mb=config['mem_mb']
    conda:
        "envs/python.yaml"
    log:
        "logs/workflow/merge_dataframes.log"
    script:
        "scripts/mergedataframes.py"


rule biallelic_sites:
    input:
        csv="workflow/merged_df.csv",
        sn1=config['Sample_names'],
        psc=config['Pseudogenome_codes']
    output:
        "workflow/biallelic.csv"
    params:
        extra="",
        java_opts=config['java_opts'],
    threads: config['threads']
    resources:
        mem_mb=config['mem_mb']
    conda:
        "envs/python.yaml"
    log:
        "logs/workflow/biallelic_sites.log"
    script:
        "scripts/biallelicsites.py"


rule average:
    input:
        csv="workflow/biallelic.csv",
        sn1=config['Sample_names'],
        psc=config['Pseudogenome_codes']
    output:
        "workflow/average.csv"
    params:
        extra="",
        java_opts=config['java_opts'],
    threads: config['threads']
    resources:
        mem_mb=config['mem_mb']
    conda:
        "envs/python.yaml"
    log:
        "logs/workflow/average.log"
    script:
        "scripts/average.py"



rule AD10:
    input:
        csv="workflow/average.csv",
        sn1=config['Sample_names']
    output:
        "workflow/ad_valid_genome.csv"
    params:
        extra="",
        java_opts=config['java_opts'],
    threads: config['threads']
    resources:
        mem_mb=config['mem_mb']
    conda:
        "envs/python.yaml"
    log:
        "logs/workflow/ad10.log"
    script:
        "scripts/ad10.py"


rule genotype:
    input:
        csv="workflow/ad_valid_genome.csv",
        sn1=config['Sample_names']
    output:
        result1="results/Uniform_to_validate_SNPs_valid_genome.csv",
        result2="results/SNP_panel_valid_genome.csv", #This SNP panel contains the SNPs found expressed in the experimental popilation
        result3="results/Genotype_models_valid_genome.csv"
    params:
        extra="",
        java_opts=config['java_opts'],
    threads: config['threads']
    resources:
        mem_mb=config['mem_mb']
    conda:
        "envs/python.yaml"
    log:
        "logs/workflow/genotype.log"
    script:
        "scripts/genotype.py"

rule mae:
    input:
        result1="results/Uniform_to_validate_SNPs_valid_genome.csv",
        names=config['Sample_names'],
        mae=config['Samples_MAE']
    output:
        result_mae="results/SNPs_ready_mae.csv", # SNPs that follow monoallelic expression for one tissue
        result_no_mae="results/SNPs_ready_no_mae.csv" # SNPs that do not follow monoallelic expression
    params:
        extra="",
        java_opts=config['java_opts'],
    threads: config['threads']
    resources:
        mem_mb=config['mem_mb']
    conda:
        "envs/python.yaml"
    log:
        "logs/workflow/mae.log"
    script:
        "scripts/ASE_workflow_MAE.py"




# In case you have the SNPs from the genotype, here you will calculate the error

rule valids:
    input:
        result1="results/Uniform_to_validate_SNPs_valid_genome.csv",
    output:
        result4="results/SNPs_ready_valid_genome.csv",
        error = "results/Error_in_SNP_calling_valid_genome.csv"
    params:
        extra="",
        java_opts=config['java_opts'],
    threads: config['threads']
    resources:
        mem_mb=config['mem_mb']
    conda:
        "envs/python.yaml"
    log:
        "logs/python/ASE_data_wrangling/valids.log"
    script:
        "scripts/valid.py"
