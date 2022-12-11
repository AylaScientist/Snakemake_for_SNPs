# Quality control of t 773he output data in the workflow
rule QC:
    input:
        i1="results/SNPs_ready_no_mae.csv",
        i2=config['Experimental_groups'],
        sn1=config['Sample_names']
    output:
        results="results/QC_valid_genome.csv"
    params:
        extra="",
        java_opts=config['java_opts'],
    threads: config['threads']
    resources:
        mem_mb=config['mem_mb']
    conda:
        "envs/python.yaml"
    log:
        "logs/workflow/qc.log"
    script:
        "scripts/qc.py"


rule stats:
    input:
        i1=config['Experimental_design'],
        i2=config['Experimental_groups'],
        i3=config['Sample_names'],
        i4="results/QC_valid_genome.csv"
    output:
        o1="stats/Chi_test_valid_genome.csv",
        o2="stats/Fisher_test_valid_genome.csv",
        o3="results/SNPs_analysed_valid_genome.csv"
    params:
        extra="",
        java_opts=config['java_opts'],
    threads: config['threads']
    resources:
        mem_mb=config['mem_mb']
    conda:
        "envs/python.yaml"
    log:
        "logs/python/stats/SNPs_analysed.log"
    script:
        "scripts/Stats.py"


rule stats_mae:
    input:
        i1=config['Experimental_design'],
        i2=config['Experimental_groups'],
        i3=config['Sample_names'],
        i4="results/QC_mae_valid_genome.csv"
    output:
        o1="stats/Chi_test_mae_valid_genome.csv",
        o2="stats/Fisher_test_mae_valid_genome.csv",
        o3="results/SNPs_analysed_mae_valid_genome.csv"
    params:
        extra="",
        java_opts=config['java_opts'],
    threads: config['threads']
    resources:
        mem_mb=config['mem_mb']
    conda:
        "envs/python.yaml"
    log:
        "logs/python/stats/SNPs_mae_analysed.log"
    script:
        "scripts/Stats.py"
