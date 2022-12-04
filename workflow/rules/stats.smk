rule stats:
    input:
        i1="complementary/Experimental_design.csv",
        i2="complementary/Experimental_groups.csv",
        i3="complementary/Sample_names.csv",
        i4="results/QC_no_mae_valid_genome.csv"
    output:
        o1="stats/Chi_test_valid_genome.csv",
        o2="stats/Fisher_test_valid_genome.csv",
        o3="results/SNPs_analysed_valid_genome.csv"
    log:
        "logs/python/stats/SNPs_analysed.log"
    conda:
        "envs/python.yaml"
    resources:
        mem_mb=200000
    threads:20
    script:
        "scripts/Stats.py"

rule stats_mae:
    input:
        i1="complementary/Experimental_design.csv",
        i2="complementary/Experimental_groups.csv",
        i3="complementary/Sample_names.csv",
        i4="results/QC_mae_valid_genome.csv"
    output:
        o1="stats/Chi_test_mae_valid_genome.csv",
        o2="stats/Fisher_test_mae_valid_genome.csv",
        o3="results/SNPs_analysed_mae_valid_genome.csv"
    log:
        "logs/python/stats/SNPs_mae_analysed.log"
    conda:
        "envs/python.yaml"
    resources:
        mem_mb=200000
    threads:20
    script:
        "scripts/Stats.py"
