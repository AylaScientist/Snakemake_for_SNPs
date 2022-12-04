rule table_step1_PSG:
    input:
        table1="variants/AD_GT_counts_bi_{pseudo}.table"
    output:
        table1="variants/AD_GT_counts_bi_{pseudo}_step1.csv"
    shell:
        "sed 's/|/\//g' {input.table1} > {output.table1}"


rule table_step2_PSG:
    input:
        table="variants/AD_GT_counts_bi_{pseudo}_step1.csv",
        tb_colnames = "complementary/tb1_colnames.csv",
    output:
        csv="workflow/AD_GT_counts_bi_{pseudo}.csv"
    resources:
        mem_mb=200000
    threads: 20
    conda:
        "envs/python.yaml"
    log:
        "logs/python/table_step2/AD_GT_counts_bi_{pseudo}.log"
    script:
        "scripts/table2df_step2.py"
