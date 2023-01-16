rule table_step1_PSG:
    input:
        table1="variants/AD_GT_counts_bi_{pseudo}.table"
    output:
        table1="variants/AD_GT_counts_bi_{pseudo}_step1.csv"
    params:
        char = "characters{pseudo}.table",
        multi = "delete_multiallelic{pseudo}.table"
    run:
        shell("sed 's/|/\//g' {input.table1} > {params.char} ")
        shell("awk 'seen[$1,$2]++' {params.char} > {params.multi}") #Collect SNPs with same chromosome and position thus are multiallelic by duplication of CHROM and POS columns 1 and 2
        shell("awk -F'\t' 'NR==FNR{{a[$1,$2]++;next}} !(a[$1,$2])' {params.multi} {params.char}  > {output.table1} ") # Eliminate the rows collected as multiallelic previously
        shell("rm {params.char} {params.multi}") # Eliminate intermediary files



rule table_step2_PSG:
    input:
        table="variants/AD_GT_counts_bi_{pseudo}_step1.csv",
        tb_colnames = "../config/tb{pseudo}_colnames.csv",
    output:
        csv="workflow/AD_GT_counts_bi_{pseudo}.csv"
    resources:
        mem_mb=config['mem_mb_combine']
    threads: config['threads_combine']
    conda:
        "envs/python.yaml"
    log:
        "logs/python/table_step2/AD_GT_counts_bi_{pseudo}.log"
    script:
        "scripts/table2df_step2.py"
