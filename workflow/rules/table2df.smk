rule table_step1_PSG:
    input:
        table1="variants/AD_GT_counts_bi_{pseudo}.table"
    output:
        table1="variants/AD_GT_counts_bi_{pseudo}_step1.csv"
    run:
        shell("sed 's/|/\//g' {input.table1} > characters.table ")
        shell("awk 'seen[$1,$2]++' characters.table > delete_multiallelic.table") #Collect SNPs with same chromosome and position thus are multiallelic by duplication of CHROM and POS columns 1 and 2
        #shell("awk '$7~/^(\d+,\d+,\d+)*$/' characters.table > delete_multiallelic.table")
        shell("awk '{print $1 $2}' delete_multiallelic.table > coordinates_multiallelic.table")
        shell("grep -Fvf coordinates_multiallelic.table characters.table  > {output.table1} ") #Eliminate multiallelics, will output contents in file1 not in file2 
        #shell("rm characters.table delete_multiallelic.table")




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
