rule gatk_variantstotable_PSG:
    input:
        o1="annotated/annotated_all_snps_{pseudo}"+config['params']['annotation']['output'],
        ref=config['params']['var2table']
    output:
        vcf="variants/AD_GT_counts_bi_{pseudo}.table"
    conda:
        "envs/gatk.yaml"
    log:
        "logs/gatk/var2table/snvs_{pseudo}.log"
    params:
        o1="annotated/annotated_all_snps_{pseudo}"+config['params']['annotation']['output'],
        extra= config['params']['annotation']['extra'],  # optional filter arguments, see GATK docs
        java_opts=config['java_opts_combine'],
    threads: config['threads_combine']
    resources:
        mem_mb=config['mem_mb_combine']
    script:
        "scripts/variantstotable.py"
