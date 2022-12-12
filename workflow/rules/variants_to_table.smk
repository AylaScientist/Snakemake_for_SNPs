rule gatk_variantstotable_PSG:
    input:
        vcf=config['params']['annotation']['o1'],
        ref=config['params']['var2table']
    output:
        vcf="variants/AD_GT_counts_bi_{pseudo}.table"
    conda:
        "envs/gatk.yaml"
    log:
        "logs/gatk/var2table/snvs_{pseudo}.log"
    params:
        o1= config['params']['annotation']['o1'],
        extra= config['params']['annotation']['extra'],  # optional filter arguments, see GATK docs
        java_opts=config['java_opts_combine'],
    threads: config['threads_combine']
    resources:
        mem_mb=config['mem_mb_combine']
    script:
        "scripts/variantstotable.py"
