rule gatk_variantstotable_PSG:
    input:
        vcf="annotated/annotated_all_snps_{pseudo}.ON_multianno.vcf",
        ref="pseudogenomes/{pseudo}/{pseudo}_GCF_001858045.2.fa"
    output:
        vcf="variants/AD_GT_counts_bi_{pseudo}.table"
    conda:
        "envs/gatk.yaml"
    log:
        "logs/gatk/var2table/snvs_Nile_{pseudo}.log"
    params:
        o1= config['params']['annotation']['o1'],
        extra= config['params']['annotation']['extra'],  # optional filter arguments, see GATK docs
        java_opts=config['java_opts_combine'],
    threads: config['threads_combine']
    resources:
        mem_mb=config['mem_mb_combine']
    script:
        "scripts/variantstotable.py"
