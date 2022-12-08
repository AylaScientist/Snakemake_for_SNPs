rule gatk_select_PSG:
    input:
        vcf="calls/all_{pseudo}.vcf",
        ref="pseudogenomes/14FW20-7/14FW20-7_GCF_001858045.2.fa",
    output:
        vcf="calls/selected_{pseudo}.vcf"
    conda:
        "envs/gatk.yaml"
    log:
        "logs/gatk/select/snvs_Nile_{pseudo}.log"
    params:
        extra="--restrict-alleles-to BIALLELIC",  # optional filter arguments, see GATK docs
        java_opts=config['java_opts_combine'],
    threads: config['threads_combine']
    resources:
        mem_mb=config['mem_mb_combine']
    script:
        "scripts/selectvariants.py"
