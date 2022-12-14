rule genotype_gvcfs_PSG:
    input:
        gvcf="calls/all_g_{pseudo}.vcf",  # combined gvcf over multiple samples
        ref=config['ref']['genome'],
    output:
        vcf="calls/all_{pseudo}.vcf",
    conda:
        "envs/gatk.yaml"
    log:
        "logs/gatk/genotypegvcfs/genotypegvcfs_{pseudo}.log"
    params:
        extra="",  # optional
        java_opts=config['java_opts_parallel'],
    threads: config['threads_parallel']
    resources:
        mem_mb=config['mem_mb_parallel']
    script:
        "scripts/genotypegvcfs.py"
