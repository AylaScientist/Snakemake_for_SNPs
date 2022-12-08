rule haplotype_caller_PSG:
    input:
        # single or list of bam files
        bam="recal/{sample}_{pseudo}.bam",
        ref="pseudogenomes/14FW20-7/14FW20-7_GCF_001858045.2.fa",
    output:
        gvcf=temp("calls/{sample}_{pseudo}.g.vcf")
    conda:
        "envs/gatk.yaml"
    log:
        "logs/gatk/haplotypecaller/{sample}_{pseudo}.log"
    params:
        extra="",  # optional
        java_opts=config['java_opts_parallel'],
    threads: config['threads_parallel']
    resources:
        mem_mb=config['mem_mb_parallel']
    script:
        "scripts/haplotypecaller.py"
