rule haplotype_caller:
    input:
        # single or list of bam files
        bam="split/{sample}_ref.bam",
        ref=config['ref']['genome']
    output:
        gvcf="calls/{sample}_ref.g.vcf"
    conda:
        "envs/gatk.yaml"
    log:
        "logs/gatk/haplotypecaller/{sample}.log"
    params:
        extra="",  # optional
        java_opts=config['java_opts_parallel'],
    threads: config['threads_parallel']
    resources:
        mem_mb=config['mem_mb_parallel']
    script:
        "scripts/haplotypecaller.py"
