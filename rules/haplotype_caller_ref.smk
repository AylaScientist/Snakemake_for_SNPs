rule haplotype_caller:
    input:
        # single or list of bam files
        bam="recal/{sample}_ref.bam",
        ref="genome/Tilapia_header_GCF_001858045.2.fa"
    output:
        gvcf="calls/{sample}_ref.g.vcf"
    conda:
        "envs/gatk.yaml"
    log:
        "logs/gatk/haplotypecaller/{sample}.log"
    params:
        extra="",  # optional
        java_opts="-XX:MinRAMPercentage=80.0 -Xms100G -XX:-UseConcMarkSweepGC -XX:ParallelGCThreads=10 -XX:+UseTLAB",
    threads: 10
    resources:
        mem_mb=100000
    script:
        "scripts/haplotypecaller.py"
