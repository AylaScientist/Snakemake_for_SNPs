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
    params:
        extra="",  # optional
        java_opts="-XX:MinRAMPercentage=80.0 -Xms100G -XX:-UseConcMarkSweepGC -XX:ParallelGCThreads=10 -XX:+UseTLAB",
    threads: 10
    resources:
        mem_mb=100000
    script:
        "scripts/haplotypecaller.py"
