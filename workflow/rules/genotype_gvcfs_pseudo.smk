rule genotype_gvcfs_PSG:
    input:
        gvcf="calls/all_g_{pseudo}.vcf",  # combined gvcf over multiple samples
        ref="pseudogenomes/14FW20-7/14FW20-7_GCF_001858045.2.fa",
    output:
        vcf="calls/all_{pseudo}.vcf",
    conda:
        "envs/gatk.yaml"
    log:
        "logs/gatk/genotypegvcfs/genotypegvcfs_{pseudo}.log"
    params:
        extra="",  # optional
        java_opts="-XX:MinRAMPercentage=80.0 -Xms200G -XX:-UseConcMarkSweepGC -XX:ParallelGCThreads=20 -XX:+UseTLAB",
    threads: 20
    resources:
        mem_mb=200000
    script:
        "scripts/genotypegvcfs.py"
