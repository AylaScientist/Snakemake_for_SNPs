rule combine_gvcfs_PSG:
    input:
        gvcfs=expand("calls/{sample}_{pseudo}.g.vcf", sample = samples, pseudo = pseudos),
        ref="pseudogenomes/14FW20-7/14FW20-7_GCF_001858045.2.fa",
    output:
        gvcf=temp("calls/all_g_{pseudo}.vcf")
    conda:
        "envs/gatk.yaml"
    log:
        "logs/gatk/combinegvcfs/combinegvcfs_{pseudo}.log"
    params:
        extra = "",
    params:
        extra="",  # optional
        java_opts="-XX:MinRAMPercentage=80.0 -Xms200G -XX:-UseConcMarkSweepGC -XX:ParallelGCThreads=20 -XX:+UseTLAB",
    threads: 20
    resources:
        mem_mb=200000
    script:
        "scripts/combinegvcfs.py"
