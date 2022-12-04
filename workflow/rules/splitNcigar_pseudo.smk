rule splitncigarreads_PSG:
    input:
        bam="marked_dedup/{sample}_{pseudo}.bam",
        ref="genome/Tilapia_header_GCF_001858045.2.fa"
    output:
        bam=temp("split/{sample}_{pseudo}.bam")
    conda:
        "envs/gatk.yaml"
    log:
        "logs/gatk/splitNCIGARreads/{sample}_{pseudo}.log"
    params:
        extra="",  # optional
        java_opts="-XX:MinRAMPercentage=80.0 -Xms100G -XX:-UseConcMarkSweepGC -XX:ParallelGCThreads=10 -XX:+UseTLAB",
    threads: 10
    resources:
        mem_mb=100000
    script:
        "scripts/splitncigarreads.py"
