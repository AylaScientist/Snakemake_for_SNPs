rule mark_duplicates:
    input:
        "fixed-rg/{sample}.bam"
    output:
        bam=temp("marked_dedup/{sample}_ref.bam"),
        metrics="marked_dedup/{sample}_ref.metrics.txt"
    conda:
        "envs/picard.yaml"
    log:
        "logs/picard/dedup/{sample}.log"
    params:
        extra = "REMOVE_DUPLICATES=true", #Duplicates can also be removed with UMI-tools
        java_opts="-XX:MinRAMPercentage=80.0 -Xms100G -XX:-UseConcMarkSweepGC -XX:ParallelGCThreads=10 -XX:+UseTLAB",
    threads: 10
    resources:
        mem_mb=100000
    script:
        "scripts/mark_duplicates.py"
