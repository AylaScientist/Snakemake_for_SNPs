# Mark and eliminate Duplicates:
rule mark_duplicates_pseudo:
    input:
        "fixed-rg/{sample}_{pseudo}.bam"
    output:
        bam=temp("marked_dedup/{sample}_{pseudo}.bam"),
        metrics="marked_dedup/{sample}_{pseudo}.metrics.txt"
    conda:
        "envs/picard.yaml"
    log:
        "logs/picard/dedup/{sample}_{pseudo}.log"
    params:
        extra = "REMOVE_DUPLICATES=true", #Duplicates can also be removed with UMI-tools
        java_opts="-XX:MinRAMPercentage=80.0 -Xms100G -XX:-UseConcMarkSweepGC -XX:ParallelGCThreads=10 -XX:+UseTLAB",
    threads: 10
    resources:
        mem_mb=100000
    script:
        "scripts/mark_duplicates.py"
