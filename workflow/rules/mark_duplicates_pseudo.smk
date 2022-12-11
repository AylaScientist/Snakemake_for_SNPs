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
        java_opts=config['java_opts_parallel'],
    threads: config['threads_parallel']
    resources:
        mem_mb=config['mem_mb_parallel']
    script:
        "scripts/mark_duplicates.py"
