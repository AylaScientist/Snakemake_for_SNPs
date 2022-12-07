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
        java_opts=config['java_opts_parallel'],
    threads: config['threads_parallel']
    resources:
        mem_mb=config['mem_mb_parallel']
    script:
        "scripts/mark_duplicates.py"
