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
        java_opts=config['java_opts_parallel'],
    threads: config['threads_parallel']
    resources:
        mem_mb=config['mem_mb_parallel']
    script:
        "scripts/splitncigarreads.py"
