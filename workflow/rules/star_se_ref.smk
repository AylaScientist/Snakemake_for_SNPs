rule star_se:
    input:
        fastq1="trimmed/{sample}.1.fastq",
        # path to STAR reference genome index
        index=config['ref']['dict'],
        params=config['ref']['SAindex']idx="index",
    output:
        # see STAR manual for additional output files
        file=temp("{sample}_ref_Aligned.sortedByCoord.out.bam") #""star/{sample}/Aligned.sortedByCoord.out.sam""
    log:
        "logs/star/{sample}_ref.log",
    params:
        # path to STAR reference genome index
        # optional parameters
        # For Java memory specifications, please only use resources.mem_mb.
        extra="",
        filename="{sample}_ref_",
        threads=config['threads_parallel'],
        java_opts=config['java_opts_parallel'],
    threads: config['threads_parallel']
    resources:
        mem_mb=config['mem_mb_parallel']
    conda:
        "envs/star.yaml"
    script:
        "scripts/star.py"
