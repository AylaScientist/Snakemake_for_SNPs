rule star_se_PSG:
    input:
        fastq1="trimmed/{sample}.1.fastq",
        # path to STAR reference genome index
        index="pseudogenomes/{pseudo}/{pseudo}_"+config['ref']['release']+".dict",
        index2="pseudogenomes/{pseudo}/SAindex"
    output:
        # see STAR manual for additional output files
        file="{sample}_{pseudo}_Aligned.sortedByCoord.out.bam" #""star/{sample}/Aligned.sortedByCoord.out.sam""
    log:
        "logs/star/{sample}_{pseudo}.log"
    params:
        # path to STAR reference genome index
        # optional parameters
        extra="",
        filename="{sample}_{pseudo}_",
        threads=config['threads_parallel'],
        java_opts=config['java_opts_parallel'],
    threads: config['threads_parallel']
    resources:
        mem_mb=config['mem_mb_parallel']
    conda:
        "envs/star.yaml"
    script:
        "scripts/star_se.py"

