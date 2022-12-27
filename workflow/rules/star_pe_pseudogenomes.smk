rule star_pe_PSG:
    input:
        # use a list for multiple fastq files for one sample
        # usually technical replicates across lanes/flowcells example: "reads/{sample}._R2.1.fastq", "reads/{sample}._R2.2.fastq"
        fastq1="trimmed/{sample}.1.fastq",
        # paired end reads needs to be ordered so each item in the two lists match
        fastq2="trimmed/{sample}.2.fastq",
        index="pseudogenomes/{pseudo}/{pseudo}_"+config['ref']['release']+".dict",
        index2="pseudogenomes/{pseudo}/SAindex"
    output:
        # see STAR manual for additional output files
        file=temp("{sample}_{pseudo}_Aligned.sortedByCoord.out.bam") #""star/{sample}/Aligned.sortedByCoord.out.sam""
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
        "scripts/star.py"