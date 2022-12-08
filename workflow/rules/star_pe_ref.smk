rule star_pe:
    input:
        # use a list for multiple fastq files for one sample
        # usually technical replicates across lanes/flowcells example: "reads/{sample}._R2.1.fastq", "reads/{sample}._R2.2.fastq"
        fastq1="trimmed/{sample}.1.fastq",
        # paired end reads needs to be ordered so each item in the two lists match
        fastq2="trimmed/{sample}.2.fastq", #optional
        index=config['ref']['dict'],
        params=config['ref']['SAindex']
    output:
        # see STAR manual for additional output files
        file=temp("{sample}_ref_Aligned.sortedByCoord.out.bam") #""star/{sample}/Aligned.sortedByCoord.out.sam""
    log:
        "logs/star/{sample}_ref.log"
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
