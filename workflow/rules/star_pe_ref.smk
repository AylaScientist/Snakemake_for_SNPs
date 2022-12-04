rule star_pe:
    input:
        # use a list for multiple fastq files for one sample
        # usually technical replicates across lanes/flowcells example: "reads/{sample}._R2.1.fastq", "reads/{sample}._R2.2.fastq"
        fastq1="trimmed/{sample}.1.fastq",
        # paired end reads needs to be ordered so each item in the two lists match
        fastq2="trimmed/{sample}.2.fastq",
        index="genome/Tilapia_header_GCF_001858045.2.dict",
        params="genome/SAindex"
    output:
        # see STAR manual for additional output files
        file=temp("{sample}__ref_Aligned.sortedByCoord.out.bam") #""star/{sample}/Aligned.sortedByCoord.out.sam""
    log:
        "logs/star/{sample}__ref.log"
    params:
        # path to STAR reference genome index
        # optional parameters
        # "-Xmx100g -XX:ConcGCThreads=4 -XX:ParallelGCThreads=4" #""#"-Xmx180g -XX:ParallelGCThreads=18"  # For Java memory specifications, please only use resources.mem_mb.
        extra="",
        filename="{sample}__ref_",
        threads="10",
        java_opts="-XX:MinRAMPercentage=80.0 -Xms100G -XX:-UseConcMarkSweepGC -XX:ParallelGCThreads=10 -XX:+UseTLAB",
    threads: 10
    resources:
        mem_mb=100000
    conda:
        "envs/star.yaml"
    script:
        "scripts/star.py"
