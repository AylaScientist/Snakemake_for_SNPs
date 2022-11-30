rule star_pe_PSG1:
    input:
        # use a list for multiple fastq files for one sample
        # usually technical replicates across lanes/flowcells example: "reads/{sample}._R2.1.fastq", "reads/{sample}._R2.2.fastq"
        fastq1="trimmed/{sample}.1.fastq",
        # paired end reads needs to be ordered so each item in the two lists match
        fastq2="trimmed/{sample}.2.fastq",
        index="pseudogenomes/14FW20-7/14FW20-7_GCF_001858045.2.dict",
        index2="pseudogenomes/14FW20-7/SAindex"
    output:
        # see STAR manual for additional output files
        file=temp("{sample}_UMI_PSG1_Aligned.sortedByCoord.out.bam") #""star/{sample}/Aligned.sortedByCoord.out.sam""
    log:
        "logs/star/{sample}_UMI_PSG1.log"
    params:
        # path to STAR reference genome index
        # optional parameters
        extra="",
        filename="{sample}_UMI_PSG1_",
        threads="10",
        java_opts="-XX:MinRAMPercentage=80.0 -Xms100G -XX:-UseConcMarkSweepGC -XX:ParallelGCThreads=10 -XX:+UseTLAB",
    threads: 10
    resources:
        mem_mb=100000
    conda:
        "envs/star.yaml"
    script:
        "scripts/star.py"

rule star_pe_PSG2:
    input:
        # use a list for multiple fastq files for one sample
        # usually technical replicates across lanes/flowcells example: "reads/{sample}._R2.1.fastq", "reads/{sample}._R2.2.fastq"
        fastq1="trimmed/{sample}.1.fastq",
        # paired end reads needs to be ordered so each item in the two lists match
        fastq2="trimmed/{sample}.2.fastq",
        index="pseudogenomes/4SW2-9/4SW2-9_GCF_001858045.2.dict",
        index1="pseudogenomes/4SW2-9/SAindex"
    output:
        # see STAR manual for additional output files
        file=temp("{sample}_UMI_PSG2_Aligned.sortedByCoord.out.bam") #""star/{sample}/Aligned.sortedByCoord.out.sam""
    log:
        "logs/star/{sample}_UMI_PSG2.log"
    params:
        # path to STAR reference genome index
        # optional parameters
        extra="",
        filename="{sample}_UMI_PSG2_",
        threads="10",
        java_opts="-XX:MinRAMPercentage=80.0 -Xms100G -XX:-UseConcMarkSweepGC -XX:ParallelGCThreads=10 -XX:+UseTLAB",
    threads: 10
    resources:
        mem_mb=100000
    conda:
        "envs/star.yaml"
    script:
        "scripts/star.py"
