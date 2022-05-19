rule trimmomatic_pe:
    input:
        r1="fastq_merged/{sample}_R1.fastq",
        r2="fastq_merged/{sample}_R3.fastq"
    output:
        r1=("trimmed/{sample}.1.fastq"),
        r2=("trimmed/{sample}.2.fastq"),
        # reads where trimming entirely removed the mate
        r1_unpaired=temp("trimmed/{sample}.1.unpaired.fastq"),
        r2_unpaired=temp("trimmed/{sample}.2.unpaired.fastq")
    conda:
        "envs/trimmomatic.yaml"
    log:
        "logs/trimmomatic/{sample}.log"
    params:
        # list of trimmers (see manual)
        trimmer="ILLUMINACLIP:/home/fish_station/anaconda3/envs/snakemake/share/trimmomatic-0.39-2/adapters/TruSeq3-PE-2.fa:2:30:10:2:keepBothReads LEADING:3 TRAILING:3 MINLEN:38", # optional parameters
        #"ILLUMINACLIP:TruSeq3-PE-2.fa:2:30:10:2:keepBothReads SLIDINGWINDOW:4:5 LEADING:10 TRAILING:10 MINLEN:38",
        extra="",
        java_opts="-XX:MinRAMPercentage=80.0 -XX:-UseConcMarkSweepGC -XX:ParallelGCThreads=10 -XX:+UseTLAB",
    threads: 10
    resources:
        mem_mb=100000
    script:
        "scripts/trimmomatic_pe.py"
