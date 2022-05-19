rule validate_bqsr_PSG:
    input:
        bam="recal/{sample}_{pseudo}.bam",
        recal_table="recal_tables/{sample}_after_{pseudo}.table"
    output:
        "logs/picard/ValidateSamFile_{pseudo}/Validate_bqsr_{sample}_{pseudo}.log"
    conda:
        "envs/gatk.yaml"
    params:
        extra = "",
        mode = "SUMMARY",
        ref="genome/Tilapia_header_GCF_001858045.2.fa"
    params:
        extra="",  # optional
        java_opts="-XX:MinRAMPercentage=80.0 -Xms100G -XX:-UseConcMarkSweepGC -XX:ParallelGCThreads=10 -XX:+UseTLAB",
    threads: 10
    resources:
        mem_mb=100000
    log:
        "logs/picard/ValidateSamFile_{pseudo}/Validate_bqsr_{sample}_{pseudo}.log"
    shell:
        "validatesamfile.py"
