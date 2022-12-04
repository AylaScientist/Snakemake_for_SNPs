rule validate_bqsr:
    input:
        bam="split/{sample}_ref.bam"
    output:
        "logs/picard/ValidateSamFile/Validate_bqsr_{sample}.log"
    conda:
        "envs/gatk.yaml"
    params:
        extra = "",
        mode = "SUMMARY",
        ref="genome/Tilapia_header_GCF_001858045.2.fa",
        java_opts="-XX:MinRAMPercentage=80.0 -Xms100G -XX:-UseConcMarkSweepGC -XX:ParallelGCThreads=10 -XX:+UseTLAB",
    threads: 10
    resources:
        mem_mb=100000
    log:
        "logs/picard/ValidateSamFile/Validate_bqsr_{sample}.log"
    shell:
        "validatesamfile.py"
