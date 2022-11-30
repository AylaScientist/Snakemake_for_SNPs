rule ref_index:
    input:
        "genome/Tilapia_header_GCF_001858045.2.fa"
    output:
        "genome/SAindex"
    params:
        threads= "10",
        dir = "genome/",
        annotation = "genome/Tilapia_GCF_001858045.2_annotation.gtf",
        read_length = "148", #Read length -1
        java_opts="-XX:MinRAMPercentage=80.0 -Xms100G -XX:-UseConcMarkSweepGC -XX:ParallelGCThreads=10 -XX:+UseTLAB",
    threads: 10
    resources:
        mem_mb=100000
    conda:
        "envs/star.yaml"
    log:
        "logs/genomes/REF_index.log"
    script:
        "scripts/star_gi.py"
