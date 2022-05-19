rule combine_gvcfs:
    input:
        gvcfs=expand("calls/{sample}_ref.g.vcf", sample = samples),
        ref="genome/Tilapia_header_GCF_001858045.2.fa"
    output:
        gvcf=temp("calls/all_ref_g.vcf")
    conda:
        "envs/gatk.yaml"
    log:
        "logs/gatk/combinegvcfs/combinegvcfs.log"
    params:
        extra = "",
        java_opts="-XX:MinRAMPercentage=80.0 -Xms200G -XX:-UseConcMarkSweepGC -XX:ParallelGCThreads=20 -XX:+UseTLAB",
    threads: 20
    resources:
        mem_mb=200000
    script:
        "scripts/combinegvcfs.py"
