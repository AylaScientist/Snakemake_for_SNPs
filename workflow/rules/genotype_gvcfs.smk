rule genotype_gvcfs:
    input:
        gvcf="calls/all_ref_g.vcf",  # combined gvcf over multiple samples
        ref="genome/Tilapia_header_GCF_001858045.2.fa"
    output:
        vcf="calls/all_ref.vcf",
    conda:
        "envs/gatk.yaml"
    log:
        "logs/gatk/genotypegvcfs/genotypegvcfs.log"
    params:
        extra="",  # optional
        java_opts="-XX:MinRAMPercentage=80.0 -Xms200G -XX:-UseConcMarkSweepGC -XX:ParallelGCThreads=20 -XX:+UseTLAB",
    threads: 20
    resources:
        mem_mb=200000
    script:
        "scripts/genotypegvcfs.py"
