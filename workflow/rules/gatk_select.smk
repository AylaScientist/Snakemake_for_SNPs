rule gatk_select:
    input:
        vcf="calls/all_ref_g.vcf",
        ref="genome/Tilapia_header_GCF_001858045.2.fa"
    output:
        vcf="calls/selected_ref_RNA.vcf"
    conda:
        "envs/gatk.yaml"
    log:
        "logs/gatk/select/snvs_Nile.log"
    params:
        extra="--restrict-alleles-to BIALLELIC",  # optional filter arguments, see GATK docs
        java_opts="-XX:MinRAMPercentage=80.0 -Xms200G -XX:-UseConcMarkSweepGC -XX:ParallelGCThreads=20 -XX:+UseTLAB",
    threads: 20
    resources:
        mem_mb=200000
    script:
        "scripts/selectvariants.py"
