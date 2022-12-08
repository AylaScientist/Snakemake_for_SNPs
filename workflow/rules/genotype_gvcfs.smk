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
        java_opts=config['java_opts'],
    threads: config['threads']
    resources:
        mem_mb= config['mem_mb']
    script:
        "scripts/genotypegvcfs.py"
