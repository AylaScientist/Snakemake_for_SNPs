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
        java_opts=config['java_opts'],
    threads: config['threads']
    resources:
        mem_mb=config['mem_mb']
    script:
        "scripts/combinegvcfs.py"
