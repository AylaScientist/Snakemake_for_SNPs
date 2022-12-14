rule combine_gvcfs:
    input:
        gvcfs=expand("calls/{sample}_ref.g.vcf", sample = samples),
        ref=config['ref']['genome'],
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
