rule combine_gvcfs_PSG:
    input:
        gvcfs=expand("calls/{sample}_{pseudo}.g.vcf", sample = samples, pseudo = pseudos),
        ref=config['ref']['genome'],
    output:
        gvcf=temp("calls/all_g_{pseudo}.vcf")
    conda:
        "envs/gatk.yaml"
    log:
        "logs/gatk/combinegvcfs/combinegvcfs_{pseudo}.log"
    params:
        extra = "",
    params:
        extra="",  # optional
        java_opts=config['java_opts_combine'],
    threads: config['threads_combine']
    resources:
        mem_mb=config['mem_mb_combine']
    script:
        "scripts/combinegvcfs.py"
