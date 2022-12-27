#configfile: '../config/config.yaml'

rule combine_gvcfs_PSG0:
    input:
        gvcfs=expand("calls/{sample}_"+pseudos[0]+".g.vcf", sample = samples),
        ref=config['ref']['genome'],
    output:
        gvcf="calls/all_g_"+pseudos[0]+".vcf"
    conda:
        "envs/gatk.yaml"
    log:
        "logs/gatk/combinegvcfs/combinegvcfs_"+pseudos[0]+".log"
    params:
        extra="",  # optional
        java_opts=config['java_opts_combine'],
    threads: config['threads_combine']
    resources:
        mem_mb=config['mem_mb_combine']
    script:
        "scripts/combinegvcfs.py"

rule combine_gvcfs_PSG1:
    input:
        gvcfs=expand("calls/{sample}_"+pseudos[1]+".g.vcf", sample = samples),
        ref=config['ref']['genome'],
    output:
        gvcf="calls/all_g_"+pseudos[1]+".vcf"
    conda:
        "envs/gatk.yaml"
    log:
        "logs/gatk/combinegvcfs/combinegvcfs_"+pseudos[0]+".log"
    params:
        extra="",  # optional
        java_opts=config['java_opts_combine'],
    threads: config['threads_combine']
    resources:
        mem_mb=config['mem_mb_combine']
    script:
        "scripts/combinegvcfs.py"
