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
        ref=config['ref']['genome'],
        java_opts=config['java_opts_parallel'],
    threads: config['threads_parallel']
    resources:
        mem_mb=config['mem_mb_parallel']
    log:
        "logs/picard/ValidateSamFile/Validate_bqsr_{sample}.log"
    shell:
        "validatesamfile.py"
