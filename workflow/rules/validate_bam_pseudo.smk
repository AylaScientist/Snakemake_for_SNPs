rule validate_bqsr_PSG:
    input:
        bam="recal/{sample}_{pseudo}.bam",
        recal_table="recal_tables/{sample}_after_{pseudo}.table"
    output:
        "logs/picard/ValidateSamFile_{pseudo}/Validate_bqsr_{sample}_{pseudo}.log"
    conda:
        "envs/gatk.yaml"
    params:
        extra = "",
        mode = "SUMMARY",
        ref=config['ref']['genome'],java_opts=config['java_opts_parallel'],
    threads: config['threads_parallel']
    resources:
        mem_mb=config['mem_mb_parallel']
    log:
        "logs/picard/ValidateSamFile_{pseudo}/Validate_bqsr_{sample}_{pseudo}.log"
    shell:
        "validatesamfile.py"
