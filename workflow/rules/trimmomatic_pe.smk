rule trimmomatic_pe:
    input:
        r1=config['input_files']['r1'],
        r2=config['input_files']['r2']
    output:
        r1=("trimmed/{sample}.1.fastq"),
        r2=("trimmed/{sample}.2.fastq"),
        # reads where trimming entirely removed the mate
        r1_unpaired=temp("trimmed/{sample}.1.unpaired.fastq"),
        r2_unpaired=temp("trimmed/{sample}.2.unpaired.fastq")
    conda:
        "envs/trimmomatic.yaml"
    log:
        "logs/trimmomatic/{sample}.log"
    params:
        path = config['params']['trimmomatic']['path'], 
        # list of trimmers (see manual)
        trimmer=config['params']['trimmomatic']['pe']['trimmer'],
        extra="",
        java_opts=config['java_opts_parallel']
    threads: config['threads_parallel']
    resources:
        mem_mb=config['mem_mb_parallel']
    script:
        "scripts/trimmomatic_pe.py"
