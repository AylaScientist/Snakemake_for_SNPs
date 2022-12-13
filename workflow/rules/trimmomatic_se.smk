rule trimmomatic_se:
    input:
        r1=config["input_files"]["r1"],
    output:
        r1="trimmed/{sample}.1.fastq",
    conda:
        "envs/trimmomatic.yaml"
    log:
        "logs/trimmomatic/{sample}.log"
    params:
        path = config['params']['trimmomatic']['path'], 
        # list of trimmers (see manual)
        trimmer=config['params']['trimmomatic']['se']['trimmer'],
        extra="",
        java_opts=config['java_opts_parallel']
    threads: config['threads_parallel']
    resources:
        mem_mb=config['mem_mb_parallel']
    wrapper:
        "scripts/trimmomatic_se.py"