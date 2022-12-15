rule ref_index:
    input:
        config['ref']['genome']
    output:
        config['ref']['SAindex'],
    params:
        threads= config['threads'],
        dir = "genome/",
        annotation = config['ref']['annotation'],
        read_length = config['params']['star']['read_length'], #Read length -1
        java_opts=config['java_opts'],
    threads: config['threads']
    resources:
        mem_mb=config['mem_mb']
    conda:
        "envs/star.yaml"
    log:
        "logs/genomes/REF_index.log"
    script:
        "scripts/star_gi.py"

rule create_dict:
    input:
        config['ref']['genome'],
    output:
        config['ref']['dict'],
    log:
        "logs/picard/create_dict.log",
    params:
        extra="",  
        java_opts=config['java_opts'],
    threads: config['threads']
    resources:
        mem_mb=config['mem_mb']
    conda:
        "envs/picard.yaml"
    script:
        "scripts/createsequencedictionary.py"