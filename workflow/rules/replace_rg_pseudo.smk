rule replace_rg_PSG:
    input:
        "{sample}_{pseudo}_Aligned.sortedByCoord.out.bam"
    output:
        temp("fixed-rg/{sample}_{pseudo}.bam")
    conda:
        "envs/picard.yaml"
    log:
        "logs/picard/replace_rg/{sample}_{pseudo}.log"
    params:
        extra = "SORT_ORDER=coordinate RGID=NextSeq RGLB=idp RGPL=illumina RGPU={sample} RGSM={sample} CREATE_INDEX=True",
        java_opts=config['java_opts_parallel'],
    threads: config['threads_parallel']
    resources:
        mem_mb=config['mem_mb_parallel']
    script:
        "scripts/readgroups.py"
