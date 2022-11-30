rule replace_rg_PSG:
    input:
        "{sample}_UMI_{pseudo}_Aligned.sortedByCoord.out.bam"
    output:
        temp("fixed-rg/{sample}_{pseudo}.bam")
    conda:
        "envs/picard.yaml"
    log:
        "logs/picard/replace_rg/{sample}_{pseudo}.log"
    params:
        extra = "SORT_ORDER=coordinate RGID=NextSeq RGLB=idp RGPL=illumina RGPU={sample} RGSM={sample} CREATE_INDEX=True",
        java_opts="-XX:MinRAMPercentage=80.0 -Xms100G -XX:-UseConcMarkSweepGC -XX:ParallelGCThreads=10 -XX:+UseTLAB",
    threads: 10
    resources:
        mem_mb=100000
    script:
        "scripts/readgroups.py"
