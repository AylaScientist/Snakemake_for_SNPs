rule htseq:
    input:
        bam = "split/{sample}_ref.bam",
        annotation ="genome/Tilapia_GCF_001858045.2_annotation.gff3"
    output:
        config["path"]+"gene_counts/{sample}.count_ref.txt"
    params:
        extra=" -f bam -r pos -s no -i Parent", # -n is the number of samples that process in parallel, may need to change in every experimental design.
        java_opts=config['java_opts'],
    threads: config['threads']
    resources:
        mem_mb=config['mem_mb']
    conda:
        "envs/htseq.yaml"
    log:
        "logs/htseq/htseq_{sample}.log"
    script:
        "scripts/htseq.py"
