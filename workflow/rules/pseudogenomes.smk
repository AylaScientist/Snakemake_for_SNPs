rule pseudogenome:
    input:
        vcf="pseudogenomes/subset_vcf_file_wo_indels.recode.vcf.gz",
        ref="genome/Tilapia_header_GCF_001858045.2.fa",
        index="pseudogenomes/subset_vcf_file_wo_indels.recode.vcf.gz.tbi"
    output:
        pseudo = "pseudogenomes/{pseudo}/{pseudo}_GCF_001858045.2.fa"
    params:
        extra = "{pseudo}"
    resources:
        mem_mb=config['mem_mb_combine']
    threads: config['threads_combine']
    conda:
        "envs/bcftools.yaml"
    log:
        "logs/pseudogenomes/{pseudo}.log"
    script:
        "scripts/pseudogenomes.py"



# Make index for the pseudogenomes for further mapping
rule pseudo_index:
    input:
        "pseudogenomes/{pseudo}/{pseudo}_GCF_001858045.2.fa"
    output:
        "pseudogenomes/{pseudo}/SAindex"
    params:
        threads= "16",
        dir = "pseudogenomes/{pseudo}/",
        annotation = "genome/Tilapia_GCF_001858045.2_annotation.gtf",
        read_length = "148" #Read length -1
    resources:
        mem_mb=config['mem_mb_combine']
    threads: config['threads_combine']
    conda:
        "envs/star.yaml"
    log:
        "logs/pseudogenomes/{pseudo}_index.log"
    script:
        "scripts/star_gi.py"


rule pseudo_dict:
    input:
        "pseudogenomes/{pseudo}/{pseudo}_GCF_001858045.2.fa"
    output:
        "pseudogenomes/{pseudo}/{pseudo}_GCF_001858045.2.dict"
    resources:
        mem_mb=config['mem_mb_combine']
    threads: config['threads_combine']
    conda:
        "envs/picard.yaml"
    log:
        "logs/pseudogenomes/{pseudo}_dictionary.log"
    script:
        "scripts/create_genome_dictionary.py"



rule pseudo_fai:
    input:
        "pseudogenomes/{pseudo}/{pseudo}_GCF_001858045.2.fa"
    output:
        "pseudogenomes/{pseudo}/{pseudo}_GCF_001858045.2.fa.fai"
    conda:
        "envs/samtools.yaml"
    shell:
        "samtools faidx pseudogenomes/{pseudo}/{pseudo}_GCF_001858045.2.fa"
