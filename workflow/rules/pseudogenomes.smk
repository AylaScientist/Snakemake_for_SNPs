
rule pseudogenome:
    input:
        vcf="pseudogenomes/subset_vcf_file_wo_indels.recode.vcf.gz",
        ref=config['ref']['genome'],
        index="pseudogenomes/subset_vcf_file_wo_indels.recode.vcf.gz.tbi"
    output:
        pseudo = "pseudogenomes/{pseudo}/{pseudo}_"+config['ref']['release']+".fa"
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
        "pseudogenomes/{pseudo}/{pseudo}_"+config['ref']['release']+".fa"
    output:
        "pseudogenomes/{pseudo}/SAindex"
    params:
        threads= "16",
        dir = "pseudogenomes/{pseudo}/",
        annotation = config['ref']['gtf'],
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
        "pseudogenomes/{pseudo}/{pseudo}_"+config['ref']['release']+".fa"
    output:
        "pseudogenomes/{pseudo}/{pseudo}_"+config['ref']['release']+".dict"
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
        "pseudogenomes/{pseudo}/{pseudo}_"+config['ref']['release']+".fa"
    output:
        "pseudogenomes/{pseudo}/{pseudo}_"+config['ref']['release']+".fa.fai"
    conda:
        "envs/samtools.yaml"
    shell:
        "samtools faidx {input}"
