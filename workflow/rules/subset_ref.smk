rule subset_ref_vcf:
    priority: 1
    input:
        vcf="calls/selected_ref_RNA.vcf",
    output:
        vcf="pseudogenomes/subset_vcf_file_wo_indels.recode.vcf"
    resources:
        mem_mb=200000
    threads: 20
    conda:
        "envs/vcftools.yaml"
    log:
        "logs/pseudogenomes/subset_wo_indels.log"
    script:
        "scripts/subset_ref_vcf.py"

rule bgzip:
    input:
        vcf="pseudogenomes/subset_vcf_file_wo_indels.recode.vcf",
    output:
        vcf="pseudogenomes/subset_vcf_file_wo_indels.recode.vcf.gz"
    conda:
        "envs/bgzip.yaml"
    log:
        "logs/pseudogenomes/bgzip.log"
    shell:
        # gzip and index the vcf file created from the reference genome
        "bgzip {input.vcf} "

rule backup:
    input:
        vcf="pseudogenomes/subset_vcf_file_wo_indels.recode.vcf.gz"
    output:
        vcf="pseudogenomes/subset_vcf_file_wo_indels.recode.vcf.gz.copy"
    conda:
        "envs/bgzip.yaml"
    log:
        "logs/pseudogenomes/backup.log"
    shell:
        "cp {input.vcf} {output.vcf}"

rule tabix:
    priority: 1
    input:
        vcf="pseudogenomes/subset_vcf_file_wo_indels.recode.vcf.gz"
    output:
        "pseudogenomes/subset_vcf_file_wo_indels.recode.vcf.gz.tbi"
    conda:
        "envs/bgzip.yaml"
    log:
        "logs/pseudogenomes/bgzip.log"
    shell:
        # gzip and index the vcf file created from the reference genome
        "tabix {input.vcf}"
