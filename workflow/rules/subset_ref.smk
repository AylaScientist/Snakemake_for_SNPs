rule subset_ref_vcf:
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


rule tabix:
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