rule gatk_select:
    input:
        vcf="calls/all_ref_g.vcf",
        ref=config['ref']['genome'],
    output:
        vcf="calls/selected_ref_RNA.vcf"
    conda:
        "envs/gatk.yaml"
    log:
        "logs/gatk/select/snvs_Nile.log"
    params:
        extra="--restrict-alleles-to BIALLELIC",  # optional filter arguments, see GATK docs
        java_opts=config['java_opts'],
    threads: config['threads']
    resources:
        mem_mb= config['mem_mb']
    script:
        "scripts/selectvariants.py"
