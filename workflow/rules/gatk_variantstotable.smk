rule gatk_variantstotable:
    input:
        vcf="calls/selected_ref_RNA.vcf",
        ref=config['ref']['genome'],
    output:
        vcf="variants/AD_GT_counts_bi.table"
    conda:
        "envs/gatk.yaml"
    log:
        "logs/gatk/var2table/snvs_Nile.log"
    params:
        o1="annotated/annotated_all_snps.ON_multianno.vcf",
        extra="-SMA TRUE -F CHROM -F POS -F Gene.refGene -F Func.refGene -F ExonicFunc.refGene -F AF -GF AD -GF GT",  # optional filter arguments, see GATK docs
        java_opts=config['java_opts'],
    threads: config['threads']
    resources:
        mem_mb=config['mem_mb']
    script:
        "scripts/variantstotable.py"
        
