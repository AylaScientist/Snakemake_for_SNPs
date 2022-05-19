rule gatk_variantstotable_PSG:
    input:
        vcf="annotated/annotated_all_snps_{pseudo}.ON_multianno.vcf",
        ref="pseudogenomes/{pseudo}/{pseudo}_GCF_001858045.2.fa"
    output:
        vcf="variants/AD_GT_counts_bi_{pseudo}.table"
    conda:
        "envs/gatk.yaml"
    log:
        "logs/gatk/var2table/snvs_Nile_{pseudo}.log"
    params:
        o1="annotated/annotated_all_snps_{pseudo}.ON_multianno.vcf",
        extra="-SMA TRUE -F CHROM -F POS -F Gene.refGene -F Func.refGene -F ExonicFunc.refGene -F AF -GF AD -GF GT",  # optional filter arguments, see GATK docs
        java_opts="-XX:MinRAMPercentage=80.0 -Xms200G -XX:-UseConcMarkSweepGC -XX:ParallelGCThreads=20 -XX:+UseTLAB",
    threads: 20
    resources:
        mem_mb=200000
    script:
        "scripts/variantstotable.py"
