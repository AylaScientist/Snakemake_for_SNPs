# Convert th vcf file to annovar format thus extracting the SNPs
rule convert_to_annovar:
    input:
        gvcf1="calls/selected_{pseudo}.vcf"
    output:
        o1="annotated/annovar_Nile_{pseudo}"
    resources:
        mem_mb=200000
    threads: 20
    conda:
        "envs/perl.yaml"
    shell:
        "perl scripts/convert2annovar.pl {input.gvcf1} -format vcf4 -allsample -withfreq -withfilter -context -out {output.o1}"



#Make tokens for the annotation script
rule token_annotation:
    input:
        "annotated/annovar_Nile_{pseudo}"
    output:
        o1 = "annotated/annotated_Nile_snps_{pseudo}"
    shell:
        "touch {output.o1} "



# Collect the annotations in the db
rule annotate:
    input:
        i1="annotated/annovar_Nile_{pseudo}",
        i2="annotated/annotated_Nile_snps_{pseudo}",
        buildver="ON/ON"
    output:
        o2="annotated/annotated_Nile_snps_{pseudo}.variant_function"
    resources:
        mem_mb=200000
    threads: 20
    conda:
        "envs/perl.yaml"
    shell:
        "perl scripts/annotate_variation.pl -geneanno {input.i1} -buildver {input.buildver} ./ -outfile {input.i2}"



# Make tokens for the tables with the annotation:
rule token_table:
    input:
        "annotated/annotated_Nile_snps_{pseudo}.variant_function"
    output:
        "annotated/annotated_all_snps_{pseudo}"
    shell:
        "touch {output}"


# Annotate the file. This file can be used for the creation of the pseudogenomes
rule an_table:
    input:
        i1="annotated/annotated_Nile_snps_{pseudo}",
        i2="annotated/annotated_Nile_snps_{pseudo}.variant_function",
        i3="annotated/annotated_all_snps_{pseudo}",
        gvcf1="calls/all_{pseudo}.vcf",
        path="ON/", #Path to the database (buildver)
        buildver="ON"
    output:
        o1=temp("annotated/annotated_all_snps_{pseudo}.ON_multianno.vcf")
    resources:
        mem_mb=200000
    threads: 20
    conda:
        "envs/perl.yaml"
    shell:
        "perl scripts/table_annovar.pl {input.gvcf1} {input.path} -buildver {input.buildver} -out {input.i3} -remove -protocol refGene -operation g -nastring . -vcfinput"
