# Convert th vcf file to annovar format thus extracting the SNPs
rule convert_to_annovar:
    input:
        gvcf1="calls/selected_{pseudo}.vcf"
    output:
        o1= config['params']['annotation']['convert']+'{pseudo}'
    resources:
        mem_mb=config['mem_mb']
    threads: config['threads']
    conda:
        "envs/perl.yaml"
    shell:
        "perl scripts/convert2annovar.pl {input.gvcf1} -format vcf4 -allsample -withfreq -withfilter -context -out {output.o1}"



#Make tokens for the annotation script
rule token_annotation:
    input:
        config['params']['annotation']['convert']+'{pseudo}'
    output:
        o1 = config['params']['annotation']['output_annotate']+'{pseudo}'
    shell:
        "touch {output.o1} "



# Collect the annotations in the db
rule annotate:
    input:
        i1=config['params']['annotation']['convert']+'{pseudo}',
        buildver= config['params']['annotation']['pathbuild'],
        i2=config['params']['annotation']['output_annotate']+'{pseudo}'
    output:
        o2=config['params']['annotation']['output_annotate']+"{pseudo}.variant_function"
    resources:
        mem_mb=config['mem_mb']
    threads: config['threads']
    conda:
        "envs/perl.yaml"
    shell:
        "perl scripts/annotate_variation.pl -geneanno {input.i1} -buildver {input.buildver} ./ -outfile {input.i2}"



# Make tokens for the tables with the annotation:
rule token_table:
    input:
        config['params']['annotation']['output_annotate']+"{pseudo}.variant_function"
    output:
        "annotated/annotated_all_snps_{pseudo}"
    shell:
        "touch {output}"


# Annotate the file. This file can be used for the creation of the pseudogenomes
rule an_table:
    input:
        i1=config['params']['annotation']['convert']+'{pseudo}',
        i2=config['params']['annotation']['output_annotate']+"{pseudo}.variant_function",
        i3=config['params']['annotation']['output_annotate']+'{pseudo}',
        gvcf1="calls/all_{pseudo}.vcf",
        path=config['params']['annotation']['path'], #Path to the database (buildver)
        buildver=config['params']['annotation']['buildver']
    output:
        o1="annotated/annotated_all_snps_{pseudo}"+config['params']['annotation']['output']
    resources:
        mem_mb=config['mem_mb']
    threads: config['threads']
    conda:
        "envs/perl.yaml"
    shell:
        "perl scripts/table_annovar.pl {input.gvcf1} {input.path} -buildver {input.buildver} -out {input.i3} -remove -protocol refGene -operation g -nastring . -vcfinput"
