rule results:
    input:
        i1=config['species_dict'],
        i2=config['DO_dict'],
        i3="results/SNPs_analysed_valid_genome.csv",
        i4=config['Sample_names'],
        i5=config['Experimental_groups'],
        i6=config['Experimental_design'],
        i7=config['GeneID_to_Genebank']
    output:
        o1=config['Results']['o1'],
        o2=config['Results']['o2'],
        o3=config['Results']['o3'],
        o4=config['Results']['o4'],
        o5=config['Results']['o5'],
    params:
        extra="",
        java_opts=config['java_opts'],
    threads: config['threads']
    resources:
        mem_mb=config['mem_mb']
    conda:
        "envs/python.yaml"
    log:
        "logs/python/results/tables.log"
    script:
        "scripts/Tables.py"
