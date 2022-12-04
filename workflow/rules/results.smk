rule results:
    input:
        i1="config/Tilapia gene transcript protein GO KEGG uniprot.csv",
        i2="config/GO_dictionary.csv",
        i3="results/SNPs_analysed_valid_genome.csv",
        i4="config/Sample_names.csv",
        i5="config/Experimental_groups.csv",
        i6="config/Experimental_design.csv",
        i7="config/GeneID_to_Genebank.csv"
    output:
        o1="results/SNPs_function.csv",
        o2="results/SNP_dictionary.csv",
        o3="results/Treatment_SNPs_sig_all_for_Venn.csv",
        o4="results/Summary_of_polymorphisms.csv",
        o5="results/Intronic_SNPs_caryotype.csv"
    log:
        "logs/python/results/tables.log"
    conda:
        "envs/python.yaml"
    resources:
        mem_mb=200000
    threads:20
    script:
        "scripts/Tables.py"
