rule results_mae:
    input:
        i1="config/Tilapia gene transcript protein GO KEGG uniprot.csv",
        i2="config/GO_dictionary.csv",
        i3="results/results/QC_mae_valid_genome.csv",
        i4="config/Sample_names.csv",
        i5="config/Experimental_groups.csv",
        i6="config/Experimental_design.csv",
        i7="config/GeneID_to_Genebank.csv"
    output:
        o1="results/SNPs_function_mae.csv",
        o2="results/SNP_dictionary_mae.csv",
        o3="results/Treatment_SNPs_sig_all_for_Venn_mae.csv",
        o4="results/Summary_of_polymorphisms_mae.csv",
        o5="results/Intronic_SNPs_caryotype_mae.csv"
    log:
        "logs/python/results/tables.log"
    conda:
        "envs/python.yaml"
    resources:
        mem_mb=200000
    threads:20
    script:
        "scripts/Tables.py"
