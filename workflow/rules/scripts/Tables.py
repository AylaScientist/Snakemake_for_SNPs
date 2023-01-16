# python 3.7
# Tables.py


"""
python 3.7

    @autor : Ayla
    @version : 0.1
    @date : August 2020

In this script we will process one table or dataframe containing SNPs that have been processed with the script
"ASE_workflow.py" and whose quality has been checked for mapping bias with the script "QC.py". These data was treated by
the script "Stats.py" in order to determine the significance of each SNP for treatment (CHI² test), for tissue (Fisher
exact test) or for sample (Binomial test).

The annotation of these SNPs is been performed with ANNOVAR and the present dataframe took the names of the columns
according to the annotation file (.gff3 or gft) provided for the reference genome of the species. Therefore, the names
of the columns such as the kind of SNP and the function may vary from genome to genome and from species to species.

These SNPs to be worked here are unbiased and biallelic, and the total monoallelic expression (MAE) in the same
individual has been discarded.
This script produces the next tables as output:

    1 Treatment table: one column for the SNP_IDs from each treatment (significant for each CHI²). The number of columns
    will vary with each experimental setup
        This table will be used for a Venn diagram

    2 Types of polymorphism by test: the columns will be the type of polymorphism, the GO terms and the relative values of
    these annotations. One row correspond to a different Chi² tests.

    3 Tables for DEG and ASE SNPs: a table that will correlate the significant DEG (differentially expressed genes)
    with the ASE SNPs for the same test, the annotation of these genes and the function of the SNP.

    4 Table for SNPs in DEG and SNPs in ASE: One table with two columns for the SNP_ID corresponding to total SNPs in ASE
    and another column with the SNPs included in the significant DEGs.
        This table will be used for a Venn diagram

    5 Compare the % of intronic and exonic SNPs (from the total SNPs) in Scaffolds and in Chromosomes. The columns will be
    intronic and exonic and there will be two rows: Chromosomes and Scaffolds.
        This table is related to the Heatmap (also drawn) and to the Manhattan plot

    6 Physiological function: This table will include a column for each treatment (CHI²) and a row for each GO function.
    The cell content will be the % of SNPs in each GO function for treatment.

    7 Results overview: This table includes the number of SNPs and genes in ASE for a particular test

    8 Results: Out of these SNPs in ASE for a particular test, the SNPs that are in AI for a level of a factor.
    Comes from the binomial tests.
    (table with these SNPs location – exon, 5’, 3’…)


NOTE: The next dictionaries must be placed on the same folder:
    dict = pd.read_csv('GeneID_to_Genebank.csv')
    mega_dict = pd.read_excel ( 'Tilapia gene transcript protein GO KEGG uniprot.xlsx' )
"""

# Import the libraries
import pandas as pd
import numpy as np
import os
from os import path
from matplotlib import pyplot as plt
from math import pi
from scipy.stats import uniform
from scipy.stats import randint
import seaborn as sns
from snakemake.shell import shell
from snakemake_wrapper_utils.java import get_java_opts
shell.executable("bash")

extra = snakemake.params.get("extra")
java_opts = snakemake.params.get("java_opts")
log = snakemake.log_fmt_shell(stdout=False, stderr=True)


def GO_graph (df_stats,experiments):
    # This function plots the GO terms that show SNPs in ASE. It also produces a table with the genes and SNPs in ASE
    #tests=["test1","test2","test3","test4"]
    tests=experiments
    for test in tests :
        test_name = str(test) + "_CHI_p-fdr"
        if test_name in df_stats:
            # Table 7
            index_exp = df_stats[(df_stats[test_name] >= 0.05)].index
            index_nan = df_stats[df_stats[test_name].isnull ()].index
            df_sig = df_stats.drop ( index=index_exp )  # Collect the significant SNPs in a table for this experiment
            df_exp = df_sig.drop ( index=index_nan )  # Drop the empty values for the tests
            genes_sig = df_exp["Gene_ID"].unique ()
            results = { 'SNPs in ASE':  [len(df_exp)],
                        'Genes containing ASE': [len(genes_sig)]}
            df_results = pd.DataFrame (results, columns = ['SNPs in ASE','Genes containing ASE'])
            PATH = os.getcwd()
            os.chdir("results")
            df_results.to_csv(test_name + "_basic_results_chi_square.csv")
            df_exp.to_csv(test_name + "_ASE_description_chi_square.csv")
            os.chdir(PATH)

            # GO terms analysis
            mega_dict = pd.read_csv ( snakemake.input.get("i1"), low_memory=False )
            df_func = pd.merge ( df_stats, mega_dict, on='Gene_ID', how='right' )
            df_func.to_csv ( snakemake.output.get("o1") )

            GO_dict = pd.read_csv ( snakemake.input.get("i2"), low_memory=False )

            df_GO = pd.merge(df_exp,GO_dict, on = "Genebank", how = "left")
            GO_counts = df_GO['GO term accession'].value_counts ()
            GO_counts = pd.DataFrame(GO_counts)
            GO_counts.columns = ['Frequency']
            GO_most = GO_counts.head(15)

            # Pie chart
            sizes = []
            labels = GO_most.index.values
            sizes.append (GO_most['Frequency'].tolist() )
            fig1, ax1 = plt.subplots ()
            ax1.pie ( sizes[0], labels=labels, autopct='%1.1f%%', shadow=True, startangle=90 )
            ax1.axis ( 'equal' )  # Equal aspect ratio ensures that pie is drawn as a circle.
            PATH = os.getcwd ()
            os.chdir("results")
            filename = test_name + "_GO_pie_chart_chi_square.svg"
            plt.savefig ( filename )
            plt.clf()
            os.chdir(PATH)

        else:
            test_name = str(test) + "_Fisher_p-fdr"
                 # Table 7
            index_exp = df_stats[(df_stats[test_name] >= 0.05)].index
            index_nan = df_stats[df_stats[test_name].isnull ()].index
            df_sig = df_stats.drop ( index=index_exp )  # Collect the significant SNPs in a table for this experiment
            df_exp = df_sig.drop ( index=index_nan )  # Drop the empty values for the tests
            genes_sig = df_exp["Gene_ID"].unique ()
            results = { 'SNPs in ASE':  [len(df_exp)],
                        'Genes containing ASE': [len(genes_sig)]}
            df_results = pd.DataFrame (results, columns = ['SNPs in ASE','Genes containing ASE'])
            PATH = os.getcwd()
            os.chdir("results")
            df_results.to_csv(test_name + "_basic_results_fisher_test.csv")
            df_exp.to_csv(test_name + "_ASE_description_fisher_test.csv")
            os.chdir(PATH)

            # GO terms analysis
            mega_dict = pd.read_csv ( snakemake.input.get("i1"), low_memory=False )
            df_func = pd.merge ( df_stats, mega_dict, on='Gene_ID', how='right' )
            df_func.to_csv ( snakemake.output.get("o1") )

            GO_dict = pd.read_csv ( snakemake.input.get("i2"), low_memory=False )

            df_GO = pd.merge(df_exp,GO_dict, on = "Genebank", how = "left")
            GO_counts = df_GO['GO term accession'].value_counts ()
            GO_counts = pd.DataFrame(GO_counts)
            GO_counts.columns = ['Frequency']
            GO_most = GO_counts.head(15)

            # Pie chart
            sizes = []
            labels = GO_most.index.values
            sizes.append (GO_most['Frequency'].tolist() )
            fig1, ax1 = plt.subplots ()
            ax1.pie ( sizes[0], labels=labels, autopct='%1.1f%%', shadow=True, startangle=90 )
            ax1.axis ( 'equal' )  # Equal aspect ratio ensures that pie is drawn as a circle.
            PATH = os.getcwd ()
            os.chdir("results")
            filename = test_name + "_GO_pie_chart_fisher_test.svg"
            plt.savefig ( filename )
            plt.clf()
            os.chdir(PATH)       

        

def tables_1_2(df_stats, experiments):
    # The columns will be the type of polymorphism, the GO terms and the relative values of these annotations. One row
    # correspond to a different test.
    # Initialize the tables and the columns that will be in table 1 and table 2:
    all_treatment_snp = pd.DataFrame ()
    df_polymorphism = pd.DataFrame ()
    df_stats = df_stats.fillna ( 1 )  # to drop rows with empty values

    for experiment in experiments:
        # Initalize the intermediate table:
        col_name = str ( experiment ) + "_CHI_p-fdr"
        if col_name in df_stats:
            index_exp = df_stats[(df_stats[col_name] >= 0.05)].index
            df_exp = df_stats.drop ( index=index_exp )  # Collect the significant SNPs in a table for this experiment
            df_new = pd.DataFrame ()
            df_new[col_name] = df_exp["SNP_ID"]

            # Table 1
            all_treatment_snp = pd.concat ( [all_treatment_snp, df_new], axis=1 )

            # Table 2
            func_refgene = df_exp['Func.refGene_x'].value_counts ()
            exonicfunc_refgene = df_exp['ExonicFunc.refGene_x'].value_counts ()
            polymorphism_exp = pd.Series ( func_refgene.append ( exonicfunc_refgene ) )
            df_polymorphism = pd.concat ( [df_polymorphism, polymorphism_exp.to_frame ().T] )
            df_polymorphism.rename ( index={0: col_name}, inplace=True )  # add row name that is the column of the origin
        else:
            col_name = str ( experiment ) + "_Fisher_p-fdr"
            index_exp = df_stats[(df_stats[col_name] >= 0.05)].index
            df_exp = df_stats.drop ( index=index_exp )  # Collect the significant SNPs in a table for this experiment
            df_new = pd.DataFrame ()
            df_new[col_name] = df_exp["SNP_ID"]

            # Table 1
            all_treatment_snp = pd.concat ( [all_treatment_snp, df_new], axis=1 )

            # Table 2
            func_refgene = df_exp['Func.refGene_x'].value_counts ()
            exonicfunc_refgene = df_exp['ExonicFunc.refGene_x'].value_counts ()
            polymorphism_exp = pd.Series ( func_refgene.append ( exonicfunc_refgene ) )
            df_polymorphism = pd.concat ( [df_polymorphism, polymorphism_exp.to_frame ().T] )
            df_polymorphism.rename ( index={0: col_name}, inplace=True )  # add row name that is the column of the origin

     
    # Table 1:
    all_treatment_snp.to_csv ( snakemake.output.get("o3") )
    
    # Table 2:
    df_polymorphism = df_polymorphism.fillna ( 0 )
    df_polymorphism.to_csv ( snakemake.output.get("o4") )
   
    # Radar chart for all treatents :
    # subset SNP Function and Exonic SNP Function
    categories_exo = list ( df_polymorphism )[11:20]  # It avoids the category "." that is non-informative
    print("Categories exo: ", categories_exo)
    df_pol_fun = df_polymorphism[2:9]
    df_pol_exo = df_polymorphism[categories_exo]
    

    # Polymorphism function by treatment:
    categories = list ( df_pol_fun )[0:]  # For function of SNP 1:9. Exonic function is from 10:20
    print(categories)
    values_list = []
    angles_list = []

    treatments = df_pol_fun.index.values
    fig, ax = plt.subplots ( nrows=1, ncols=1, figsize=(8, 8), subplot_kw=dict ( polar=True ) )
    for i in range ( len ( df_pol_fun ) ):
        values = df_pol_fun[i:(i + 1)].values.flatten ().tolist ()
        values += values[:1]  # repeat the first value to close the circular graph
        values_list.append ( values )
        angles = [n / float ( len ( categories ) ) * 2 * pi for n in range ( len ( categories ) )]
        angles += angles[:1]
        angles_list.append ( angles )
        ax.plot ( angles_list[i], values_list[i], linewidth=1, linestyle='solid', label=treatments[i] )
        ax.fill ( angles_list[i], values_list[i], alpha=0.4 )
    plt.xticks ( angles[:-1], categories, color='grey', size=12 )
    plt.yticks ( np.arange ( 1, 6 ), ['1', '2', '3', '4', '5'], color='grey', size=12 )
    plt.ylim ( 0, 800 )
    ax.set_rlabel_position ( 30 )
    plt.legend ( loc='upper right', bbox_to_anchor=(0.1, 0.1) )
    #plt.show ()
    PATH = os.getcwd ()
    os.chdir("results")
    plt.savefig("Polymorphism_function_by_treatment_radar_chart.svg")
    os.chdir(PATH)

    plt.clf ()
    

    # Polymorphism exonic function:
    categories = list ( df_pol_exo )[0:]  # For function of SNP 1:9. Exonic function is from 10:20
    values_list = []
    angles_list = []
    treatments = df_pol_exo.index.values
    fig, ax = plt.subplots ( nrows=1, ncols=1, figsize=(8, 8), subplot_kw=dict ( polar=True ) )
    for i in range ( len ( df_pol_exo ) ):
        values = df_pol_exo[i:(i + 1)].values.flatten ().tolist ()
        values += values[:1]  # repeat the first value to close the circular graph
        values_list.append ( values )
        angles = [n / float ( len ( categories ) ) * 2 * pi for n in range ( len ( categories ) )]
        angles += angles[:1]
        angles_list.append ( angles )
        ax.plot ( angles_list[i], values_list[i], linewidth=1, linestyle='solid', label=treatments[i] )
        ax.fill ( angles_list[i], values_list[i], alpha=0.4 )
    plt.xticks ( angles[:-1], categories, color='grey', size=12 )
    plt.yticks ( np.arange ( 1, 6 ), ['1', '2', '3', '4', '5'], color='grey', size=12 )
    plt.ylim ( 0, 400 )
    ax.set_rlabel_position ( 30 )
    plt.legend ( loc='upper right', bbox_to_anchor=(0.1, 0.1) )
    PATH = os.getcwd ()
    os.chdir("results")
    plt.savefig("Polymorphism_exonic_function_by_treatment.svg")
    os.chdir(PATH)

    plt.clf ()

    # Pie charts :
    labels = df_polymorphism.columns
    styles = list ( plt.style.available )
    for i in range ( (df_polymorphism.shape[0]) ):
        sizes = []
        # Using iloc to access the values of
        # the current row denoted by "i"
        plt.style.use ( styles[5] )
        sizes.append ( list ( df_polymorphism.iloc[i, :] ) )
        fig1, ax1 = plt.subplots ()
        ax1.pie ( sizes[0], labels=labels, autopct='%1.1f%%', shadow=True, startangle=90 )
        ax1.axis ( 'equal' )  # Equal aspect ratio ensures that pie is drawn as a circle.
        #plt.show ()
        PATH = os.getcwd ()
        os.chdir("results")
        filename = df_polymorphism.index.values[i] + "_pie_chart.svg"
        plt.savefig(filename)
        os.chdir(PATH)

        plt.clf ()


def heatmap(df, group_names):
    # This heatmap illustrates the allele frequencies for the genomic positions for identifying regions in AI
    af_groups = []
    for group_name in group_names:
        af_group = "AF_" + group_name
        af_groups.append ( af_group )

    htmap = df[af_groups]
    sns.heatmap ( htmap, cmap='YlGnBu' )
    #plt.show ()
    PATH = os.getcwd ()
    os.chdir("results")
    plt.savefig("Heatmap.svg")
    os.chdir(PATH)

    plt.clf ()


def table5(df_stats):
    # This part of the analysis have to be adapted to the different strings on each genome identifying chromosomes and scaffolds
    index_scf = df_stats[df_stats['CHROM'].str.match ( r'^NC_' )].index # Substitute the chain NC_ by the pattern on the genome of the organism you work with
    df_scf = df_stats.drop ( index_scf )
    intronic_scf = df_scf["Func.refGene_x"].value_counts ()
    exonic_scf = df_scf["ExonicFunc.refGene_x"].value_counts ()

    index_chrom = df_stats[df_stats['CHROM'].str.match ( r'^NW_' )].index # Substitute the chain NW_ by the pattern on the genome of the organism you work with
    df_chrom = df_stats.drop ( index_chrom )
    intronic_chrom = df_chrom["Func.refGene_x"].value_counts ()
    exonic_chrom = df_chrom["ExonicFunc.refGene_x"].value_counts ()
    print("Intronic chrom: ",intronic_chrom)
    print("Exonic chrom: ", exonic_chrom)


    # Create pandas DataFrame.
    data = {
        'Intronic': [sum ( intronic_chrom ) - intronic_chrom["exonic"], sum ( intronic_scf ) - intronic_scf["exonic"]],
        'Exonic': [sum ( exonic_chrom ) - exonic_chrom["."], sum ( exonic_scf ) - exonic_scf["."]]}
    table5 = pd.DataFrame ( data, index=["Chromosome", "Scaffold"] )
    table5["Intronic %"] = table5["Intronic"] * 100 / (table5["Intronic"] + table5["Exonic"])
    table5["Exonic %"] = table5["Exonic"] * 100 / (table5["Intronic"] + table5["Exonic"])
    table5.to_csv ( snakemake.output.get("o5") )


def Manhattan_plot(df_stats, experiments):
    # This function builds a Manhattan plot
    #plt.rcParams.update ( {'font.size': 10} )
    for experiment in experiments:
        index_chrom = df_stats[df_stats['CHROM'].str.match ( r'^NW_' )].index
        df_chrom = df_stats.drop ( index_chrom )

        df = pd.DataFrame ( {"SNP": df_chrom["SNP_ID"],
                             "pvalue": df_chrom[str ( experiment ) + "_CHI_p-val"],
                             "Chromosome": df_chrom["CHROM"]} )
        # -log_10(pvalue)
        df['minuslog10pvalue'] = -np.log10 ( df.pvalue )
        df.Chromosome = df.Chromosome.astype ( 'category' )
        df = df.sort_values ( 'Chromosome' )

        # How to plot gene vs. -log10(pvalue) and colour it by chromosome?
        df['ind'] = range ( len ( df ) )
        df_grouped = df.groupby ( 'Chromosome' )

        fig = plt.figure ()
        ax = fig.add_subplot ( 111 )
        colors = (
        "#a7414a", "#696464", "#00743f", "#563838", "#6a8a82", "#a37c27", "#5edfff", "#282726", "#c0334d", "#c9753d")
        # colors = ['red', 'green', 'blue', 'yellow']
        x_labels = []
        x_labels_pos = []
        for num, (name, group) in enumerate ( df_grouped ):
            group.plot ( kind='scatter', x='ind', y='minuslog10pvalue', color=colors[num % len ( colors )], ax=ax )
            x_labels.append ( name )
            x_labels_pos.append ( (group['ind'].iloc[-1] - (group['ind'].iloc[-1] - group['ind'].iloc[0]) / 2) )
        ax.set_xticks ( x_labels_pos )
        ax.set_xticklabels ( x_labels )
        ax.set_xlim ( [0, len ( df )] )
        ax.set_ylim ( [0, 2] )
        ax.set_xlabel ( 'Chromosome' )
        #plt.show ()
        PATH = os.getcwd ()
        os.chdir("results")
        plt.savefig(experiment + "_Manhattan_plot.svg")
        os.chdir(PATH)

        plt.clf ()



def main():
    """
    Read the files
    """

    df_stats = pd.DataFrame ( pd.read_csv ( snakemake.input.get("i3") ) )
    sample_names = pd.DataFrame ( pd.read_csv ( snakemake.input.get("i4")) )
    groups_df = pd.DataFrame ( pd.read_csv ( snakemake.input.get("i5") ) )
    exp = pd.DataFrame ( pd.read_csv ( snakemake.input.get("i6") ) )

    # Rename the key Gene ID and add the Genebank accessions for further analysis
    df_stats.rename ( columns={"Gene.refGene_x": "Gene_ID"}, inplace=True )
    # Read the dictionaries and merge with the df:
    dict = pd.read_csv ( snakemake.input.get("i7") )
    df_stats = df_stats.merge ( dict, on='Gene_ID', how='left' )


    # Create a numpy array of arrays with the samples of each group and the total samples of the experiment
    groups = groups_df.values

    # Create a numpy array with the name of each group
    group_names = list ( groups_df.columns[0] )

    # Create arrays of the sample names
    samples = sample_names["Sample_name"].values

    # Create an array with the name of each experimental test
    experiments = list ( exp["Test_ID"] )



    """
    SNP_ID and dictionary
    """
    # Create a column for SNP_ID
    df_stats["SNP_ID"] = df_stats.index
    SNP_dictionary = df_stats[['SNP_ID', 'CHROM', 'POS', 'Gene_ID', 'Func.refGene_x', 'ExonicFunc.refGene_x']]
    SNP_dictionary.to_csv ( snakemake.output.get("o2") )

    """
        Physiological function: This table will include a column for each treatment (CHI²) and a row for each GO function.
        The cell content will be the % of SNPs in each GO function for treatment.
        Table 7
        """
    GO_graph ( df_stats, experiments )


    """
    1 Treatment table: one column for the SNP_IDs from each treatment (significant for each CHI²). The number of columns
        will vary with each experimental setup
            This table will be used for a Venn diagram

    2 Types of polymorphism by test: the columns will be the type of polymorphism, the GO terms and the relative values of
    these annotations. One row correspond to a different Chi² test.
    """

    tables_1_2 ( df_stats, experiments )


    """
    5 Compare the % of intronic and exonic SNPs (from the total SNPs) in Scaffolds and in Chromosomes. The columns will
    be intronic and exonic and there will be two rows: Chromosomes and Scaffolds.

    PLOT HEATMAP
    ------------
    Use the csv file produced in the control quality for plotting the allele frequencies of each experimental group.
    The matrix has to be specific for allele frequencies per group and Chromosome position.
    Now make a general comparison between treatments. Most of the colors in these graphs are or 0 or 1, since they
    correspond to the allelic imbalance with significance.
    """

    heatmap ( df_stats, group_names )
    table5 ( df_stats )

    """
    MANHATTAN PLOT
    --------------
    A Mannhattan plot for each test for ASE SNPs will be performed.
    """
    Manhattan_plot ( df_stats, experiments )




if __name__ == '__main__':
    main ()
