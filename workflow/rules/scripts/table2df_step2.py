# python 3.7

"""
python 3.7

    @autor : Ayla
    @version : 0.1
    @date : November 2020

This script prepares the tables that come from the server after the call of SNPs to fit the format of the dataframe.

The format of the table is the same for both dataframes. Columns from 1 to 6 are:
    CHROM: Chromosome,
    POS: position,
    Gene.refGene: Gene_ID
    Func.refGene: Annotation of the function for this SNP
    ExonicFunc.refGene: Annotation of the function if this SNP is exonic
    AF: Allele frequency estimated for all the samples

From column 7th till the end there will be 4 columns for each sample. The structure of the column name is
as follows:
    NAME (determined by the sample name)
    _R or _A (reference or alternative allele)
    _??? (The pseudogenome code)
    .AD (allele depth)
    .GT (Genotype)
"""
import pandas as pd
import numpy as np
import os
from os import path
from snakemake import shell


log = snakemake.log_fmt_shell(stdout=False, stderr=True)

def main():
    """
    Input the files and convert them into data frames
    """

    # Read the dataframe and convert it into pandas dataframe. tb stands for table:

    tb1_colnames = pd.read_csv ( snakemake.input.get("tb_colnames") )
    tb1_cols = tb1_colnames['Col_names'].values
    print("Cols read")

    #tb1 = pd.read_csv(snakemake.input.get("table"), sep = "\,|/|\t", header=None, engine = "python")
    tb1 = pd.read_csv(snakemake.input.get("table"), sep = "\,|/|\t", names = tb1_cols, engine = "python")
    print("Table read")

    tb1 = pd.DataFrame ( tb1 )

    tb1.drop ( tb1.index[:1], inplace=True )

    tb1.to_csv ( snakemake.output.get("csv") )
    print("Table with correct column names is created for each pseudogenome in the folder called \"variants\".")
    print("Please verify the order of the samples on the new .csv table is the same as in the header of the .table file")

if __name__ == '__main__':
    main()

shell(
    "{log}"
)
