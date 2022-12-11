__author__ = "Aurora Campo"
__copyright__ = "Copyright 2022, Aurora Campo"
__email__ = "ayla.bcn@gmail.com"
__license__ = "MIT"


import os
import pandas as pd
import numpy as np
from os import path
from typing import List, Any, Generator
from snakemake.shell import shell
from snakemake_wrapper_utils.java import get_java_opts

shell.executable("bash")


"""
ELIMINATE MONOALLELIC EXPRESSION FROM THE ALLELEIC IMBALANCE
------------------------------------------------------------
The alleles that are expressed only in one of the tissues analysed
will be sent to a different path for studing the monoallelic expression
Then, the alleles expressed in both tissues can be studied separately
"""




def frequencies (df, samples):
    # Calculate the allele frequencies
    for sample in samples:
        sample_r_ad = str(sample + "_R_.AD")
        sample_a_ad = str(sample + "_A_.AD")
        af_sample = str("AF_" + sample)
        df[af_sample] = df[sample_r_ad]/(df[sample_r_ad] + df[sample_a_ad])
    df = df.fillna(0)
    print ("Allele frequencies calculated")

    return df




def MAE (df,samples):
    # Collect the indexes whose allele frequencies are 0 or 1 in both tissues of the same sample. These indexes correspond
    # to the SNPs that are MAE for at least one sample
    if str(path.exists("config/Samples_MAE.csv"))==True:
        mae = pd.DataFrame() # Create an empty dataframe
        tissues = pd.read_csv(samples_mae)
        tissues = pd.DataFrame(tissues)
        first_samples = tissues.iloc[:, 0]
        second_samples = tissues.iloc[:, 1]
        i = 0
        for first_sample in first_samples:
            print("MAE loop in sample ", first_sample, " and sample ", second_samples[i])
            af1 = str("AF_" + first_sample)
            af2 = str("AF_" + second_samples[i])
            index_af = df[((df[af1] == 0)  & (df[af2] == 0)) | (df[af1] == 1)  & (df[af2] == 1)].index
            df.drop(index=index_af, inplace=True)
            mae = pd.merge(mae, index_af, on=('CHROM', 'POS'))
            i = i + 1
    else:
        mae = pd.DataFrame()
        for sample in samples:
            af = str("AF_" + sample)
            index_af = df[(df[af] == 0) | (df[af] == 1)].index
            df.drop(index=index_af, inplace=True)
            mae = pd.merge(mae, index_af, on=('CHROM', 'POS'))

    return df, mae




def main():

    # Read the files from snakemake:
    extra = snakemake.params.get("extra")
    java_opts = snakemake.params.get("java_opts")

    log = snakemake.log_fmt_shell(stdout=False, stderr=True)


    PATH = os.getcwd()
    os.chdir("temp")
    uniform = snakemake.input.get("result1")
    df_uni = pd.read_csv(uniform, low_memory=False)
    df_uni = pd.DataFrame(df_uni)
    os.chdir(PATH)
    names = snakemake.input.get("names")
    sample_names = pd.read_csv(names)
    samples_mae = snakemake.input.get("mae")

    results_mae = snakemake.input.get("resut_mae")

    # Create arrays of the sample names and the pseudogenome codes
    samples = sample_names["Sample_name"].values

    # Mine the data
    df_af = frequencies(df_uni, samples)
    df_no_mae = MAE(df_af, samples)
    df_no_mae.to_csv()
    df_mae.to_csv()


if __name__ == '__main__':
    main()


shell(
    "{log}"
)
