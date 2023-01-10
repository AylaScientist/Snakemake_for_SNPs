# python 3.7
# QC.py

__author__ = "Aurora Campo"
__copyright__ = "Copyright 2021, Aurora Campo"
__email__ = "ayla.bcn@gmail.com"
__license__ = "MIT"


"""
python 3.7

    @autor : Ayla
    @version : 1.0
    @date : September 2021

This script is a quality control of the data filtered by the script "ASE_workflow.py".
These SNPs are unbiased and biallelic, and the monoallelic expression (MAE) has been discarded.
This script will plot the distribution of the data in general and for each experimental group.
It needs a file with the name of the experimental group as the acronyms used in the sample name. For example:
    The sample number 1 corresponds to the gills exposed to freshwater
    The two experimental factors are "gills" (G) and "freshwater" (F)
    Then the sample is called "1GF" and the group name will be called "GF"
This script also needs a csv file with the first column naming the experimental groups and the next columns with the
sample names of each sample in a group
"""

import pandas as pd
import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import StrMethodFormatter
from snakemake.shell import shell
from snakemake_wrapper_utils.java import get_java_opts
shell.executable("bash")

extra = snakemake.params.get("extra")
java_opts = snakemake.params.get("java_opts")
log = snakemake.log_fmt_shell(stdout=False, stderr=True)


"""
    ALLELE FREQUENCIES
    ------------------
    Compute allele frequencies for each allele and later for each experimental group.
    """

def frequencies(df, samples):
    """
    :param df: dataframe with unbiased SNP counts without allele frequencies
    :param samples: list of strings including the sample names
    :return: dataframe with the allele frequencies by sample
    """
    # Calculate the allele frequencies
    for sample in samples:
        sample_r_ad = str(sample + "_R_.AD")
        sample_a_ad = str(sample + "_A_.AD")
        af_sample = str("AF_" + sample)
        df[af_sample] = df[sample_r_ad] / (df[sample_r_ad] + df[sample_a_ad])
    df = df.fillna(0)  # To fill cases where there are no counts and therefore AF is divided by 0
    print("Allele frequencies by sample calculated")

    return df

def allele_freqs(group_names, df):
    # Calculate the allele frequency for each experimental group in order to make a plot of each group:

    for group_name in group_names:
        af_group = "AF_" + group_name
        x_name = af_group + "_x"
        y_name = af_group + "_y"
        print(x_name)
        print(y_name)
        x = df.loc[:, df.columns.str.contains("([0-9])+" + group_name + "_R_.AD", regex=True)]
        y = df.loc[:, df.columns.str.contains("([0-9])+" + group_name + "_A_.AD", regex=True)]
        df[x_name] = x.sum(axis=1)
        df[y_name] = y.sum(axis=1)
        df[af_group] = df[x_name]/(df[x_name] + df[y_name])
        #df.drop(columns=x_name, inplace=True)
        #df.drop(columns=y_name, inplace=True)

    df = df.fillna(0)
    print("Allele frequencies by experimental group calculated")
    return df


def plot(df, group_names):
    # Plot all the histograms with the allele frequencies

    for group_name in group_names:
        af_group = "AF_" + group_name

        # Possibility 1
        # df.hist(column=af_group, bins=100, grid=False, figsize=(8,10), layout=(1,1), sharex=True, color='#86bf91', zorder=2, rwidth=0.9)

       # Possibility 3 previous in the script
        hist_data = df[af_group].values
        plt.hist(hist_data, 100, density=False, facecolor='#86bf91')
        plt.xlabel('Allele frequency')
        plt.ylabel('Frequency')
        plt.title(af_group)
        plt.grid(True)
        axes = plt.gca()
        axes.set_ylim([0, 15000])
        # plt.show()

        # Get the path
        PATH = os.getcwd ()
        os.chdir("results")
        plt.savefig(group_name + "_AF_mapping_bias.svg")
        os.chdir(PATH)

        plt.clf ()


def main():
    # Import the files
    df = pd.read_csv(snakemake.input.get("i1"))
    sample_names = pd.read_csv(snakemake.input.get("sn1"))
    # Create arrays of the sample names
    samples = sample_names['Sample_name'].values
    groups_df = pd.read_csv(snakemake.input.get("i2"))
    # Create a numpy array with the name of each group
    group_names = list(groups_df.columns[0])

    # Mine the data
    df_af = frequencies(df, samples)
    df_qc = allele_freqs(group_names, df_af)

    # Copy the quality control file for further calculation of the stats
    df_qc.to_csv(snakemake.output.get("results"))

    # Plot the histograms
    plot(df_qc, group_names)


if __name__ == '__main__':
    main()


shell(
    "{log}"
)
