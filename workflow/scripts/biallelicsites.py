__author__ = "Aurora Campo"
__copyright__ = "Copyright 2021, Aurora Campo"
__email__ = "ayla.bcn@gmail.com"
__license__ = "MIT"


import os
import pandas as pd
import numpy as np
from os import path
from snakemake.shell import shell
from snakemake_wrapper_utils.java import get_java_opts
shell.executable("bash")

extra = snakemake.params.get("extra")
java_opts = snakemake.params.get("java_opts")
log = snakemake.log_fmt_shell(stdout=False, stderr=True)


"""
DEL MULTIALLELIC SITES
----------------------

Compare the genotypes for assessing the counts on each allele
Set the data for the mapping against the pseudogenome of sample PSG1 because it
has the variants of the reference genotype. The SNPs called by this sample are
set as reference and alternative alleles for all the other samples.

The tables containing the SNPs from mapping to PSG1 and PSG2 contained only
biallelic sites. However multiallelic possibilities can appear if one SNP was
called for a different genotpe in each of the pseudogenomes. The multiallelic
sites will be deleted from this dataframe to continue the analysis with
bialleleic sites only.

Delete multiallelic sites:
Let's start droping the multiallelic. The criteria compares if the reference
allele in the SNPs resulting from mapping against PSG1 pseudogenome is different
form the reference and the alternative alleles from the mapping against PSG2
pseudogenome. This is later repeated for the alternative allele. Note that the
SNPs mapped against PSG1 pseudogenome have been pre-filtered for biallelic sites
for all the samples.

Collect the indexes where the reference allele in PSG1 is different from both
alleles in the other mappings or the alternative allele in PSG1 is different
from both alleles in the PSG2 mapping.
"""



# Make an iterator for the collection of multiallelic sites:
def multiallelic_sample(df, sample, PSGs):
    # For obtaining the column names:
    r = "_R_"
    a = "_A_"
    gt = ".GT"
    PSG1_code: str = PSGs[0]
    PSG2_code: str = PSGs[1]
    rPSG1_GT = int(df.columns.get_loc(str(sample + r + PSG1_code + gt)))
    aPSG1_GT = int(df.columns.get_loc(str(sample + a + PSG1_code + gt)))
    rPSG2_GT = int(df.columns.get_loc(str(sample + r + PSG2_code + gt)))
    aPSG2_GT = int(df.columns.get_loc(str(sample + a + PSG2_code + gt)))
    print("Multiallelic loop in sample: ", sample)

    sample_index = []
    for i in range(len(df)):
        if (((df.iloc[i, rPSG1_GT]) != df.iloc[i, rPSG2_GT]) & (df.iloc[i, rPSG1_GT] != df.iloc[i, aPSG2_GT])) | (
                (df.iloc[i, aPSG1_GT] != df.iloc[i, rPSG2_GT]) & (df.iloc[i, aPSG1_GT] != df.iloc[i, aPSG2_GT])):
            sample_index.append(i)
    return sample_index


# Eliminate the multiallelic sites
def multiallelic(df, samples, PSGs, multi_index):
    for sample in samples:
        multiallelic_sample_index = multiallelic_sample(df, sample, PSGs)
        multi_index = multi_index + multiallelic_sample_index
    df_bi = df.drop(index=multi_index)
    return df_bi

def main():
    
    # Add files
    mdf = pd.read_csv(snakemake.input.get("csv"), low_memory=False)
    mdf = pd.DataFrame(mdf)
    sample_names = pd.read_csv(snakemake.input.get("sn1"))
    PSG_codes = pd.read_csv(snakemake.input.get("psc"))

    # Create arrays of the sample names and the pseudogenome codes
    samples = sample_names['Sample_name'].values
    PSGs = PSG_codes['PSGs'].values

    # Create an empty array for collect the indexes with multiallelic sites:
    multi_index = []

    # Mine the datafiles
    df_bi = multiallelic(mdf, samples, PSGs, multi_index)
    df_bi.to_csv(snakemake.output[0])


if __name__ == '__main__':
    main()


shell(
    "{log}"
)
