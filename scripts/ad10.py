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


"""
DEL SNPs AD<10
---------------
This new function cleans the SNPs that do not have enough counts and are considered possible poor quality reads.
The SNPs with AD < 3 in one of the alleles are also cleaned.
"""


extra = snakemake.params.get("extra")
java_opts = snakemake.params.get("java_opts")
log = snakemake.log_fmt_shell(stdout=False, stderr=True)

av_df=snakemake.input.get("csv")
ad10_df=snakemake.output[0]

df_average = pd.read_csv(av_df, low_memory=False)
df_average = pd.DataFrame(df_average)
sample_names = pd.read_csv(snakemake.input.get("sn1"))
# Create arrays of the sample names
samples = sample_names['Sample_name'].values


def AD10(df, samples):
    # The iterator collects the index were AD from both alleles is below 10 or the AD for one allele is below 3
    # Keep all this analysis for the same individual including two tissues at a time
    i = 0
    indexes_ad10 = 0
    # There is only one tissue
    for sample in samples:
        print("AD10 filter in sample ", sample)
        sample_r_ad = str(sample + "_R_.AD")
        sample_a_ad = str(sample + "_A_.AD")
        index_sample = df[((df[sample_a_ad] + df[sample_r_ad]) < 10)|(df[sample_a_ad] <3) | (df[sample_r_ad] <3)].index
        """ Uncomment this if you have a way to verify the genotypes
        if i == 0: #For the first sample
            indexes_ad10 = np.array(index_sample)
            #print("Number of rows to drop", len(indexes_ad10))
            i = i + 1
        elif i != 0:#We accumulate the indexes to drop in the next samples
            index_sample = np.array(index_sample) # We need to convert into an array for compare with the array ad10_index
            # The next is the array of the common and unique elements in both arrays (common indexes or common rows mean the indexes of the rows to drop)
            indexes_ad10 = np.intersect1d(indexes_ad10, index_sample)
            print("Number of rows to drop", len(indexes_ad10))
            i = i + 1
        """
    df.drop(index=indexes_ad10, inplace=True)
    return df


def main():
    df_AD10 = AD10(df_average, samples)
    # Write the temporary file
    df_AD10.to_csv(ad10_df)



if __name__ == '__main__':
    main()


shell(
    "{log}"
)
