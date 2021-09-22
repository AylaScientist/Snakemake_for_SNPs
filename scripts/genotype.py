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
ASSIGN REFERENCE AND ALTERNATIVE ALLELES UNIFORMLY IN ALL SAMPLES
-----------------------------------------------------------------
Reference and alternative alleles have been assigned by sample. This function will correct it and will assign the
reference and alternative alleles according to the pattern established in sample 1.
"""


extra = snakemake.params.get("extra")
java_opts = snakemake.params.get("java_opts")

ad10_df=snakemake.input.get("csv")
geno_df=snakemake.output.get("result1")
result2=snakemake.output.get("result2")
result3=snakemake.output.get("result3")

df_ad10 = pd.read_csv(ad10_df, low_memory=False)
df_ad10 = pd.DataFrame(df_ad10)
sample_names = pd.read_csv(snakemake.input.get("sn1"))
# Create arrays of the sample names
samples = sample_names['Sample_name'].values



def compare(df, sample, R_models, A_models,SNP_panel):
    # This iterator will construct the new columns for the samples according to the reference and alternative alleles
    # assigned in sample1

    # Position of the sample columns
    sample_r_gt = int(df.columns.get_loc(str(sample + "_R_.GT")))
    sample_a_gt = int(df.columns.get_loc(str(sample + "_A_.GT")))
    sample_r_ad = int(df.columns.get_loc(str(sample + "_R_.AD")))
    sample_a_ad = int(df.columns.get_loc(str(sample + "_A_.AD")))

    # Create the empty lists to collect the average variables and genotypes
    r_AD = []
    a_AD = []
    r_GT = []
    a_GT = []

    print("Compare genotypes loop in sample: ", sample)

    # Start iterator for model genotypes
    for k in range(len(df)):
        # Set in the model the new genotypes in case the samples didn't express this allele.
        if R_models[k] == A_models[k] and R_models[k] == ".":
            # print("Previous: ", R_models[k], A_models[k], ", new: ", df.iloc[k, sample_r_gt], df.iloc[k, sample_a_gt])
            R_models[k] = df.iloc[k, sample_r_gt]
            A_models[k] = df.iloc[k, sample_a_gt]
        elif R_models[k] == A_models[k] and R_models[k] == df.iloc[k, sample_r_gt]:
            A_models[k] = df.iloc[k, sample_a_gt]
        elif R_models[k] == A_models[k] and R_models[k] == df.iloc[k, sample_a_gt]:
            A_models[k] = df.iloc[k, sample_r_gt]

    # Make a dataframe with the model genotypes
    SNP_panel["Reference allele"] = R_models
    SNP_panel["Alternative allele"] = A_models
    SNP_panel.to_csv(result2)

    # Start iterator for comparison with the new models
    for i in range(len(df)):
        if R_models[i] == df.iloc[i, sample_r_gt] and A_models[i] == df.iloc[i, sample_a_gt]:
            # 1 Reference and alternative are set in the same position for sample and sample1 (reference)
            # It applies to both samples heterozygots or both homozygots with same genotype
            r_AD.append(df.iloc[i, sample_r_ad])
            a_AD.append(df.iloc[i, sample_a_ad])
            r_GT.append(df.iloc[i, sample_r_gt])
            a_GT.append(df.iloc[i, sample_a_gt])
        elif A_models[i] == df.iloc[i, sample_r_gt] and R_models[i] == df.iloc[i, sample_a_gt]:
            # 2 Reference and alternative are set in the inverse position for sample and sample1 (reference)
            # It applies to both samples heterozygots or both homozygots with same genotype
            r_AD.append(df.iloc[i, sample_a_ad])
            a_AD.append(df.iloc[i, sample_r_ad])
            r_GT.append(df.iloc[i, sample_a_gt])
            a_GT.append(df.iloc[i, sample_r_gt])
        elif R_models[i] == df.iloc[i, sample_r_gt] and A_models[i] == df.iloc[i, sample_r_gt] and A_models[i] != \
                df.iloc[i, sample_a_gt]:
            # 3 It applies to first sample homozygot and second sample heterozygot with the alternative allele in second
            r_AD.append(df.iloc[i, sample_r_ad])
            a_AD.append(df.iloc[i, sample_a_ad])
            r_GT.append(df.iloc[i, sample_r_gt])
            a_GT.append(df.iloc[i, sample_a_gt])
        elif R_models[i] != df.iloc[i, sample_r_gt] and A_models[i] != df.iloc[i, sample_r_gt] and A_models[i] == \
                df.iloc[
                    i, sample_a_gt]:
            # 4 It applies to first sample homozygot and second sample heterozygot with the alternative allele in first
            r_AD.append(df.iloc[i, sample_a_ad])
            a_AD.append(df.iloc[i, sample_r_ad])
            r_GT.append(df.iloc[i, sample_a_gt])
            a_GT.append(df.iloc[i, sample_r_gt])
        elif R_models[i] == A_models[i] and R_models[i] != df.iloc[i, sample_r_gt] and df.iloc[i, sample_r_gt] == \
                df.iloc[i, sample_a_gt]:
            # 5 It applies to first sample homozygot for the reference allele and second sample homozygot for the
            # alternative allele
            a_AD.append(df.iloc[i, sample_r_ad] + df.iloc[i, sample_a_ad])  # One of them will be 0
            r_GT.append(df.iloc[i, sample_a_gt])
            a_GT.append(df.iloc[i, sample_a_gt])
            r_AD.append(0)
        elif R_models[i] != A_models[i] and R_models[i] == df.iloc[i, sample_r_gt] and df.iloc[i, sample_r_gt] == \
                df.iloc[i, sample_a_gt]:
            # 6 It applies to first sample heterozygot and second sample homozygot for the reference allele
            r_GT.append(df.iloc[i, sample_r_gt])
            a_GT.append(df.iloc[i, sample_r_gt])
            r_AD.append(df.iloc[i, sample_r_ad] + df.iloc[i, sample_r_ad])  # One of them will be 0
            a_AD.append(0)
        elif R_models[i] != A_models[i] and A_models[i] == df.iloc[i, sample_a_gt] and df.iloc[i, sample_r_gt] == \
                df.iloc[i, sample_a_gt]:
            # 7 It applies to first sample heterozygot and second sample homozygot for the alternative allele
            a_AD.append(df.iloc[i, sample_r_ad] + df.iloc[i, sample_a_ad])  # One of them will be 0
            r_GT.append(df.iloc[i, sample_a_gt])
            a_GT.append(df.iloc[i, sample_a_gt])
            r_AD.append(0)
        elif R_models[i] == df.iloc[i, sample_r_gt] and A_models[i] == df.iloc[i, sample_r_gt] and df.iloc[
            i, sample_r_gt] == df.iloc[i, sample_a_gt]:
            # 8 It applies to first sample homozygot for the reference allele and second sample homozygot for the
            # reference allele too
            r_AD.append(df.iloc[i, sample_r_ad] + df.iloc[i, sample_a_ad])  # One of them will be 0
            r_GT.append(df.iloc[i, sample_r_gt])
            a_GT.append(df.iloc[i, sample_r_gt])
            a_AD.append(0)
        elif R_models[i] == A_models[i] and R_models[i] == ".":
            # 9 SNP not expressed in the sample1 as reference gentoype
            r_AD.append(df.iloc[i, sample_r_ad])
            a_AD.append(df.iloc[i, sample_a_ad])
            r_GT.append(df.iloc[i, sample_r_gt])
            a_GT.append(df.iloc[i, sample_a_gt])
            #print(i)

            genotype_models = pd.DataFrame(data=[R_models, A_models]).T
            genotype_models.to_csv(result3)

        elif df.iloc[i, sample_r_gt] == df.iloc[i, sample_a_gt] and df.iloc[i, sample_r_gt] == ".":
            # 10 SNP not expressed in the sample to evaluate
            r_AD.append(0)
            a_AD.append(0)
            r_GT.append(".")
            a_GT.append(".")
        else:
            print("ERROR: Comparison of genotypes of sample ",
                    sample, " provided an error for the SNP in the row ", i, " with genotypes ",
                    df.iloc[i, sample_a_gt], ", ",
                    df.iloc[i, sample_r_gt], ". You must check if the evaluation of this case is correct")


    return r_AD, a_AD, r_GT, a_GT, R_models, A_models



def genotype(df, samples):
    # Make the new df and prepare the first columns of the SNP_panel:
    df_uni = df.iloc[:, 1:7]
    SNP_panel = df.iloc[:, 1:7]
    # Determines the genotypes for sample1 and calls the function to compare genotypes
    j = 0
    for sample in samples:
        j = j + 1
        if j == 1:  # If we are in the first sample we set the reference
            R_models = df[str(sample + "_R_.GT")].to_numpy()
            A_models = df[str(sample + "_A_.GT")].to_numpy()
            # Collection the name of the samples
            sample_r_gt = str(sample + "_R_.GT")
            sample_a_gt = str(sample + "_A_.GT")
            sample_r_ad = str(sample + "_R_.AD")
            sample_a_ad = str(sample + "_A_.AD")
            # Assign the list to the column name
            df_uni[sample_r_gt] = df[str(sample + "_R_.GT")]
            df_uni[sample_a_gt] = df[str(sample + "_A_.GT")]
            df_uni[sample_r_ad] = df[str(sample + "_R_.AD")]
            df_uni[sample_a_ad] = df[str(sample + "_A_.AD")]
        else:  # If we already passed the first sample
            (r_ad, a_ad, r_gt, a_gt, R_models, A_models) = compare(df, sample, R_models, A_models, SNP_panel)
            # Collection the name of the samples
            sample_r_gt = str(sample + "_R_.GT")
            sample_a_gt = str(sample + "_A_.GT")
            sample_r_ad = str(sample + "_R_.AD")
            sample_a_ad = str(sample + "_A_.AD")
            # Assign the list to the column name
            df_uni[sample_r_gt] = r_gt
            df_uni[sample_a_gt] = a_gt
            df_uni[sample_r_ad] = r_ad
            df_uni[sample_a_ad] = a_ad

    return df_uni

def main():
    df_uni = genotype(df_AD10, samples)
    df_uni.to_csv(geno_df)

if __name__ == '__main__':
main()


shell(
    "{log}"
)
