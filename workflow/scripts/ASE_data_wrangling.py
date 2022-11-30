# Python 3.7
# ASE_workflow.py

"""
python 3.7

    @autor : Ayla
    @version : 0.1
    @date : July 2020

In this script we will process two tables containing SNPs that have been called against two sample-based pseudogenomes
from the same species. The input is produced after annotation of the vcf file with annovar and extraction of the table
using the "Variant to table" tool from the GATK package. The variants are biallelic and the structure of the csv file
includes four columns for each sample: allele depths from recessive and alternative alleles and genotypes from recessive
and alternative alleles.

NOTE_ PSG stands for pseudogenome
"""

import pandas as pd
import numpy as np
import os
from os import path





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


def multiallelic(df, samples, PSGs, multi_index):
    for sample in samples:
        multiallelic_sample_index = multiallelic_sample(df, sample, PSGs)
        multi_index = multi_index + multiallelic_sample_index
    df_bi = df.drop(index=multi_index)
    return df_bi


def evaluation(df, sample, PSGs):
    # Evaluation of genotype and average. This is an iterator
    r = "_R_"
    a = "_A_"
    gt = ".GT"
    ad = ".AD"
    PSG1_code: str = PSGs[0]
    PSG2_code: str = PSGs[1]
    rPSG1_GT = int(df.columns.get_loc(str(sample + r + PSG1_code + gt)))
    aPSG1_GT = int(df.columns.get_loc(str(sample + a + PSG1_code + gt)))
    rPSG2_GT = int(df.columns.get_loc(str(sample + r + PSG2_code + gt)))
    aPSG2_GT = int(df.columns.get_loc(str(sample + a + PSG2_code + gt)))
    rPSG1_AD = int(df.columns.get_loc(str(sample + r + PSG1_code + ad)))
    aPSG1_AD = int(df.columns.get_loc(str(sample + a + PSG1_code + ad)))
    rPSG2_AD = int(df.columns.get_loc(str(sample + r + PSG2_code + ad)))
    aPSG2_AD = int(df.columns.get_loc(str(sample + a + PSG2_code + ad)))

    # Create the empty lists to collect the average variables and genotypes
    r_AD = []
    a_AD = []
    r_GT = []
    a_GT = []

    print("Evaluation loop in sample: ", sample)

    # Start iterator
    for i in range(len(df)):
        # Case 1 homozygots for the two pseudogenoes
        if (df.iloc[i, rPSG1_GT] == df.iloc[i, aPSG1_GT]) and (df.iloc[i, rPSG2_GT] == df.iloc[i, aPSG2_GT]) and (
                df.iloc[i, rPSG1_GT] == df.iloc[i, rPSG2_GT]) and (df.iloc[i, aPSG1_GT] == df.iloc[i, aPSG2_GT]):
            #print("Ref GT: ", df.iloc[i, rPSG1_GT], " Counts: ", df.iloc[i, rPSG1_AD], " ", df.iloc[i, aPSG1_AD], " ", df.iloc[i, rPSG2_AD],  " ", df.iloc[i, aPSG2_AD]," averaged by 2")
            r_GT.append(df.iloc[i, rPSG1_GT])
            if (df.iloc[i, rPSG1_AD] != 0 and df.iloc[i, aPSG1_AD] !=0):
                x = ((((df.iloc[i, rPSG1_AD] + df.iloc[i, aPSG1_AD])/2) + df.iloc[i, rPSG2_AD] + df.iloc[i, aPSG2_AD]) / 2)
            elif (df.iloc[i, rPSG2_AD] != 0 and df.iloc[i, aPSG2_AD] !=0):
                x = ((((df.iloc[i, rPSG2_AD] + df.iloc[i, aPSG2_AD]) / 2) + df.iloc[i, rPSG1_AD] + df.iloc[i, aPSG1_AD]) / 2)
            elif (df.iloc[i, rPSG1_AD] != 0 and df.iloc[i, aPSG1_AD] !=0 and df.iloc[i, rPSG2_AD] != 0 and df.iloc[i, aPSG2_AD] != 0):
                x = ((((df.iloc[i, rPSG2_AD] + df.iloc[i, aPSG2_AD]) / 2) + (df.iloc[i, rPSG1_AD] + df.iloc[
                    i, aPSG1_AD])/2) / 2)
            else:
                x = ((df.iloc[i, rPSG1_AD] + df.iloc[i, rPSG2_AD] + df.iloc[i, aPSG1_AD] + df.iloc[i, aPSG2_AD]) / 2)
            r_AD.append(x)
            a_GT.append(df.iloc[i, rPSG2_GT])
            y = 0
            a_AD.append(y)

        # Case 2 Same situation Ref and Alt in both for df
        elif (df.iloc[i, rPSG1_GT] == df.iloc[i, rPSG2_GT]) and (df.iloc[i, aPSG1_GT] == df.iloc[i, aPSG2_GT]):
            r_GT.append(df.iloc[i, rPSG1_GT])
            x = (df.iloc[i, rPSG1_AD] + df.iloc[i, rPSG2_AD]) / 2
            r_AD.append(x)
            a_GT.append(df.iloc[i, aPSG1_GT])
            y = (df.iloc[i, aPSG1_AD] + df.iloc[i, aPSG2_AD]) / 2
            a_AD.append(y)

        # Case 3 Inverse situation Ref and Alt in both for df
        elif (df.iloc[i, rPSG1_GT] == df.iloc[i, aPSG2_GT]) and (df.iloc[i, aPSG1_GT] == df.iloc[i, rPSG2_GT]):
            r_GT.append(df.iloc[i, rPSG1_GT])
            x = (df.iloc[i, rPSG1_AD] + df.iloc[i, aPSG2_AD]) / 2  # averaged by 2 (for each mapping method)
            r_AD.append(x)
            a_GT.append(df.iloc[i, aPSG1_GT])
            y = (df.iloc[i, aPSG1_AD] + df.iloc[i, rPSG2_AD]) / 2  # averaged by 2 (for each mapping method)
            a_AD.append(y)

        # Case 4 First homozygot, second with alternative allele in second position (as alternative)
        elif (df.iloc[i, rPSG1_GT] == df.iloc[i, aPSG1_GT]) and (df.iloc[i, rPSG2_GT] != df.iloc[i, aPSG2_GT]) and (
                df.iloc[i, rPSG1_GT] == df.iloc[i, rPSG2_GT]):
            # print("Ref GT: ", df.iloc[i, rPSG1_GT]," Counts: ",df.iloc[i, rPSG1_AD]," ",df.iloc[i, aPSG1_AD]," ",df.iloc[i, rPSG2_AD]," averaged by 2")
            r_GT.append(df.iloc[i, rPSG1_GT])
            if (df.iloc[i, rPSG1_AD] != 0 and df.iloc[i, aPSG1_AD] != 0 and df.iloc[i, rPSG2_AD] != 0):
                #print ("WARNING: Line ",i," in sample ",sample," has three counts for the same allele. It is corrected by averaging by 3")
                x = (df.iloc[i, rPSG1_AD] + df.iloc[i, aPSG1_AD] + df.iloc[i, rPSG2_AD]) / 3
            else:
                x = (df.iloc[i, rPSG1_AD] + df.iloc[i, aPSG1_AD] + df.iloc[i, rPSG2_AD]) / 2  # counts from 3 cells from
                # which one is 0, from 2 mapping methods (averaged by 2)
            r_AD.append(x)
            a_GT.append(df.iloc[i, aPSG2_GT])
            y = df.iloc[i, aPSG2_AD]
            a_AD.append(y)

        # Case 5 First homozygot, second with alternative allele in first position  (as reference)
        elif ((df.iloc[i, rPSG1_GT] == df.iloc[i, aPSG1_GT]) and (df.iloc[i, rPSG2_GT] != df.iloc[i, aPSG2_GT]) and (
                df.iloc[i, rPSG1_GT] == df.iloc[i, aPSG2_GT])):
            #print("Ref GT: ", df.iloc[i, rPSG1_GT], " Counts: ", df.iloc[i, rPSG1_AD], " ", df.iloc[i, aPSG1_AD], " ",df.iloc[i, aPSG2_AD], " averaged by 2")
            r_GT.append(df.iloc[i, rPSG1_GT])
            if (df.iloc[i, rPSG1_AD] != 0 and df.iloc[i, aPSG1_AD] != 0 and df.iloc[i, aPSG2_AD]  != 0) :
                #print("WARNING: Line ", i, " in sample ", sample, " has three counts for the same allele. It is corrected by averaging by 3")
                x = (df.iloc[i, rPSG1_AD] + df.iloc[i, aPSG1_AD] + df.iloc[i, aPSG2_AD]) / 3
            else:
                x = (df.iloc[i, rPSG1_AD] + df.iloc[i, aPSG1_AD] + df.iloc[i, aPSG2_AD]) / 2  # counts from 3 cells from which one is 0, from 2 mapping methods (averaged by 2)
            r_AD.append(x)
            a_GT.append(df.iloc[i, rPSG2_GT])
            y = df.iloc[i, rPSG2_AD]
            a_AD.append(y)

        # Case 6 Second homozygot, first heterozygot with alternative allele in second position (as alternative)
        elif (df.iloc[i, rPSG1_GT] != df.iloc[i, aPSG1_GT]) and (df.iloc[i, rPSG2_GT] == df.iloc[i, aPSG2_GT]) and (
                df.iloc[i, rPSG1_GT] == df.iloc[i, rPSG2_GT]):
            #print("Ref GT: ", df.iloc[i, rPSG1_GT], " Counts: ", df.iloc[i, rPSG1_AD], " ", df.iloc[i, rPSG2_AD], " ",df.iloc[i, aPSG2_AD], " averaged by 2")
            r_GT.append(df.iloc[i, rPSG1_GT])
            if (df.iloc[i, rPSG1_AD] != 0 and df.iloc[i, rPSG2_AD] != 0 and df.iloc[i, aPSG2_AD] !=0 ):
                #print("WARNING: Line ", i, " in sample ", sample, " has three counts for the same allele. It is corrected by averaging by 3")
                x = (df.iloc[i, rPSG1_AD] + df.iloc[i, rPSG2_AD] + df.iloc[i, aPSG2_AD]) / 3
            else:
                x = (df.iloc[i, rPSG1_AD] + df.iloc[i, rPSG2_AD] + df.iloc[i, aPSG2_AD]) / 2  # counts from 3 cells from which one is 0, from 2 mapping methods (averaged by 2)
            r_AD.append(x)
            a_GT.append(df.iloc[i, aPSG1_GT])
            y = df.iloc[i, aPSG1_AD]
            a_AD.append(y)

        # Case 7 Second homozygot, first heterozygot with alternative allele in first position (as reference)
        elif (df.iloc[i, rPSG1_GT] != df.iloc[i, aPSG1_GT]) and (df.iloc[i, rPSG2_GT] == df.iloc[i, aPSG2_GT]) and (
                df.iloc[i, aPSG1_GT] == df.iloc[i, aPSG2_GT]):
            #print("Alt GT: ", df.iloc[i, rPSG2_GT], " Counts: ", df.iloc[i, aPSG1_AD], " ", df.iloc[i, rPSG2_AD], " ",df.iloc[i, aPSG2_AD], " averaged by 2")
            r_GT.append(df.iloc[i, rPSG1_GT])
            x = df.iloc[i, rPSG1_AD]
            r_AD.append(x)
            a_GT.append(df.iloc[i, rPSG2_GT])
            if (df.iloc[i, aPSG1_AD] != 0 and df.iloc[i, rPSG2_AD] != 0 and df.iloc[i, aPSG2_AD] !=0):
                #print("WARNING: Line ", i, " in sample ", sample, " has three counts for the same allele. It is corrected by averaging by 3")
                y = (df.iloc[i, aPSG1_AD] + df.iloc[i, rPSG2_AD] + df.iloc[i, aPSG2_AD]) / 3
            else:
                y = (df.iloc[i, aPSG1_AD] + df.iloc[i, rPSG2_AD] + df.iloc[i, aPSG2_AD]) / 2  # counts from 3 cells from which one is 0, from 2 mapping methods (averaged by 2)
            a_AD.append(y)

        # Case 8 Both homozygot, first for the reference allele and second for the alternative allele
        elif (df.iloc[i, rPSG1_GT] == df.iloc[i, aPSG1_GT]) and (df.iloc[i, rPSG2_GT] == df.iloc[i, aPSG2_GT]) and (
                df.iloc[i, aPSG1_GT] != df.iloc[i, aPSG2_GT]):
            #print("Ref GT: ", df.iloc[i, rPSG1_GT], " Counts: ", df.iloc[i, rPSG1_AD], " ", df.iloc[i, aPSG1_AD], "averaged by 2")

            r_GT.append(df.iloc[i, rPSG1_GT])
            x = df.iloc[i, rPSG1_AD] + df.iloc[i, aPSG1_AD]  # one of these counts will be 0 so no need to average
            r_AD.append(x)
            a_GT.append(df.iloc[i, rPSG2_GT])
            y = df.iloc[i, rPSG2_AD] + df.iloc[i, aPSG2_AD]  # one of these counts will be 0 so no need to average
            a_AD.append(y)

        else:
            print("Error in the SNP on the coordinates ", df.iloc[i, "CHROM"], ", ", df.iloc[i, "POS"],
                    ", sample num ",
                    sample)

    return (r_GT, r_AD, a_GT, a_AD)


def sample_average(df, samples, PSGs):
    cols = []  # Index of columns to be drop in the end
    # Create 4 columns for each averaged sample:
    for sample in samples:
        print("Sample average in sample ", sample, "\n")
        # Average each sample according to genotype:
        (r_GT, r_AD, a_GT, a_AD) = evaluation(df, sample, PSGs)

        # Add the new columns
        sample_r_GT = str(sample + "_R_.GT")
        sample_a_GT = str(sample + "_A_.GT")
        sample_r_AD = str(sample + "_R_.AD")
        sample_a_AD = str(sample + "_A_.AD")
        df[sample_r_GT] = r_GT
        df[sample_a_GT] = a_GT
        df[sample_r_AD] = r_AD
        df[sample_a_AD] = a_AD

        # Drop the columns that won't be used
        r = "_R_"
        a = "_A_"
        gt = ".GT"
        ad = ".AD"
        PSG1_code: str = PSGs[0]
        PSG2_code: str = PSGs[1]
        rPSG1_GT = int(df.columns.get_loc(str(sample + r + PSG1_code + gt)))
        aPSG1_GT = int(df.columns.get_loc(str(sample + a + PSG1_code + gt)))
        rPSG2_GT = int(df.columns.get_loc(str(sample + r + PSG2_code + gt)))
        aPSG2_GT = int(df.columns.get_loc(str(sample + a + PSG2_code + gt)))
        rPSG1_AD = int(df.columns.get_loc(str(sample + r + PSG1_code + ad)))
        aPSG1_AD = int(df.columns.get_loc(str(sample + a + PSG1_code + ad)))
        rPSG2_AD = int(df.columns.get_loc(str(sample + r + PSG2_code + ad)))
        aPSG2_AD = int(df.columns.get_loc(str(sample + a + PSG2_code + ad)))

        cols.append(rPSG1_GT)
        cols.append(aPSG1_GT)
        cols.append(rPSG2_GT)
        cols.append(aPSG2_GT)
        cols.append(rPSG1_AD)
        cols.append(aPSG1_AD)
        cols.append(rPSG2_AD)
        cols.append(aPSG2_AD)

    return df, cols


def AD10(df, samples):
    # The iterator collects the index were AD from both alleles is below 10 or the AD for one allele is below 3
    # Keep all this analysis for the same individual including two tissues at a time
    i = 0
    indexes_ad10 = 0
    x = str(path.exists("Samples_MAE.csv"))
    if (x == True):
        # Verification that there are two tissues
        tissues = pd.read_csv(snakemake.input.get("MAE")) #Tissues to determine MAE
        tissues = pd.DataFrame(tissues)
        first_samples = tissues.iloc[:, 0]
        second_samples = tissues.iloc[:, 1]

        for first_sample in first_samples:
            print("AD10 filter in sample ", first_sample, " and sample ", second_samples[i])
            sample1_r_ad = str(first_sample + "_R_.AD")
            sample1_a_ad = str(first_sample + "_A_.AD")
            sample2_r_ad = str(second_samples[i] + "_R_.AD")
            sample2_a_ad = str(second_samples[i] + "_A_.AD")
            # Collect the index <10 for the total of the alleles in the sample and the index <3 for one of the alleles.
            index_sample = df[
                (((df[sample1_a_ad] + df[sample1_r_ad]) < 10) | ((df[sample2_a_ad] + df[sample2_r_ad]) < 10) | (df[sample1_a_ad] < 3) | (df[sample1_r_ad]  < 3) | (df[sample2_a_ad] <3) | (df[sample2_r_ad] <3))].index
            if i == 0:
                indexes_ad10 = np.array(index_sample)
                print("Number of rows to drop for sample ", first_sample, " ", len(indexes_ad10))
                i = i + 1
            elif i != 0:
                index_sample = np.array(index_sample)
                indexes_ad10 = np.intersect1d(indexes_ad10, index_sample)
                print("Number accumulated of rows to drop for sample ", first_sample, " ", len(indexes_ad10))
                i = i + 1

    else:
        # There is only one tissue
        for sample in samples:
            print("AD10 filter in sample ", sample)
            sample_r_ad = str(sample + "_R_.AD")
            sample_a_ad = str(sample + "_A_.AD")

            index_sample = df[((df[sample_a_ad] + df[sample_r_ad]) < 10)|(df[sample_a_ad] <3) | (df[sample_r_ad] <3)].index
            if i == 0:
                indexes_ad10 = np.array(index_sample)
                print("Number of rows to drop", len(indexes_ad10))
                i = i + 1
            elif i != 0:
                index_sample = np.array(index_sample)
                indexes_ad10 = np.intersect1d(indexes_ad10, index_sample)
                print("Number of rows to drop", len(indexes_ad10))
                i = i + 1

    df.drop(index=indexes_ad10, inplace=True)
    return df


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
    SNP_panel.to_csv(snakemake.output.get("result2"))

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
            a_AD.append(df.iloc[i, sample_a_ad])a.bool(), a.item(), a.any() or a.all().
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
            genotype_models.to_csv(snakemake.output.get("result3"))

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


def MAE(df, samples):
    # Collect the indexes whose allele frequencies are 0 or 1 in both tissues of the same individual. These indexes
    # correspond to the SNPs that are MAE for at least one sample
    x = str(path.exists("Samples_MAE.csv"))
    if (x == 'True'):
        # Case where there are two tissues compared
        tissues = pd.read_csv(snakemake.input("MAE"))
        tissues = pd.DataFrame(tissues)
        first_samples = tissues.iloc[:, 0]
        second_samples = tissues.iloc[:, 1]
        i = 0
        for first_sample in first_samples:
            print("MAE loop in sample ", first_sample, " and sample ", second_samples[i])
            af1 = str("AF_" + first_sample)
            af2 = str("AF_" + second_samples[i])
            index_af = df[((df[af1] == 0) & (df[af2] == 0)) | (df[af1] == 1) & (df[af2] == 1)].index
            df.drop(index=index_af, inplace=True)
            i = i + 1
    else:
        # Case where there is only one tissue
        for sample in samples:
            af = str("AF_" + sample)
            index_af = df[(df[af] == 0) | (df[af] == 1)].index
            df.drop(index=index_af, inplace=True)

    return df


def main():
    """
    Input the files and convert them into data frames
    """

    # Read the dataframe and convert it into pandas dataframe. PSG stands for pseudogenome:

    #PATH = os.getcwd ()

    """
    PSG1_name: str = input("Please enter the file name of the SNPs_for_wrangling for the first pseudogenome:\n")
    PSG2_name: str = input("Please enter the file name of the SNPs_for_wrangling for the second pseudogenome:\n")

    PSG1 = pd.read_csv(PSG1_name, low_memory=False)
    PSG2 = pd.read_csv(PSG2_name, low_memory=False)

    """
    PSG1 = pd.read_csv(snakemake.input.get("csv1"), low_memory=False)
    PSG2 = pd.read_csv(snakemake.input.get("csv2"), low_memory=False)

    PSG1 = pd.DataFrame(PSG1)
    PSG2 = pd.DataFrame(PSG2)

    # Read the sample names and the codes for each pseudogenome:
    """
    sample_names = pd.read_csv(input("Please enter the file name that has the sample names:\n"))
    PSG_codes = pd.read_csv(input("Please enter the file name that has the pseudogenomes codes:\n"))"
    """

    sample_names = pd.read_csv(snakemake.input.get("sn1"))
    PSG_codes = pd.read_csv(snakemake.input.get("psc"))

    # Create arrays of the sample names and the pseudogenome codes
    samples = sample_names['Sample_name'].values
    PSGs = PSG_codes['PSGs'].values

    """
     MERGE THE DATAFRAMES
    -----------------------

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

    The next step is to merge the two dataframes on the common chromosome and position:
    """

    df = pd.merge(PSG1, PSG2, on=('CHROM', 'POS'))
    print('Files merged')

    """
     DELETE MULTIALLELIC SITES
     -------------------------

    Compare the genotypes for assessing the counts on each allele
    Set the data for the mapping against the pseudogenome of sample PSG1 because it hast the variants of the reference genotype. The SNPs called by this sample are set as reference and alternative alleles for all the other samples.

    The tables containing the SNPs from mapping to PSG1 and PSG2 contained only biallelic sites. However multiallelic possibilities can appear if one SNP was called for a different genotpe in each of the pseudogenomes. The multiallelic sites will be deleted from this dataframe to continue the analysis with bialleleic sites only.

    Delete multiallelic sites
    Let's start droping the multiallelic. The criteria compares if the reference allele in the SNPs resulting from mapping against PSG1 pseudogenome is different form the reference and the alternative alleles from the mapping against PSG2 pseudogenome. This is later repeated for the alternative allele. Note that the SNPs mapped against PSG1 pseudogenome have been pre-filtered for biallelic sites for all the samples.

    Collect the indexes where the reference allele in PSG1 is different from both alleles in the other mappings or the alternative allele in PSG1 is different from both alleles in the PSG2 mapping.
    """

    # Call the function to clean the multiallelic sites and drop the dfs with this condition:
    multi_index = []
    df_bi = multiallelic(df, samples, PSGs, multi_index)

    # Make a temporary folder where you can copy temporary files for backups
    # os.mkdir("temp")
    #os.chdir("temp")
    df_bi.to_csv(snakemake.output.get("temp1"))
    #os.chdir(PATH)

    # From now on we work with df_bi standing for "dataframe_biallelic_SNPs"

    """
    COMPARE GENOTYPES AND MAKE THE AVERAGE OF COUNTS
    ------------------------------------------------

    Now let's make the average on the counts for the reference and the alternative alleles. For each sample, select the
    genotype of the reference allele established by the reads mapped against PSG1 and compare it to both alternative and
    reference genotypes of the values in the mapping with PSG2 genome. Then proceed to make the average of the counts. The
    possible cases are the next:

    1.- Both homozigots for the reference allele: In that case set the same genotype to both alleles and sum all the
    counts from both sites, making the average for two individuals.

    2.- Reference and alternative allele in the SNP from the PSG1 pseudogenome are placed in the same position for the SNP
    from the PSG2 pseudogenome. Set the genotype to the reference and alternative alleles to the ones from the SNP in PSG1
    and make the average of the counts for reference and alternative genome

    3.- The reference allele in the reference pseudogenome is called as alternative allele in the alternative
    pseudogenome. In that case, set the genotype to the reference allele and make the average of the counts for the
    corresponding variant.

    4.- While the mapping against the reference genome resulted in an homozygot for the reference allele, the mapping
    against the alternative genome resulted in an heterozygot with the alternative allele in the alternative position.
    The genotype of the reference allele will be set from genotype of the reference allele in the reference mapping and
    the alternative allele will be set from the genotype of the alternative allele in the alternative mapping. The counts
    of each allele will be set accordingly to this distribution, thus dividing the counts by two individuals in the
    reference allele.

    5.- The mapping against the reference genome resulted in an homozygot for the reference allele and the mapping
    against the alternative genome is heterozygot. In that case, the genotype of the alternative allele is set as the
    reference allele on the mapping against the alternative genome. The genotypes must be set accordingly,
    being the reference allele the one found in the reference genome and the alternative allele the one found in the
    alternative genome. The average of the counts will follow this pattern, being the average in the erference allele.

    6.- The mapping against the reference genome produced an heterozygot with reference and alternative alleles. The
    mapping against the alternative pseudogenome produced an homozygot for the reference allele. In that case,
    reference and alternative allele must be set according to the reference pseudogenome and the counts must be averaged
    for both individuals in the alternative allele.

    7.- The mapping against the reference genome produced an heterozygot with reference and alternative alleles. The
    mapping against the alternative pseudogenome produced an homozygot for the alternative allele. In that case,
    reference and alternative allele must be set according to the reference pseudogenome and the counts of the
    alternative allele must be averaged for both individuals.

    8.- The mapping against the reference genome produced an homozygot for the reference allele. The mapping against the
    alterantive pseudogenome produced an homozygot for the alternative allele. In that case, the reference allele is set
    according to the first pseudogenome and summing the counts of this reference allele and the alternative allele is set
    following the second pseudogenome and the counts of this alternative allele are summed.

    This analysis is repeated for each sample and will assign reference and alternative alleles independently. The
    assignation of reference and alternative alleles uniformly in all the samples will be developed in a further step of
    the workflow.
    """

    df_average, cols = sample_average(df_bi, samples, PSGs)

    # Drop the columns that are not needed
    df_average.drop(df_average.columns[cols], axis=1, inplace=True)

    #os.chdir("temp")
    df_average.to_csv(snakemake.output.get("temp2"))
    #os.chdir(PATH)

    """
    DELETE SNPs WITH AD<10
    ----------------------
    This new function cleans the SNPs that do not have enough counts and are considered possible poor quality reads.
    The SNPs with AD < 3 in one of the alleles are also cleaned.
    """

    df_AD10 = AD10(df_average, samples)

    # Write the temporary file
    #os.chdir("temp")
    df_AD10.to_csv(snakemake.output.get("temp3"))
    #os.chdir(PATH)

    """
        ASSIGN REFERENCE AND ALTERNATIVE ALLELES UNIFORMLY IN ALL SAMPLES
        -----------------------------------------------------------------
        Reference and alternative alleles have been assigned by sample. This function will correct it and will assign the
        reference and alternative alleles according to the pattern established in sample 1.
        """

    df_uni = genotype(df_AD10, samples)

    # Write the temporary file
    #os.chdir("temp")
    df_uni.to_csv(snakemake.output.get("result1"))
    #os.chdir(PATH)

    """
    ELIMINATE MONOALLELIC EXPRESSION
    --------------------------------
    The monoallelic expression (MAE) can result from homozygots as well as imprinted genes. In order to distinguish each
    case, we need to know the genotype. However, that's not possible for the amount of SNPs that we are working with.
    Therefore we drop all the SNPs whose allele frequency is 0 or 1.

    If the samples are taken from two different organs from the same individual, the system needs a file called
    "Samples_MAE.csv" where the sample name from an individual are in the same row and there are two columns,
    one for each tissue/organ we need to compare. The name of the columns must be in plural for further iterations.
    """

    df_af = frequencies(df_uni, samples)

    df_mae = MAE(df_af, samples)

    df_mae.to_csv(snakemake.output.get("temp4"))


if __name__ == '__main__':
    main()
