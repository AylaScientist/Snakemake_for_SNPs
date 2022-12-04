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


extra = snakemake.params.get("extra")
java_opts = snakemake.params.get("java_opts")
log = snakemake.log_fmt_shell(stdout=False, stderr=True)


df_bi=snakemake.input.get("csv")
df_bi = pd.read_csv(df_bi, low_memory=False)
df_bi = pd.DataFrame(df_bi)


sample_names = pd.read_csv(snakemake.input.get("sn1"))
PSG_codes = pd.read_csv(snakemake.input.get("psc"))
# Create arrays of the sample names and the pseudogenome codes
samples = sample_names['Sample_name'].values
PSGs = PSG_codes['PSGs'].values

df_av=snakemake.output[0]


# Definitions:
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
            y = df.iloc[i, aPSG2_AD] / 2
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
            y = df.iloc[i, rPSG2_AD] / 2
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
            y = df.iloc[i, aPSG1_AD] / 2
            a_AD.append(y)

        # Case 7 Second homozygot, first heterozygot with alternative allele in first position (as reference)
        elif (df.iloc[i, rPSG1_GT] != df.iloc[i, aPSG1_GT]) and (df.iloc[i, rPSG2_GT] == df.iloc[i, aPSG2_GT]) and (
                df.iloc[i, aPSG1_GT] == df.iloc[i, aPSG2_GT]):
            #print("Alt GT: ", df.iloc[i, rPSG2_GT], " Counts: ", df.iloc[i, aPSG1_AD], " ", df.iloc[i, rPSG2_AD], " ",df.iloc[i, aPSG2_AD], " averaged by 2")
            r_GT.append(df.iloc[i, rPSG1_GT])
            x = df.iloc[i, rPSG1_AD] / 2
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

def main():
    df_average, cols = sample_average(df_bi, samples, PSGs)
    # Drop the columns that are not needed
    df_average.drop(df_average.columns[cols], axis=1, inplace=True)
    df_average.to_csv(df_av)


if __name__ == '__main__':
    main()


shell(
    "{log}"
)
