# python 3.7
# Stats.py


"""
python 3.7

    @autor : Ayla
    @version : 0.1
    @date : July 2020

In this script we will process one table containing SNPs that have been processed with the script "ASE_workflow.py" and
whose quality has been checked for mapping bias with the script "QC.py".
These SNPs are unbiased and biallelic, and the monoallelic expression (MAE) has been discarded.
This script will process the data for testing the factors described in the experimental design.
It will also draw a heatmap with the allele frequencies of the SNPs.
It needs a file with the name of the experimental group as the acronyms used in the sample name. For example:
    The sample number 1 corresponds to the gills exposed to freshwater
    The two experimental factors are "gills" (G) and "freshwater" (F)
    Then the sample is called "1GF" and the group name will be called "GF"
Another file is also the experimental design, including the groups that will be tested against each other:
    The first column must include the test_ID (Test_1, Test_2,...)
    The second column must include the acronyms for the control group in the test
    The next columns must include the acronyms for the other groups to test against the control

This script can be modified to adapt for other type of statistical analysis such as LME or RL.
"""

# Import the libraries
import pandas as pd
from scipy import stats
from scipy.stats import chisquare
import scipy as sc
import matplotlib.pyplot as plt

import os
import statsmodels.stats as smt
from statsmodels.stats import multitest
from snakemake.shell import shell
from snakemake_wrapper_utils.java import get_java_opts
shell.executable("bash")

extra = snakemake.params.get("extra")
java_opts = snakemake.params.get("java_opts")
log = snakemake.log_fmt_shell(stdout=False, stderr=True)


def chi_test(df, af_samples, experiment, exp):
    # Collect the row of the experiment where the design is described
    design = exp[exp["Test_ID"].str.match(experiment)].iloc[0]
    group1 = design["Group_1"]
    group2 = design["Group_2"]
    # Collect the columns of the dataframe that will go into the test:
    group1_samples = [i for i in af_samples if group1 in i]
    group2_samples = [i for i in af_samples if group2 in i]
    if len(group1_samples) > len(group2_samples):
        print("WARNING: CHI2 Test for ", experiment,
              " not possible with the group distribution.\nThe number of samples will adjust for performing the test")

        # Pop the extra samples:
        to_pop = len(group1_samples) - len(group2_samples)
        for j in range(to_pop):
            group1_samples.pop(j)

        cols_control = df[df.columns[df.columns.isin(group1_samples)]]
        cols_exp = df[df.columns[df.columns.isin(group2_samples)]]
        # Initialize the lists for collecting the statstics
        chi_group2 = []
        p_chi_group2 = []
        stat_chi_group2 = []

        # Perform the test
        for i in range(len(df)):
            chi_group2.append(chisquare(cols_exp.iloc[i,], f_exp=cols_control.iloc[i,]))
            p_chi_group2.append(chi_group2[i].pvalue)
            stat_chi_group2.append(chi_group2[i].statistic)

        # add the multitest also
        multi = smt.multitest.multipletests(p_chi_group2, alpha=0.05, method='bonferroni', is_sorted=False,
                                            returnsorted=False)
        # Append the p-values
        pvalue = experiment + "_CHI_p-val"
        pfdr = experiment + "_CHI_p-fdr"
        df[pvalue] = p_chi_group2
        df[pfdr] = multi[1]

    elif len(group1_samples) < len(group2_samples):
        print("WARNING: CHI2 Test for ", experiment,
              " not possible with the group distribution.\nThe number of samples will adjust for performing of the test")

        # Pop the extra samples
        to_pop = len(group2_samples) - len(group1_samples)
        for j in range(to_pop):
            group2_samples.pop(j)
        cols_control = df[df.columns[df.columns.isin(group1_samples)]]
        cols_exp = df[df.columns[df.columns.isin(group2_samples)]]
        # Initialize the lists for collecting the statstics
        chi_group2 = []
        p_chi_group2 = []
        stat_chi_group2 = []

        # Perform the test
        for i in range(len(df)):
            chi_group2.append(chisquare(cols_exp.iloc[i,], f_exp=cols_control.iloc[i,]))
            p_chi_group2.append(chi_group2[i].pvalue)
            stat_chi_group2.append(chi_group2[i].statistic)

        # add the multitest also
        multi = smt.multitest.multipletests(p_chi_group2, alpha=0.05, method='bonferroni', is_sorted=False,
                                            returnsorted=False)

        # Append the p-values
        pvalue = experiment + "_CHI_p-val"
        pfdr = experiment + "_CHI_p-fdr"
        df[pvalue] = p_chi_group2
        df[pfdr] = multi[1]

    elif (len(group1_samples) < 5) | (len(group2_samples) < 5):
        print("The number of counts for ", experiment," is too low and the chi_square test cannot be performed")
        print("The number of samples in each group must be at least 5")

    else:
        cols_control = df[df.columns[df.columns.isin(group1_samples)]]
        cols_exp = df[df.columns[df.columns.isin(group2_samples)]]
        # Initialize the lists for collecting the statstics
        chi_group2 = []
        p_chi_group2 = []
        stat_chi_group2 = []

        # Perform the test
        for i in range(len(df)):
            chi_group2.append(chisquare(cols_exp.iloc[i,], f_exp=cols_control.iloc[i,]))
            p_chi_group2.append(chi_group2[i].pvalue)
            stat_chi_group2.append(chi_group2[i].statistic)

        # add the multitest also
        multi = smt.multitest.multipletests(p_chi_group2, alpha=0.05, method='bonferroni', is_sorted=False,
                                            returnsorted=False)
        # Append the p-values
        pvalue = experiment + "_CHI_p-val"
        pfdr = experiment + "_CHI_p-fdr"
        df[pvalue] = p_chi_group2
        df[pfdr] = multi[1]

    return df


def chi_square(df, samples, exp, experiments):
    # Create the strings corresponding to the allele frequency in each sample
    global df_chi
    af_samples = []
    for sample in samples:
        af_sample = "AF_" + sample
        af_samples.append(af_sample)

    # Elements for the test
    for experiment in experiments:
        df_chi = chi_test(df, af_samples, experiment, exp)
    print("Chi^2 test performed if possible on all the experiments\nOtherwise, Fisher test will be used for the final report of the results.")

    return df_chi




def fisher_test(df, samples, exp, experiment):
    # Prepare the strings for the column names
    design = exp[exp["Test_ID"].str.match(experiment)].iloc[0]
    group1 = design["Group_1"]
    group2 = design["Group_2"]

    # Collect the columns of the dataframe that will go into the test:
    group1_samples = [i for i in samples if group1 in i]
    group2_samples = [i for i in samples if group2 in i]
    r = "_R_.AD"
    a = "_A_.AD"


    if len(group1_samples) == len(group2_samples):

        # Perform the test
        k = 0  # Iterator for the samples in group 2
        for group1_sample in group1_samples:
            p_fisher = []
            stat_fisher = []
            for i in range(len(df)):
                # sample names and positions
                control_r = df.iloc[i, int(df.columns.get_loc(group1_sample + r))]
                control_a = df.iloc[i, int(df.columns.get_loc(group1_sample + a))]
                exp_r = df.iloc[i, int(df.columns.get_loc(group2_samples[k] + r))]
                exp_a = df.iloc[i, int(df.columns.get_loc(group2_samples[k] + a))]

                # test
                oddsratio, pvalue = stats.fisher_exact([[exp_r, exp_a], [control_r, control_a]])
                stat_fisher.append(oddsratio)
                p_fisher.append(pvalue)

            # add the multitest also
            multi = smt.multitest.multipletests(p_fisher, alpha=0.05, method='fdr_bh', is_sorted=False,
                                                returnsorted=False)

            # Make columns of the dataframe with the results
            fisher_pval = experiment + "_Fisher_p-val"
            fisher_fdr = experiment + "_Fisher_p-fdr"
            df[fisher_pval] = pd.Series(p_fisher, index = df.index)
            df[fisher_fdr] = pd.Series(multi[1], index = df.index)
            k = (k + 1)


    elif len(group1_samples) > len(group2_samples):
        print("WARNING: Fisher Test for ", experiment,
              " not possible with the group distribution.\nThe number of samples will adjust for performing of the test")

        # Pop the extra samples
        to_pop = len(group1_samples) - len(group2_samples)
        for j in range(to_pop):
            group1_samples.pop(j)

        # Perform the test
        k = 0  # Iterator for the samples in group 2
        for group1_sample in group1_samples:
            p_fisher = []
            stat_fisher = []
            for i in range(len(df)):
                # sample names and positions
                control_r = df.iloc[i, int(df.columns.get_loc(group1_sample + r))]
                control_a = df.iloc[i, int(df.columns.get_loc(group1_sample + a))]
                exp_r = df.iloc[i, int(df.columns.get_loc(group2_samples[k] + r))]
                exp_a = df.iloc[i, int(df.columns.get_loc(group2_samples[k] + a))]

                # test
                oddsratio, pvalue = stats.fisher_exact([[exp_r, exp_a], [control_r, control_a]])
                stat_fisher.append(oddsratio)
                p_fisher.append(pvalue)


            # add the multitest also
            multi = smt.multitest.multipletests(p_fisher, alpha=0.05, method='fdr_bh', is_sorted=False,
                                                returnsorted=False)

            # Make columns of the dataframe with the results
            fisher_pval = experiment + "_" + group1_sample + group2_samples[k] + "_Fisher_p-val"
            fisher_fdr = experiment + "_" + group1_sample + group2_samples[k] + "_Fisher_p-fdr"
            df[fisher_pval] = pd.Series(p_fisher, index=df.index)
            df[fisher_fdr] = pd.Series(multi[1], index=df.index)
            k = k + 1

    else:
        print("WARNING: Fisher Test for ", experiment,
              " not possible with the group distribution.\nThe number of samples will adjust for performing of the test")

        # Pop the extra samples
        to_pop = len(group2_samples) - len(group1_samples)
        for j in range(to_pop):
            group2_samples.pop(j)

        # Perform the test
        k = 0  # Iterator for the samples in group 2
        for group1_sample in group1_samples:
            p_fisher = []
            stat_fisher = []
            for i in range(len(df)):
                # sample names and positions
                control_r = df.iloc[i, int(df.columns.get_loc(group1_sample + r))]
                control_a = df.iloc[i, int(df.columns.get_loc(group1_sample + a))]
                exp_r = df.iloc[i, int(df.columns.get_loc(group2_samples[k] + r))]
                exp_a = df.iloc[i, int(df.columns.get_loc(group2_samples[k] + a))]

                # test
                oddsratio, pvalue = stats.fisher_exact([[exp_r, exp_a], [control_r, control_a]])
                stat_fisher.append(oddsratio)
                p_fisher.append(pvalue)

            # add the multitest also
            multi = smt.multitest.multipletests(p_fisher, alpha=0.05, method='fdr_bh', is_sorted=False,
                                                returnsorted=False)

            # Make columns of the dataframe with the results
            fisher_pval = experiment + "_" + group1_sample + group2_samples[k] + "_Fisher_p-val"
            fisher_fdr = experiment + "_" + group1_sample + group2_samples[k] + "_Fisher_p-fdr"
            df[fisher_pval] = pd.Series(p_fisher, index=df.index)
            df[fisher_fdr] = pd.Series(multi[1], index=df.index)
            k = k + 1



    return df


def Fisher(df, samples, exp, experiments):
    for experiment in experiments:
        df_fisher = fisher_test(df, samples, exp, experiment)
    print("Fisher exact test successfully performed on  all the experiments")
    return df_fisher


def binomial(df,samples):
    r = "_R_.AD"
    a = "_A_.AD"
    for sample in samples:
        pval = []
        fdr_pval = []

        for i in range(len(df)):
            # Get the values of the rows
            ref = df.iloc[i, int(df.columns.get_loc(sample + r))]
            alt = df.iloc[i, int(df.columns.get_loc(sample + a))]
            total_reads = ref + alt

            # Perform the test
            binom = (sc.stats.binom_test(ref, total_reads, 0.5, alternative='two-sided'))
            pval.append(binom)

        # add the multitest also
        multi = smt.multitest.multipletests(pval, alpha=0.05, method='fdr_bh', is_sorted=False,
                                            returnsorted=False)
        col_name_binom = sample + "_Binomial_pvalue"
        col_name_multi = sample + "Binomial_fdr_pvalue"
        df[col_name_binom] = pd.Series(pval, index=df.index)
        df[col_name_multi] = pd.Series(multi[1], index=df.index)
    print("Binomial test successfully performed on all the experiments")
    return df


def main ():
    
    # Import the files
    exp = pd.read_csv(snakemake.input.get("i1"))
    sample_names = pd.read_csv(snakemake.input.get("i3"))
    df_qc = pd.read_csv(snakemake.input.get("i4"))

    groups_df = pd.read_csv(snakemake.input.get("i2"))
    # Create a numpy array of arrays with the samples of each group and the total samples of the experiment
    groups = groups_df.values
    # Create a numpy array with the name of each group
    group_names = list(groups_df["Group"])

    # Create arrays of the sample names
    samples = sample_names["Sample_name"].values

    # Create an array with the name of each experimental test
    experiments = list(exp["Test_ID"])


    """
    STATISTICAL TESTS
    -----------------
    We implement:
        Chi² test for allele specific expression in each experimental group based in allele frequencies
        Binomial test for allelic imbalance in each sample based in allele counts
        Fisher exact test for allele specific expression in samples from the same individual based in allele counts
    A Bonferroni one-step correction will be applied to the p-values of the chi-square test and a Benjamini-Hochberg
    non-negative test will be applied to the Fisher and Binomial tests respectively.


    Chi-square test for conditions based on allele frequencies
    ----------------------------------------------------------
    In this test we compare the allele frequencies in allelic imbalance according to the experimental design described in
    the file "Experimental_design.csv". It is possible to test more than one factor but it will compare them by pairs, being
    the "group1" the control group and the "group2" the experimental.

    The function used is "chisquare" from the package scipy:
        scipy.stats.chisquare(f_obs, f_exp=None, ddof=0, axis=0)[source]¶

        Calculate a one-way chi-square test.

        The chi-square test tests the null hypothesis: the categorical data has the expected frequencies.

        Parameters

            f_obsarray_like
                Observed frequencies in each category.

            f_exparray_like, optional
                Expected frequencies in each category. By default the categories are assumed to be equally likely.

            ddofint, optional
                “Delta degrees of freedom”: adjustment to the degrees of freedom for the p-value. The p-value is computed using
                a chi-squared distribution with k - 1 - ddof degrees of freedom, where k is the number of observed frequencies.
                The default value of ddof is 0.

            axisint or None, optional
                The axis of the broadcast result of f_obs and f_exp along which to apply the test. If axis is None, all values
                in f_obs are treated as a single data set. Default is 0.


        Returns

            chisqfloat or ndarray
                The chi-squared test statistic. The value is a float if axis is None or f_obs and f_exp are 1-D.
            pfloat or ndarray
                The p-value of the test. The value is a float if ddof and the return value chisq are scalars.

    Chi-square for each described test
    ----------------------------------
    This test will compare the differences between the two experimental groups described in the test. For this, the
    dataframe can't have allele frequencies equal to 0 for all the members of the experimental group. The formula needs the
    next variables:
        sampling size
        expected probability (reference allele frequency from control group)
        observed probability (reference allele frequency from the experimental group)

    NOTE: The number of samples in the control and experimental groups have to be the same. Otherwise the script will
    pop the extra samples in the biggest group. Also the number of samples in each group must be at least 5

    Now let's perform the test.
    """

    #os.chdir("temp")
    df_chi = chi_square(df_qc, samples, exp, experiments)

    df_chi.to_csv(snakemake.output.get("o1"))

    """
    Fisher exact test for ASE
    -------------------------
    This test is based on the read counts of each allele, including alternative and reference alleles.

    """

    df_fisher = Fisher(df_chi, samples, exp, experiments)

    df_fisher.to_csv(snakemake.output.get("o2"))

    """
    Binomial test
    -------------
    The binomial test will determine the allelic imbalance in each sample by comparing the counts of the reference allele
    with the total counts of the SNP. This test does not compare experimental groups, but each sample.

    Perform the binomial tests in the read counts of the reference allele (it can be done by
    experimental group or by individual)(x) over the total number of reads in this group (as before, it can be done by
    experimental group or by individual) (n). For this analysis, the expected frequency is 0.5 and the observed allele
    frequency is determined for the reference allele in each individual or for the control group (p). This test will
    provide the p-values for the significant differences of the groups. The formula is extracted from Wikipedia and
    adapted to each test. In python the formula is scipy.stats.binom_test (x,n,p, alternative "two-sided").
    """

    df_binomial = binomial(df_fisher, samples)

    #os.chdir(cwd)
    df_binomial.to_csv(snakemake.output.get("o3"))
    print("Analyzed data available in " + snakemake.output.get("o3"))

if __name__ == '__main__' :
    main()
