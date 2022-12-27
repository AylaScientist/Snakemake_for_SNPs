__author__ = "Aurora Campo"
__copyright__ = "Copyright 2021, Aurora Campo"
__email__ = "ayla.bcn@gmail.com"
__license__ = "MIT"


import os
import pandas as pd
import numpy as np
from os import path
from os.path import exists
from snakemake.shell import shell
from snakemake_wrapper_utils.java import get_java_opts

shell.executable("bash")

"""
FILTER SNPs that are wrongly called by using the SNPs from the genome call.
"""

extra = snakemake.params.get("extra")
java_opts = snakemake.params.get("java_opts")
log = snakemake.log_fmt_shell(stdout=False, stderr=True)




def main():
    # Input the files and convert them into data frames
    SNPs2val = snakemake.input.get("result1")
    valid = snakemake.output.get("result4")
    error = snakemake.output.get("error")

    SNPs2val = pd.read_csv(SNPs2val, low_memory=False)
    SNPs2val = pd.DataFrame(SNPs2val)

    # 
    if str(os.path.exists("config/AD_GT_counts_bi_DNA.csv"))==True:
        geno_df = pd.read_csv("config/AD_GT_counts_bi_DNA.csv", low_memory=False)
        geno_df = pd.DataFrame(geno_df)
        coincident = pd.merge(SNPs2val, geno_df, how="inner", on=["CHROM", "POS"])
        print("SNPs retrieved from the pipeline: ",len(SNPs2val),"; Valid SNPs: ", len(coincident))
        coincident.to_csv(valid)
        error = ((len(SNPs2val) - len(coincident)) * 100) / len(SNPs2val)
        print("Error: ", error)
        data = np.array([["SNPs retrived from the pipeline",len(SNPs2val)],
                        ["Valid SNPs",len(coincident)],
                        ["Error in SNP calling",error]])
        # {"SNPs retrived from the pipeline":len(SNPs2val),"Valid SNPs":len(coincident),"Error in SNP calling":error}
        df = pd.DataFrame(data)
        df.to_csv(error)
    else:
        SNPs2val=coincident
        coincident.to_csv(valid)
        print("Error: You need to call the SNPs from DNA in order to get the data for validation.")
        text_file = open(error, "w")
        text_file.write("Error: You need to call the SNPs from DNA in order to get the data for validation.")
        text_file.close()


if __name__ == '__main__':
    main()


shell(
    "{log}"
)
