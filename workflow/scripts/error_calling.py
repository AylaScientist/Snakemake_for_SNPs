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
Calculate SNPs that are wrongly called by using the SNPs from the genome call for each ID.
"""


extra = snakemake.params.get("extra")
java_opts = snakemake.params.get("java_opts")
log = snakemake.log_fmt_shell(stdout=False, stderr=True)

"""
Input the files and convert them into data frames
"""
SNPs2val=pd.read_csv(snakemake.input.get("result2"))
geno_df=pd.read_csv(snakemake.input.get("geno"))
sample_names = pd.read_csv(snakemake.input.get("sn1"))

error=snakemake.output.get("error_calling")


SNPs2val = pd.DataFrame(SNPs2val)
geno_df = pd.DataFrame(geno_df)
samples = sample_names['Sample_name'].values



def main():

    coincident_2_9 = pd.merge(SNPs2val, geno_df, how="inner", on=["CHROM", "POS"])
    print("SNPs retrieved from the pipeline: ",len(SNPs2val),"; Valid SNPs: ", len(coincident_2_9))
    coincident_2_9.to_csv(valid)
    error2_9 = ((len(SNPs2val) - len(coincident_2_9)) * 100) / len(SNPs2val)
    print("Error: ", error2_9)
    data = np.array([["SNPs retrived from the pipeline",len(SNPs2val)],
                    ["Valid SNPs",len(coincident_2_9)],
                    ["Error in SNP calling",error2_9]])
    # {"SNPs retrived from the pipeline":len(SNPs2val),"Valid SNPs":len(coincident_2_9),"Error in SNP calling":error2_9}

    df = pd.DataFrame(data)
    df.to_csv(error)

if __name__ == '__main__':
    main()


shell(
    "{log}"
)
