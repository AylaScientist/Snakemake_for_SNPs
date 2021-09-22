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

extra = snakemake.params.get("extra")
java_opts = snakemake.params.get("java_opts")

csv1=snakemake.input.get("csv1")
csv2=snakemake.input.get("csv2")
merged_df=snakemake.output[0]

PSG1 = pd.read_csv(csv1, low_memory=False)
PSG2 = pd.read_csv(csv2, low_memory=False)
PSG1 = pd.DataFrame(PSG1)
PSG2 = pd.DataFrame(PSG2)

df = pd.merge({PSG1}, {PSG2}, on=('CHROM', 'POS'))
df.to_csv({merged_df}, sept=',')
print('Files merged')

shell(
    "{log}"
)
