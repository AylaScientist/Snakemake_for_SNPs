__author__ = "Aurora Campo"
__copyright__ = "Copyright 2021, Aurora Campo"
__email__ = "ala.bcn@gmail.com"
__license__ = "MIT"


import os
from snakemake.shell import shell

log = snakemake.log_fmt_shell(stdout=False, stderr=True)
#params = snakemake.params
input = snakemake.input[0]
threads = snakemake.threads

shell(
    "fastqc --quiet "
    "-t {threads} "
    "--outdir ./qc/fastqc/ "
    "-f fastq "
    "{input} "
    "{log}"
)
