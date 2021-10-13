"""Snakemake wrapper for picard SortSam."""

__author__ = "Julian de Ruiter"
__copyright__ = "Copyright 2017, Julian de Ruiter"
__email__ = "julianderuiter@gmail.com"
__license__ = "MIT"


from snakemake.shell import shell
#from snakemake_wrapper_utils.java import get_java_opts


extra = snakemake.params.get("extra", "")
sample = snakemake.params.get("sample")
#java_opts = get_java_opts(snakemake)

log = snakemake.log_fmt_shell(stdout=False, stderr=True)

fastq1 = snakemake.input.fastq1
fastq2 = snakemake.input.get("fastq2", None)
#fastq_unpaired = snakemake.output.get("unpaired_fastq", None)

if not isinstance(fastq1, str):
    raise ValueError("f1 needs to be provided")

input = " F1=" + fastq1

if isinstance(fastq2, str):
    input += " F2=" + fastq2

"""
if isinstance(fastq_unpaired, str):
    if not isinstance(fastq2, str):
        raise ValueError("f2 is required if fastq_unpaired is set")

    input += " UNPAIRED_FASTQ=" + fastq_unpaired
"""

shell(
    "picard FastqToSam "
    "{extra} "
    "{input} "
    "O={snakemake.output[0]} "
    "SM={sample} "
)