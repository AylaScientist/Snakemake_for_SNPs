"""Snakemake wrapper for picard SortSam."""

__author__ = "Aurora Campo"
__copyright__ = "Copyright 2021, Aurora Campo"
__email__ = "ayla.bcn@gmail.com"
__license__ = "MIT"


from snakemake.shell import shell
#from snakemake_wrapper_utils.java import get_java_opts


extra = snakemake.params.get("extra", "")
#sample = snakemake.params.get("sample")
#java_opts = get_java_opts(snakemake)


genome = snakemake.input[0]
dictionary = snakemake.output[0]
log = snakemake.log_fmt_shell(stdout=False, stderr=True)

shell(
    "picard CreateSequenceDictionary "
    "REFERENCE={genome} "
    "OUTPUT= {dictionary} "
    "{log}"
)
