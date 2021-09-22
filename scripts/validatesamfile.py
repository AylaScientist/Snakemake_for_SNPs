"""Snakemake wrapper for picard SortSam."""

__author__ = "Julian de Ruiter"
__copyright__ = "Copyright 2017, Julian de Ruiter"
__email__ = "julianderuiter@gmail.com"
__license__ = "MIT"


from snakemake.shell import shell
#from snakemake_wrapper_utils.java import get_java_opts


extra = snakemake.params.get("extra", "")
#sample = snakemake.params.get("sample")
#java_opts = get_java_opts(snakemake)

bam = snakemake.input.get("bam")

if not isinstance(bam, str):
    raise ValueError("bam file needs to be provided")

input = " I=" + bam

output = " O=" + snakemake.output[0]
log = snakemake.log_fmt_shell(stdout=False, stderr=True)
ref = " R=" + snakemake.params.get("ref")
mode = " M=" + snakemake.params.get("mode")

shell(
    "picard ValidateSamFile "
    "{input} "
    #"{extra} "
    #"{ref} "
    #"IGNORE_WARNINGS=true "
    "{mode} "
    #"{output} "
    "{log}"
)
