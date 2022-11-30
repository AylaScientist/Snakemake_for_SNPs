"""Snakemake wrapper for create pseudogenomes"""

__author__ = "Aurora Campo"
__copyright__ = "Copyright 2021, Aurora Campo, modified from git"
__git__= "https://github.com/johanzi/make_pseudogenome"
__email__ = "ayla.bcn@gmail.com"
__license__ = "MIT"


from snakemake.shell import shell
#from snakemake_wrapper_utils.java import get_java_opts

input = snakemake.input.get("vcf")
extra = snakemake.params.get("extra")
ref = snakemake.input.get("ref")
pseudo = snakemake.output.get("pseudo")
log = snakemake.log_fmt_shell(stdout=False, stderr=True)


shell(
    "bcftools consensus "
    "{input} "
    "--sample {extra} "#merged_1GF_filtered_ref
    "--fasta-ref {ref} "
    "> {pseudo} "
    "{log}"
)
