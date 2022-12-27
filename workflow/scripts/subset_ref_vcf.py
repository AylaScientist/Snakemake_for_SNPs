"""Snakemake wrapper for create pseudogenomes"""

__author__ = "Aurora Campo"
__copyright__ = "Copyright 2021, Aurora Campo, modified from git"
__git__= "https://github.com/johanzi/make_pseudogenome"
__email__ = "ayla.bcn@gmail.com"
__license__ = "MIT"


from snakemake.shell import shell
#from snakemake_wrapper_utils.java import get_java_opts

vcf = snakemake.input.get("vcf")
out = snakemake.output.get("vcf")
outvcf = out.split(".")[0]
outvcfgz = outvcf + ".gz"
log = snakemake.log_fmt_shell(stdout=False, stderr=True)

shell(
    # Subset by keeping positions with GQ >= 30 and DP >=5
    # Not keeping indels
    "vcftools "
    "--vcf {vcf} "
    "--remove-indels "
    "--minGQ 30 "
    "--minDP 5 "
    "--recode "
    "--recode-INFO-all "
    "--out {outvcf} "
    "{log}"
)
