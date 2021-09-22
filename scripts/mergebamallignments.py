"""Snakemake wrapper for picard MergeBamAllignments."""

__author__ = "Aurora Campo"
__copyright__ = "Copyright 2021, Aurora Campo"
__email__ = "ayla.bcn@gmail.com"
__license__ = "MIT"


from snakemake.shell import shell
from snakemake_wrapper_utils.java import get_java_opts

extra = snakemake.params.get("extra", "")
#java_opts = get_java_opts(snakemake)

unmmapped = snakemake.input.get("ubam")
mapped = snakemake.input.get("mbam")
out = snakemake.output.get("merged")

assert out is not None, "output-> Output is a required parameter"
print(out)

ref = snakemake.input.get("ref")
log = snakemake.log_fmt_shell(stdout=False, stderr=True)

shell(
    "picard"
    " MergeBamAlignment"
    " UNMAPPED={unmmapped}"
    " ALIGNED={mapped}"
    " O={out}"
    " R={ref}"
    " {extra}"
    " {log}"
)
