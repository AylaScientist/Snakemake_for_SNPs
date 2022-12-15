__author__ = "Aurora Campo"
__copyright__ = "Copyright 2021, Aurora Campo, modified from Copyright 2019, Patrik Smeds"
__email__ = "ayla.bcn@gmail.com"
__license__ = "MIT"


import os
from snakemake.shell import shell

log = snakemake.log_fmt_shell(stdout=False, stderr=True)
extra = snakemake.params.get("extra", "")

bam_input = snakemake.input.get("bam")

"""
if bam_input is None:
    raise ValueError("Missing bam input file!")
elif not len(snakemake.input) == 1:
    raise ValueError("Only expecting one input file: " + str(snakemake.input) + "!")
"""

output_file = snakemake.output.get("bam")
groups = snakemake.output.get("hist")

"""
if output_file is None:
    raise ValueError("Missing output file")
elif not len(snakemake.output) == 1:
    raise ValueError("Only expecting one output file: " + str(snakemake.output) + "!")
"""

shell(
    "umi_tools "
    "group "
    "-I {bam_input} "
    "--group-out {groups} "
    "--output-bam "
    "--stdout {output_file} "
    "--umi-group-tag=RX "
    "{extra}"
    "--log={log} "
    "--paired "
)
