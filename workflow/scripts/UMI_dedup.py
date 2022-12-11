__author__ = "Aurora Campo"
__copyright__ = "Copyright 2021, Aurora Campo, modified from Copyright 2019, Patrik Smeds"
__email__ = "ayla.bcn@gmail.com"
__license__ = "MIT"


import os
from snakemake.shell import shell

log = snakemake.log_fmt_shell(stdout=False, stderr=True)
extra = snakemake.params.get("extra", "")

bam_input = snakemake.input[0]

if bam_input is None:
    raise ValueError("Missing bam input file!")
elif not len(snakemake.input) == 1:
    raise ValueError("Only expecting one input file: " + str(snakemake.input) + "!")

output_file = snakemake.output[0]

if output_file is None:
    raise ValueError("Missing output file")
elif not len(snakemake.output) == 1:
    raise ValueError("Only expecting one output file: " + str(output_file) + "!")


shell(
    "umi_tools "
    "dedup "
    "-I {bam_input} "
    "--paired "
    "-S {output_file} "
    "{log}"
)
