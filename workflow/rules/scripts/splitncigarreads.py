__author__ = "Aurora Campo"
__copyright__ = "Copyright 2021, Aurora Campo"
__email__ = "ayla.bcn@gmail.com"
__license__ = "MIT"

import os

from snakemake.shell import shell
from snakemake_wrapper_utils.java import get_java_opts

extra = snakemake.params.get("extra")
java_opts = snakemake.params.get("java_opts")
input = snakemake.input.get("bam")
ref = snakemake.input.get("ref")
output = snakemake.output.get("bam")
log = snakemake.log_fmt_shell(stdout=True, stderr=True)

shell(
    "gatk --java-options '{java_opts}' SplitNCigarReads {extra} "
    "-R {ref} "
    "-I {input} "
    "-O {output} "
    "{log}"
)
