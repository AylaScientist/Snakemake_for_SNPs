__author__ = "Aurora Campo"
__copyright__ = "Cpoyright 2020, Aurora Campo, modified from Copyright 2018, Johannes Köster"
__email__ = "ayla.bcn@gmail.com; johannes.koester@protonmail.com"
__license__ = "MIT"


import os

from snakemake.shell import shell
#from snakemake_wrapper_utils.java import get_java_opts

extra = snakemake.params.get("extra")
java_opts = snakemake.params.get("java_opts")

known = snakemake.input.get("known", "")
if known:
    known = "--dbsnp " + known

bams = snakemake.input.bam
if isinstance(bams, str):
    bams = [bams]
bams = list(map("-I {}".format, bams))

log = snakemake.log_fmt_shell(stdout=True, stderr=True)
shell(
    "gatk --java-options '{java_opts}' HaplotypeCaller {extra} "
    "-R {snakemake.input.ref} {bams} "
    "-ERC GVCF "
    "-O {snakemake.output.gvcf} {known} {log}"
)
