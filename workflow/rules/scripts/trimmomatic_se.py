"""
bio/trimmomatic/pe
Snakemake wrapper to trim reads with trimmomatic in SE mode with help of pigz.
pigz is the parallel implementation of gz. Trimmomatic spends most of the time
compressing and decompressing instead of trimming sequences. By using process
substitution (<(command), >(command)), we can accelerate trimmomatic a lot.
Consider providing this wrapper with at least 1 extra thread per each gzipped
input or output file.
"""

__author__ = "Aurora Campo"
__copyright__ = "Copyright 2022, Aurora Campo"
__email__ = "ayla.bcn@gmail.com"
__license__ = "MIT"


import os
from snakemake.shell import shell
from snakemake_wrapper_utils.java import get_java_opts


# Get the variables
#extra = snakemake.params.get("extra")
java_opts = get_java_opts(snakemake)
log = snakemake.log_fmt_shell(stdout=True, stderr=True)
trimmer =snakemake.params.get("trimmer")
threads = snakemake.params.get("threads")
input_r1 = snakemake.input.get("r1")
print(input_r1)
output_r1 = snakemake.output.get("r1")
path = snakemake.params.get("path")


shell(
    "java -XX:ParallelGCThreads={threads} {java_opts} -jar {path} SE "#-threads {threads} {java_opts} "#"{extra} "
    "-phred64 "
	"./{input_r1} "
    "./{output_r1} "
    "{trimmer} "#"ILLUMINACLIP:TruSeq3-PE-2.fa:2:30:10:8:TRUE"
    "{log}"
)

