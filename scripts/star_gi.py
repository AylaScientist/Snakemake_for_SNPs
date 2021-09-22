"""Snakemake wrapper for indexing a genome for STAR mapping"""

__author__ = "Aurora Campo"
__copyright__ = "Copyright 2021, Aurora Campo"
__email__ = "ayla.bcn@gmail.com"
__license__ = "MIT"

from snakemake.shell import shell
#from snakemake_wrapper_utils.java import get_java_opts

input = snakemake.input[0]
threads = snakemake.params.get("threads")
annotation = snakemake.params.get("annotation")
dir = snakemake.params.get("dir")
read_length = int(snakemake.params.get("read_length"))-1
log = snakemake.log_fmt_shell(stdout=False, stderr=True)

shell(
    "STAR "
    "--runThreadN {threads} "
    "--runMode genomeGenerate "
    "--genomeDir {dir} "
    "--genomeFastaFiles {input} "
    "--sjdbGTFfile {annotation} "
    "--sjdbOverhang {read_length} " #ReadLength-1
    "{log}"
)
