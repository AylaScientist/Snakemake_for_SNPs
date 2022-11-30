__author__ = "Aurora Campo"
__copyright__ = "Copyright 2021, Aurora Campo"
__email__ = "ayla.bcn@gmail.com"
__license__ = "MIT"


from snakemake.shell import shell
#from snakemake_wrapper_utils.java import get_java_opts

input = snakemake.input.get("bam")
ref = snakemake.input.get("annotation")
output = snakemake.output[0]
extra = snakemake.params.get("extra")
log = snakemake.log_fmt_shell(stdout=False, stderr=True)


shell(
     "htseq-count "
     "{extra} "
     "{input} "
     "{ref} "
     "> {output} "
     "{log}"
)
