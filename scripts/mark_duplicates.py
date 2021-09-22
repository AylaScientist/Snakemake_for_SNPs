__author__ = "Aurora Campo"
__copyright__ = "Copyright 2020, Aurora Campo modified from Copyright 2016, Johannes Köster"
__email__ = "ayla.bcn@gmail.com"
__license__ = "CC"


from snakemake.shell import shell
#from snakemake_wrapper_utils.java import get_java_opts

log = snakemake.log_fmt_shell(stdout=True, stderr=True)

extra = snakemake.params.get("extra")
java_opts = snakemake.params.get("java_opts")


shell(
    "picard MarkDuplicates "  # Tool and its subcommand
    "{java_opts} "  # Automatic java option
    "{extra} "  # User defined parmeters
    "INPUT={snakemake.input} "  # Input file
    "OUTPUT={snakemake.output.bam} "  # Output bam
    "METRICS_FILE={snakemake.output.metrics} "  # Output metrics
    "{log}"  # Logging
)
