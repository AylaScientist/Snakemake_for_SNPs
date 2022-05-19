__author__ = "Aurora Campo"
__copyright__ = "Copyright 2020, Aurora Campo modified from Copyright 2016, Johannes Köster"
__email__ = "ayla.bcn@gmail.com"
__license__ = "MIT"


from snakemake.shell import shell
#from snakemake_wrapper_utils.java import get_java_opts

extra = snakemake.params.get("extra")
java_opts = snakemake.params.get("java_opts")


shell(
    "picard AddOrReplaceReadGroups {java_opts} {snakemake.params} "
    "I={snakemake.input} "
    "O={snakemake.output} &> {snakemake.log}"
)
