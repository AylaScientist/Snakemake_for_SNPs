__author__ = "Aurora Campo"
__copyright__ = "Copyright 2020, Aurora Campo"
__email__ = "ayla.bcn@gmail.com"
__license__ = "CC"


from snakemake.shell import shell
from snakemake_wrapper_utils.java import get_java_opts

extra = snakemake.params.get("extra", "")
java_opts = get_java_opts(snakemake)

log = snakemake.log_fmt_shell(stdout=True, stderr=True)
shell(
    "gatk --java-options '{java_opts}' VariantsToTable -R {snakemake.input.ref} -V {snakemake.input.vcf} "
    "{extra} -O {snakemake.output.vcf} {log}"
)
