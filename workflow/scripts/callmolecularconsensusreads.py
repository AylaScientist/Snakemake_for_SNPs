__author__ = "Patrik Smeds adapted by Aurora Campo"
__copyright__ = "Copyright 2021, Aurora Campo"
__email__ = "ayla.bcn@gmail.com"
__license__ = "MIT"


from snakemake.shell import shell
#from snakemake_wrapper_utils.java import get_java_opts

shell.executable("bash")

#java_opts = get_java_opts(snakemake)
extra_params = snakemake.params.get("extra", "")
log = snakemake.log_fmt_shell(stdout=False, stderr=True)

bam_input = snakemake.input[0]

if not isinstance(bam_input, str) and len(snakemake.input[0]) != 1:
    raise ValueError("Input bam should be one bam file: " + str(bam_input) + "!")

output_file = snakemake.output[0]

if not isinstance(output_file, str) and len(snakemake.output[0]) != 1:
    raise ValueError("Output should be one bam file: " + str(output_file) + "!")

shell(
    "java -jar /data/bin/miniconda2/envs/fgbio-v1.3.0/share/fgbio/fgbio.jar "
    "CallMolecularConsensusReads "
    "-i {bam_input} "
    "-o {output_file} "
    "{extra_params} "
    "{log}"
)
