__author__ = "Patrik Smeds adapted by Aurora Campo"
__copyright__ = "Copyright 2021, Aurora Campo"
__email__ = "ayla.bcn@gmail.com"
__license__ = "MIT"


from snakemake.shell import shell
#from snakemake_wrapper_utils.java import get_java_opts

shell.executable("bash")

#java_opts = get_java_opts(snakemake)
extra_params = snakemake.params.get("extra", "")
bam_input = snakemake.input.get("bam")
umi_input = snakemake.input.get("umi")
output_file = snakemake.output.get("bam")
log = snakemake.log_fmt_shell(stdout=False, stderr=True)


if bam_input is None:
    raise ValueError("Missing bam input file!")
elif not isinstance(bam_input, str):
    raise ValueError("Input bam should be a string: " + str(bam_input) + "!")

if umi_input is None:
    raise ValueError("Missing input file with UMIs")
elif not isinstance(umi_input, str):
    raise ValueError("Input UMIs-file should be a string: " + str(umi_input) + "!")

"""
if not len(output_file) == 1:
    raise ValueError("Only one output value expected: " + str(output_file) + "!")
"""
if output_file is None:
    raise ValueError("Missing output file!")
"""
elif not isinstance(output_file, str):
    raise ValueError("Output bam-file should be a string: " + str(output_file) + "!")
"""


shell(
    "java -jar /data/bin/miniconda2/envs/fgbio-v1.3.0/share/fgbio/fgbio.jar "
    "AnnotateBamWithUmis "
    "-i {bam_input} "
    "-f {umi_input} "
    "-o {output_file} "
    "{extra_params} "
    "{log}"
)
