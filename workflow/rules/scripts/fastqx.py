"""Snakemake wrapper for fastqc."""

__author__ = "Aurora Campo"
__copyright__ = "Copyright 2020, Aurora Campo, modified from Copyright 2017, Julian de Ruiter"
__email__ = "ayla.bcn@gmail.com"
__license__ = "CC"


from os import path
from tempfile import TemporaryDirectory

from snakemake.shell import shell

log = snakemake.log_fmt_shell(stdout=False, stderr=True)


def basename_without_ext(file_path):
    """Returns basename of file path, without the file extension."""

    base = path.basename(file_path)

    split_ind = 2 if base.endswith(".fastq.gz") else 1
    base = ".".join(base.split(".")[:-split_ind])

    return base


# Run fastqc, since there can be race conditions if multiple jobs
# use the same fastqc dir, we create a temp dir.
with TemporaryDirectory() as tempdir:
    shell(
        "fastx_trimmer {snakemake.params} -t {snakemake.threads} "
        "--outdir {tempdir:q} {snakemake.input[0]:q}"
        " {log:q}"
    )

    # Move outputs into proper position.
    output_base = basename_without_ext(snakemake.input[0])
    html_path = path.join(tempdir, output_base + "_fastqc.html")
    zip_path = path.join(tempdir, output_base + "_fastqc.zip")

    if snakemake.output.html != html_path:
        shell("mv {html_path:q} {snakemake.output.html:q}")

    if snakemake.output.zip != zip_path:
        shell("mv {zip_path:q} {snakemake.output.zip:q}")
