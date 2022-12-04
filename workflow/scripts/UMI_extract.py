__author__ = "Aurora Campo"
__copyright__ = "Copyright 2021, Aurora Campo, modified from Copyright 2016, Johannes Köster"
__email__ = "ayla.bcn@gmail.com"
__license__ = "MIT"


import os
from snakemake.shell import shell

extra = snakemake.params.get("extra", "")
log = snakemake.log_fmt_shell(stdout=True, stderr=True)

fq1 = snakemake.input.get("fq1")
assert fq1 is not None, "input-> fq1 is a required input parameter"
fq1 = (
    [snakemake.input.fq1]
    if isinstance(snakemake.input.fq1, str)
    else snakemake.input.fq1
    )
fq2 = snakemake.input.get("fq2")
if fq2:
    fq2 = (
        [snakemake.input.fq2]
        if isinstance(snakemake.input.fq2, str)
        else snakemake.input.fq2
        )

assert len(fq1) == len(
        fq2
        ), "input-> equal number of files required for fq1 and fq2"

if fq1[0].endswith(".gz"):
    readcmd = "--readFilesCommand zcat"
else:
    readcmd = ""

pro1 = snakemake.output.get("pro1")
assert pro1 is not None, "output-> pro1 is a required output parameter"
pro1 = (
    [snakemake.output.pro1]
    if isinstance(snakemake.output.pro1, str)
    else snakemake.output.pro1
    )
"""
pro2 = snakemake.output.get("pro2")
pro2 = (
    [snakemake.output.pro2]
    if isinstance(snakemake.output.pro2, str)
    else snakemake.output.pro2
    )
"""

# For paired corresponds
shell(
    "umi_tools "
    "extract "
    "--extract-method=tag "
    "--bc-pattern=NNNXXXXNN "
    "--stdin={fq2} "
    "--read2-in={fq1} "
    "--read2-out={pro1} "
    "--log={log}"
)
