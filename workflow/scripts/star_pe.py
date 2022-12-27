__author__ = "Aurora Campo"
__copyright__ = "Copyright 2020, Aurora Campo, modified from Copyright 2016, Johannes Köster"
__email__ = "ayla.bcn@gmail.com"
__license__ = "MIT"


import os
from snakemake.shell import shell

extra = snakemake.params.get("extra", "")
threads = snakemake.params.get("threads")
log = snakemake.log_fmt_shell(stdout=True, stderr=True)
index = snakemake.input.get("index")
outprefix = snakemake.params.get("filename")
fq1 = snakemake.input.get("fastq1")
fq2 = snakemake.input.get("fastq2")
index_dir = index.rsplit('/', 1)[0] #retains all the string before the last character "/"

assert fq1 is not None, "input-> fq1 is a required input parameter"
fq1 = (
    [snakemake.input.fastq1]
    if isinstance(snakemake.input.fastq1, str)
    else snakemake.input.fastq1
)

if fq2:
    fq2 = (
        [snakemake.input.fastq2]
        if isinstance(snakemake.input.fastq2, str)
        else snakemake.input.fastq2
    )
    assert len(fq1) == len(
        fq2
    ), "input-> equal number of files required for fq1 and fq2"
input_str_fq1 = ",".join(fq1)
input_str_fq2 = ",".join(fq2) if fq2 is not None else ""
input_str = " ".join([input_str_fq1, input_str_fq2])


if fq1[0].endswith(".gz"):
    readcmd = "--readFilesCommand zcat"
else:
    readcmd = ""



shell(
    "STAR "
    "{extra} "
    "--runThreadN {threads} "
    "--genomeDir {index_dir} "
    "--readFilesIn {fq1} {fq2} " #{input_str} for single ends
    "{readcmd} "
    "--outFileNamePrefix {outprefix} "
    "--outStd Log "
    "--twopassMode Basic "
    "--quantMode GeneCounts "
    "--outSAMattrIHstart 0 "
    "--alignSoftClipAtReferenceEnds No "
    "--limitBAMsortRAM 300647556788 "
    "--outBAMsortingThreadN 2 "
    "--outSAMstrandField intronMotif "
    "--outFilterIntronMotifs RemoveNoncanonical "
    "--outSAMattributes All "
    "--outFilterScoreMinOverLread 0 "
    "--outFilterMatchNminOverLread 0 "
    "--outFilterMatchNmin 0 "
    "--outFilterMismatchNmax 2 "
    "--outSAMtype BAM SortedByCoordinate "
    "{log}"
)
