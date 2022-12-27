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


if fq1[0].endswith(".gz"):
    readcmd = "--readFilesCommand zcat"
else:
    readcmd = ""



shell(
    "STAR "
    "{extra} "
    "--runThreadN {threads} "
    "--genomeDir {index_dir} "
    "--readFilesIn {fq1} " 
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
