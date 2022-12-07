## Prepare the pseudogenome codes ##

configfile: '../config/config.yaml'

import pandas as pd
import os

pseudo_codes = pd.read_csv ( config['Pseudogenome_codes'])
# Create arrays of the pseudogenomes names names
pseudos = pseudo_codes['PSGs'].values


## Apply the rules ##

rule star_se_PSG1:
    input:
        fastq1="trimmed/{sample}.1.fastq",
        # path to STAR reference genome index
        index="pseudogenomes/"+config['pseudogenomes']['PSG1']+"/"+config['pseudogenomes']['PSG1']+"_GCF_001858045.2.dict",
        index2="pseudogenomes/"+config['pseudogenomes']['PSG1']+"/SAindex"
    output:
        # see STAR manual for additional output files
        file=temp("{sample}_"+pseudos[0]+"_Aligned.sortedByCoord.out.bam") #""star/{sample}/Aligned.sortedByCoord.out.sam""
    log:
        "logs/star/{sample}_PSG1.log"
    params:
        # path to STAR reference genome index
        # optional parameters
        extra="",
        filename="{sample}_"+pseudos[0]+"_",
        threads=config['threads_parallel'],
        java_opts=config['java_opts_parallel'],
    threads: config['threads_parallel']
    resources:
        mem_mb=config['mem_mb_parallel']
    conda:
        "envs/star.yaml"
    script:
        "scripts/star.py"

rule star_se_PSG2:
    input:
        fastq1="trimmed/{sample}.1.fastq",
        # path to STAR reference genome index
        index="pseudogenomes/"+config['pseudogenomes']['PSG2']+"/"+config['pseudogenomes']['PSG2']+"_GCF_001858045.2.dict",
        index2="pseudogenomes/"+config['pseudogenomes']['PSG2']+"/SAindex"
    output:
        # see STAR manual for additional output files
        file=temp("{sample}_"+pseudos[1]+"_Aligned.sortedByCoord.out.bam") #""star/{sample}/Aligned.sortedByCoord.out.sam""
    log:
        "logs/star/{sample}_PSG2.log"
    params:
        # path to STAR reference genome index
        # optional parameters
        extra="",
        filename="{sample}_"+pseudos[1]+"_",
        threads=config['threads_parallel'],
        java_opts=config['java_opts_parallel'],
    threads: config['threads_parallel']
    resources:
        mem_mb=config['mem_mb_parallel']
    conda:
        "envs/star.yaml"
    script:
        "scripts/star.py"
