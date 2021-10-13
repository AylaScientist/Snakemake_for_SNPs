#NOTE: Ask to launch it at 32 cores, then substitute the workflow resources to fit 32. If possible 64, even better (obviously)
samples = ["14FW20-7","16SW20-7","17SW20-7","18FW20-7","1FW20-7","2SW20-7","3SW20-7","4SW20-7"]
# ["14FW20-7","14FW19-8","14FW2-9","14FW5-8","16SW19-8","16SW20-7","16SW2-9",
# "16SW5-8","17SW19-8","17SW20-7","17SW2-9","17SW5-8","18FW19-8","18FW20-7","18FW2-9",
# "18FW5-8","1FW19-8","1FW20-7","1FW2-9","1FW5-8","2SW19-8","2SW20-7","2SW2-9","2SW5-8",
# "3SW19-8","3SW20-7","3SW2-9","3SW5-8","4SW19-8","4SW20-7","4SW2-9","4SW5-8"]

#["14FW20-7","14FW20-7"]


# Before start create the empty files ON/ON and ON with touch:
# touch ON/ON
# touch ON

rule all:
    input:
        o1="results/SNPs_function.csv",
        o2="resuts/SNP_dictionary.csv",
        o3="results/Treatment_SNPs_sig_all_for_Venn.csv",
        o4="results/Summary_of_polymorphisms.csv",
        o5="results/Intronic_SNPs_caryotype.csv",
        #html1=expand("qc/fastqc/{sample}_R1_fastqc.html", sample = samples),
        #html2=expand("qc/fastqc/{sample}_R3_fastqc.html", sample = samples),
        #html3=expand("qc/fastqc/{sample}.1_fastqc.html", sample = samples),
        #html4=expand("qc/fastqc/{sample}.2_fastqc.html", sample = samples),
        #valid_ubam=expand("logs/picard/ValidateSamFile/Validate_ubam{sample}.log", sample = samples),
        #bam=expand("recal/{sample}_ref.bam", sample = samples),
        #bqsr=expand("recal_tables/{sample}_after_ref.table", sample = samples),
        #valid_bqsr=expand("logs/picard/ValidateSamFile/Validate_bqsr_{sample}.log", sample = samples),
        #valid_bqsr_PSG1=expand("logs/picard/ValidateSamFile_PSG1/Validate_bqsr_{sample}_PSG1.log", sample = samples),
        #valid_bqsr_PSG2=expand("logs/picard/ValidateSamFile_PSG2/Validate_bqsr_{sample}_PSG2.log", sample = samples),
        #recal_table=expand("recal_tables/{sample}_after_ref.table", sample = samples),
        #degs=expand("gene_counts/{sample}.count_ref.txt", sample = samples),
        #PSG2_dict="pseudogenomes/4SW20-7/4SW20-7_GCF_001858045.2.dict",
        #PSG1_dict="pseudogenomes/14FW20-7/14FW20-7_GCF_001858045.2.dict",
        #PSG1_index="pseudogenomes/14FW20-7/SAindex",
        #PSG2_index="pseudogenomes/4SW20-7/SAindex",
        #table1="variants/AD_GT_counts_bi_PSG1.table",
        #table2="variants/AD_GT_counts_bi_PSG2.table",
        #vcf="annotated/snvs_Nile.vcf",
        #vcf="variants/AD_GT_counts_bi.table"




# Extract UMI
rule UMI_extract_R1:
    input:
        # use a list for multiple fastq files for one sample
        # usually technical replicates across lanes/flowcells example: "reads/{sample}._R2.1.fastq", "reads/{sample}._R2.2.fastq"
        fq1 = "fastq_merged/{sample}_R1.fastq",
        # paired end reads needs to be ordered so each item in the two lists match
        fq2 = "fastq_merged/{sample}_R2.fastq"
    output:
        # see STAR manual for additional output files
        pro1 = "fastq_umi_processed/{sample}_R1_umi.fastq"
    log:
        "logs/umi_extract/{sample}.log"
    params:
        # path to STAR reference genome index
        index="/home/ARO.local/aurorac/Projects/genomes/reference/",
        # optional parameters
        extra=""
    threads: 48
    shell:
        "awk -v FS=\"\t\" -v OFS=\"\t\" 'NR==FNR {split($1, id, \" \"); umi[id[1]]=$2;  next;} {split($1, id, \" \"); $1=id[1]\":\"umi[id[1]]\" \"id[2]; print $0}'  <(cat {input.fq2}|paste - - - -) <(cat {input.fq1}|paste - - - -)|tr \"\t\" \"\n\" > {output.pro1}.fastq"

rule UMI_extract_R3:
    input:
        # use a list for multiple fastq files for one sample
        # usually technical replicates across lanes/flowcells example: "reads/{sample}._R2.1.fastq", "reads/{sample}._R2.2.fastq"
        fq1 = "fastq_merged/{sample}_R3.fastq",
        # paired end reads needs to be ordered so each item in the two lists match
        fq2 = "fastq_merged/{sample}_R2.fastq"
    output:
        # see STAR manual for additional output files
        pro1 = "fastq_umi_processed/{sample}_R3_umi.fastq"
    log:
        "logs/umi_extract/{sample}.log"
    params:
        # path to STAR reference genome index
        index="/home/ARO.local/aurorac/Projects/genomes/reference/",
        # optional parameters
        extra=""
    threads: 48
    shell:
        "awk -v FS=\"\t\" -v OFS=\"\t\" 'NR==FNR {split($1, id, \" \"); umi[id[1]]=$2;  next;} {split($1, id, \" \"); $1=id[1]\":\"umi[id[1]]\" \"id[2]; print $0}'  <(cat {input.fq2}|paste - - - -) <(cat {input.fq1}|paste - - - -)|tr \"\t\" \"\n\" > {output.pro1}.fastq"

# Check the quality of the reads
rule fastqc_R1:
    input:
        "fastq_merged/{sample}_R1.fastq"
    output:
        "qc/fastqc/{sample}_R1_fastqc.html"
        #zip="qc/fastqc/{sample}_R1_umi_pro_fastqc.zip" # the suffix _fastqc.zip is necessary for multiqc to find the file. If not using multiqc, you are free to choose an arbitrary filename
    #params: "--quiet"
    log:
        "logs/fastqc/{sample}_R1.log"
    threads: 12
    script:
        "scripts/fastqc.py"

rule fastqc_R3:
    input:
        "fastq_merged/{sample}_R3.fastq"
    output:
        "qc/fastqc/{sample}_R3_fastqc.html"
        #zip="qc/fastqc/{sample}_R3_umi_pro_fastqc.zip" # the suffix _fastqc.zip is necessary for multiqc to find the file. If not using multiqc, you are free to choose an arbitrary filename
    params: "--quiet"
    log:
        "logs/fastqc/{sample}_R3.log"
    threads: 12
    script:
        "scripts/fastqc.py"


#In the trimming R3 input is called 2 for the output
rule trimmomatic_pe:
    input:
        r1="fastq_merged/{sample}_R1.fastq",
        r2="fastq_merged/{sample}_R3.fastq"
    output:
        r1=temp("trimmed/{sample}.1.fastq"),
        r2=temp("trimmed/{sample}.2.fastq"),
        # reads where trimming entirely removed the mate
        r1_unpaired=temp("trimmed/{sample}.1.unpaired.fastq"),
        r2_unpaired=temp("trimmed/{sample}.2.unpaired.fastq")
    log:
        "logs/trimmomatic/{sample}.log"
    params:
        # list of trimmers (see manual)
        trimmer="ILLUMINACLIP:TruSeq3-PE-2.fa:2:30:10:2:keepBothReads LEADING:3 TRAILING:3 MINLEN:38", # optional parameters
        #"ILLUMINACLIP:TruSeq3-PE-2.fa:2:30:10:2:keepBothReads SLIDINGWINDOW:4:5 LEADING:10 TRAILING:10 MINLEN:38",
        extra="",
        java_opts="-Xmx100g -XX:ConcGCThreads=4 -XX:ParallelGCThreads=4" #"-Xmx180g -XX:ParallelGCThreads=18"  # For Java memory specifications, please only use resources.mem_mb.
    threads: 16    # optional specification of memory usage of the JVM that snakemake will respect with global
    # resource restrictions (https://snakemake.readthedocs.io/en/latest/snakefiles/rules.html#resources)
    # and which can be used to request RAM during cluster job submission as `{resources.mem_mb}`:
    # https://snakemake.readthedocs.io/en/latest/executing/cluster.html#job-properties
    resources:
        mem_mb=100000
    script:
        "scripts/trimmomatic_pe.py"


rule fastqc_1:
    input:
        "trimmed/{sample}.1.fastq"
    output:
        "qc/fastqc/{sample}.1_fastqc.html"
        #zip="qc/fastqc/{sample}_1_trim_fastqc.zip" # the suffix _fastqc.zip is necessary for multiqc to find the file. If not using multiqc, you are free to choose an arbitrary filename
    params: "--quiet"
    log:
        "logs/fastqc/{sample}_trim_1.log"
    threads: 4
    script:
        "scripts/fastqc.py"

rule fastqc_2:
    input:
        "trimmed/{sample}.2.fastq"
    output:
        "qc/fastqc/{sample}.2_fastqc.html"
        #zip="qc/fastqc/{sample}.2_fastqc.zip" # the suffix _fastqc.zip is necessary for multiqc to find the file. If not using multiqc, you are free to choose an arbitrary filename
    params: "--quiet"
    log:
        "logs/fastqc/{sample}_trim_2.log"
    threads: 4
    script:
        "scripts/fastqc.py"


# fastq to uBAM (picard)
rule fastq_to_bam:
    input:
        fastq1="trimmed/{sample}.1.fastq",
        fastq2="trimmed/{sample}.2.fastq"
    output:
        temp("un_mapped/{sample}.bam")
    log:
        "logs/picard/sam_to_fastq/{sample}.log"
    params:
        sample="{sample}",
        extra="" # optional: Extra arguments for picard.
    # optional specification of memory usage of the JVM that snakemake will respect with global
    # resource restrictions (https://snakemake.readthedocs.io/en/latest/snakefiles/rules.html#resources)
    # and which can be used to request RAM during cluster job submission as `{resources.mem_mb}`:
    # https://snakemake.readthedocs.io/en/latest/executing/cluster.html#job-properties
    resources:
        mem_mb=50000
    threads: 4
    script:
        "scripts/fastqtosam.py"


# Annotate BAM wth UMIs (fgBio)
rule AnnotateBam:
    input:
        bam="un_mapped/{sample}.bam",
        umi="fastq_merged/{sample}_R2.fastq"
    output:
        bam=temp("un_mapped/{sample}.annotated.bam")
    params:
        java_opts="-Xmx100g -XX:+AggressiveOpts -XX:+AggressiveHeap",#-Xmx100g -XX:+AggressiveOpts -XX:+AggressiveHeap
        extra=""
    log:
        "logs/fgbio/annotate_bam/{sample}.log"
    resources:
        mem_mb=100000
    threads: 4
    script:
        "scripts/annotate_bam_with_UMIs.py"



# sort by name the unmapped bam files to be merged with the mapped bam files:
# samtools queryname sort (samtools sort -n) sorts readnames in a numerically aware fashion
# htsjdk expects queryname sorting to be done by lexicographical comparison of readnames
# This causes reads like HXXX:999 and HXXX:1000 to sort in different ways when using picard SortSam vs samtools to sort.
rule sort_ubam:
    input:
        "un_mapped/{sample}.annotated.bam"
    output:
        temp("sort_ubam/{sample}.ubam_sort.bam")
    log:
        "logs/picard/sort_ubam_byname/{sample}.log"
    params:
        sort_order="queryname",
        java_opts="-Xmx160g -XX:ConcGCThreads=4 -XX:ParallelGCThreads=4",#"-Xmx160g -XX:ConcGCThreads=4 -XX:ParallelGCThreads=4",  # For Java memory specifications, please only use resources.mem_mb.
        extra="VALIDATION_STRINGENCY=LENIENT" # optional: Extra arguments for picard.
    # optional specification of memory usage of the JVM that snakemake will respect with global
    # resource restrictions (https://snakemake.readthedocs.io/en/latest/snakefiles/rules.html#resources)
    # and which can be used to request RAM during cluster job submission as `{resources.mem_mb}`:
    # https://snakemake.readthedocs.io/en/latest/executing/cluster.html#job-properties
    resources:
         mem_mb=160000
    threads: 16
    script:
        "scripts/sort.py"



# Validate uBAMs before grouping reads_umi
rule validate_u_bam:
    input:
        bam="sort_ubam/{sample}.ubam_sort.bam"
    output:
        "logs/picard/ValidateSamFile/Validate_ubam{sample}.log"
    conda:
        "envs/mark_duplicates.yaml"
    params:
        extra = "",
        mode = "VERBOSE",
        ref="/home/ARO.local/aurorac/Projects/genomes/reference/Tilapia_header_GCF_001858045.2.fa"
    resources:
        mem_mb=160000
    threads: 4
    log:
        "logs/picard/ValidateSamFile/Validate_ubam{sample}.log"
    script:
        "scripts/validatesamfile.py"



# BAM to fastq (Picard)
# In this step the RX tag will be lost
rule bam_to_fastq:
    input:
        "sort_ubam/{sample}.ubam_sort.bam"
    output:
        fastq1=temp("reads_umi/{sample}.R1.fastq"),
        fastq2=temp("reads_umi/{sample}.R2.fastq")
    log:
        "logs/picard/sam_to_fastq/{sample}.log"
    params:
        extra="" # optional: Extra arguments for picard.
    # optional specification of memory usage of the JVM that snakemake will respect with global
    # resource restrictions (https://snakemake.readthedocs.io/en/latest/snakefiles/rules.html#resources)
    # and which can be used to request RAM during cluster job submission as `{resources.mem_mb}`:
    # https://snakemake.readthedocs.io/en/latest/executing/cluster.html#job-properties
    resources:
        mem_mb=50000
    threads: 4
    script:
        "scripts/samtofastq.py"



# map with STAR
rule star_pe:
    input:
        # use a list for multiple fastq files for one sample
        # usually technical replicates across lanes/flowcells example: "reads/{sample}._R2.1.fastq", "reads/{sample}._R2.2.fastq"
        fastq1="reads_umi/{sample}.R1.fastq",
        # paired end reads needs to be ordered so each item in the two lists match
        fastq2="reads_umi/{sample}.R2.fastq",
        index="/home/ARO.local/aurorac/Projects/genomes/reference/Tilapia_header_GCF_001858045.2.dict",
    output:
        # see STAR manual for additional output files
        file=temp("{sample}_ref_Aligned.sortedByCoord.out.bam") #""star/{sample}/Aligned.sortedByCoord.out.sam""
    log:
        "logs/star/{sample}_ref.log"
    params:
        # path to STAR reference genome index
        # optional parameters
        extra="",
        filename="{sample}_ref_",
        threads="4"
    resources:
         mem_mb=50000
    threads: 4
    script:
        "scripts/star.py"



# sort by name the unmapped bam files to be merged with the mapped bam files:
# samtools queryname sort (samtools sort -n) sorts readnames in a numerically aware fashion
# htsjdk expects queryname sorting to be done by lexicographical comparison of readnames
# This causes reads like HXXX:999 and HXXX:1000 to sort in different ways when using picard SortSam vs samtools to sort.
rule sort_mbam:
    input:
        "{sample}_ref_Aligned.sortedByCoord.out.bam"
    output:
        temp("sort_mbam/{sample}.mbam_sort.bam")
    log:
        "logs/picard/sort_mbam_byname/{sample}.log"
    params:
        sort_order="queryname",
        java_opts="-Xmx160g -XX:ConcGCThreads=4 -XX:ParallelGCThreads=4",#"-Xmx160g -XX:ConcGCThreads=4 -XX:ParallelGCThreads=4",  # For Java memory specifications, please only use resources.mem_mb.
        extra="VALIDATION_STRINGENCY=LENIENT" # optional: Extra arguments for picard.
    # optional specification of memory usage of the JVM that snakemake will respect with global
    # resource restrictions (https://snakemake.readthedocs.io/en/latest/snakefiles/rules.html#resources)
    # and which can be used to request RAM during cluster job submission as `{resources.mem_mb}`:
    # https://snakemake.readthedocs.io/en/latest/executing/cluster.html#job-properties
    resources:
        mem_mb=160000
    threads: 4
    script:
        "scripts/sort.py"



# Rule merge unmapped bam tag (including RX) with mapped bam (including coordinates to genome)
rule merge_bams:
    input:
        mbam="sort_mbam/{sample}.mbam_sort.bam",
        ubam="sort_ubam/{sample}.ubam_sort.bam",
        ref="/home/ARO.local/aurorac/Projects/genomes/reference/Tilapia_header_GCF_001858045.2.fa"
    output:
        merged=temp("merged_bam/{sample}.merged.bam")
    log:
        "logs/picard/{sample}_mergesamfiles.log"
    params:
        extra=""#"--CLIP_OVERLAPPING_READS=true" #"--SORT_ORDER=coordinate"# "-MAX_GAPS=-1 -ORIENTATIONS=FR"
    # optional specification of memory usage of the JVM that snakemake will respect with global
    # resource restrictions (https://snakemake.readthedocs.io/en/latest/snakefiles/rules.html#resources)
    # and which can be used to request RAM during cluster job submission as `{resources.mem_mb}`:
    # https://snakemake.readthedocs.io/en/latest/executing/cluster.html#job-properties
    resources:
        mem_mb=50000
    threads: 4
    script:
        "scripts/mergebamallignments.py"



##Now the pre-UMI data is processed. We can call the consensus reads and then we should map again:

# Group reads by UMI (fgBio)
rule GroupReads:
    input:
        "merged_bam/{sample}.merged.bam"
    output:
        bam=temp("umi/{sample}.gu.bam"),
        hist="umi/{sample}.gu.histo.tsv",
    params:
        java_opts="-Xmx100g -XX:+AggressiveOpts -XX:+AggressiveHeap",#-Xmx100g -XX:+AggressiveOpts -XX:+AggressiveHeap
        extra="-s adjacency -t RX" # For 30x use adjacency
    # identity: only reads with identical UMI sequences are grouped together. This strategy may be useful for evaluating data, but should generally be avoided as it will generate multiple UMI groups per original molecule in the presence of errors.
    # edit: reads are clustered into groups such that each read within a group has at least one other read in the group with <= edits differences and there are inter-group pairings with <= edits differences. Effective when there are small numbers of reads per UMI, but breaks down at very high coverage of UMIs.
    # adjacency: a version of the directed adjacency method described in umi_tools that allows for errors between UMIs but only when there is a count gradient.
    # paired: similar to adjacency but for methods that produce template with a pair of UMIs such that a read with A-B is related to but not identical to a read with B-A. Expects the pair of UMIs to be stored in a single tag, separated by a hyphen (e.g. ACGT-CCGG). The molecular IDs produced have more structure than for single UMI strategies, and are of the form {base}/{AB|BA}. E.g. two UMI pairs would be mapped as follows AAAA-GGGG -> 1/AB, GGGG-AAAA -> 1/BA.
    log:
        "logs/fgbio/group_reads/{sample}.log"
    resources:
        mem_mb=50000
    threads: 4
    script:
        "scripts/groupreadsbyumi.py"



# Call Molecular Consensus Reads (fgBio)
rule ConsensusReads:
    input:
        "umi/{sample}.gu.bam"
    output:
        temp("umi/{sample}.m1.bam")
    params:
        java_opts="-Xmx100g -XX:+AggressiveOpts -XX:+AggressiveHeap",#-Xmx100g -XX:+AggressiveOpts -XX:+AggressiveHeap
        extra="-M 1"#-M 2 if sequencing is at more than 10x
    log:
        "logs/fgbio/consensus_reads/{sample}.log"
    resources:
        mem_mb=50000
    threads: 4
    script:
        "scripts/callmolecularconsensusreads.py"

"""
# Filter Consensus Reads (fgBio)
rule FilterConsensusReads:
    input:
        "umi/{sample}.m1.bam"
    output:
        temp("umi/{sample}.m1.filtered.bam")
    params:
        extra="",
        min_base_quality=2,
        min_reads=[2, 2, 2],
        ref="/home/ARO.local/aurorac/Projects/genomes/reference/Tilapia_header_GCF_001858045.2.fa"
    log:
        "logs/fgbio/filterconsensusreads/{sample}.log"
    resources:
        mem_mb=50000
    threads: 4
    script:
        "scripts/filterconsensusreads.py"

"""


##Now we are ready to map the consensus reads for quantification of all of them
rule bam_to_fastq_UMI:
    input:
        "umi/{sample}.m1.bam"
    output:
        fastq1=temp("reads_umi/{sample}.R1_UMI.fastq"),
        fastq2=temp("reads_umi/{sample}.R2_UMI.fastq")
    log:
        "logs/picard/sam_to_fastq/{sample}.log"
    params:
        extra="" # optional: Extra arguments for picard.
    # optional specification of memory usage of the JVM that snakemake will respect with global
    # resource restrictions (https://snakemake.readthedocs.io/en/latest/snakefiles/rules.html#resources)
    # and which can be used to request RAM during cluster job submission as `{resources.mem_mb}`:
    # https://snakemake.readthedocs.io/en/latest/executing/cluster.html#job-properties
    resources:
        mem_mb=50000
    threads: 4
    script:
        "scripts/samtofastq.py"



# map with STAR
rule star_pe_UMI:
    input:
        # use a list for multiple fastq files for one sample
        # usually technical replicates across lanes/flowcells example: "reads/{sample}._R2.1.fastq", "reads/{sample}._R2.2.fastq"
        fastq1="reads_umi/{sample}.R1_UMI.fastq",
        # paired end reads needs to be ordered so each item in the two lists match
        fastq2="reads_umi/{sample}.R2_UMI.fastq",
        index="/home/ARO.local/aurorac/Projects/genomes/reference/Tilapia_header_GCF_001858045.2.dict"
    output:
        # see STAR manual for additional output files
        file=temp("{sample}_UMI_ref_Aligned.sortedByCoord.out.bam") #""star/{sample}/Aligned.sortedByCoord.out.sam""
    log:
        "logs/star/{sample}_UMI_ref.log"
    params:
        # path to STAR reference genome index
        # optional parameters
        extra="",
        filename="{sample}_UMI_ref_",
        threads="4"
    resources:
         mem_mb=50000
    threads: 4
    script:
        "scripts/star.py"

# Replace or add read Groups:
rule replace_rg:
    input:
        "{sample}_UMI_ref_Aligned.sortedByCoord.out.bam"
    output:
        temp("fixed-rg/{sample}.bam")
    conda:
        "envs/picard.yaml"
    log:
        "logs/picard/replace_rg/{sample}.log"
    params:
        extra = "SORT_ORDER=coordinate RGID=NextSeq RGLB=idp RGPL=illumina RGPU={sample} RGSM={sample} CREATE_INDEX=True",  # optional specification of memory usage of the JVM that snakemake will respect with global
    # resource restrictions (https://snakemake.readthedocs.io/en/latest/snakefiles/rules.html#resources)
    # and which can be used to request RAM during cluster job submission as `{resources.mem_mb}`:
    # https://snakemake.readthedocs.io/en/latest/executing/cluster.html#job-properties
        java_opts = "-Xmx50g -XX:ConcGCThreads=4 -XX:ParallelGCThreads=4"
    resources:
        mem_mb=50000
    threads: 16
    script:
        "scripts/readgroups.py"



# Mark and eliminate Duplicates:
rule mark_duplicates:
    input:
        "fixed-rg/{sample}.bam"
    output:
        bam=temp("marked_dedup/{sample}_ref.bam"),
        metrics="marked_dedup/{sample}_ref.metrics.txt"
    conda:
        "envs/picard.yaml"
    log:
        "logs/picard/dedup/{sample}.log"
    params:
        extra = "REMOVE_DUPLICATES=true BARCODE_TAG=RX", #Duplicates can also be removed with UMI-tools
        java_opts = "-Xmx50g -XX:ConcGCThreads=4 -XX:ParallelGCThreads=4"
    # optional specification of memory usage of the JVM that snakemake will respect with global
    # resource restrictions (https://snakemake.readthedocs.io/en/latest/snakefiles/rules.html#resources)
    # and which can be used to request RAM during cluster job submission as `{resources.mem_mb}`:
    # https://snakemake.readthedocs.io/en/latest/executing/cluster.html#job-properties
    resources:
        mem_mb=50000
    threads: 4
    script:
        "scripts/mark_duplicates.py"


# Split Ncigars
rule splitncigarreads:
    input:
        bam="marked_dedup/{sample}_ref.bam",
        ref="/home/ARO.local/aurorac/Projects/genomes/reference/Tilapia_header_GCF_001858045.2.fa"
    output:
        bam=temp("split/{sample}_ref.bam")
    conda:
        "envs/splitNcigars.yaml"
    log:
        "logs/gatk/splitNCIGARreads/{sample}.log"
    params:
        extra="",  # optional
        java_opts="-Xmx50g -XX:ConcGCThreads=4 -XX:ParallelGCThreads=4",  # "-Xmx160g -XX:ConcGCThreads=4 -XX:ParallelGCThreads=4" optional
    # optional specification of memory usage of the JVM that snakemake will respect with global
    # resource restrictions (https://snakemake.readthedocs.io/en/latest/snakefiles/rules.html#resources)
    # and which can be used to request RAM during cluster job submission as `{resources.mem_mb}`:
    # https://snakemake.readthedocs.io/en/latest/executing/cluster.html#job-properties
    resources:
        mem_mb=50000
    threads: 4
    script:
        "scripts/splitncigarreads.py"



# Base Recalibration (only for next step when there is already a first SNP call us ref)
rule gatk_baserecalibrator_before:
    input:
        bam="split/{sample}_ref.bam",
        ref="/home/ARO.local/aurorac/Projects/genomes/reference/Tilapia_header_GCF_001858045.2.fa",
        dict="/home/ARO.local/aurorac/Projects/genomes/reference/Tilapia_header_GCF_001858045.2.dict",
        known="complementary/snps_Nile.vcf" # optional known sites, in this case the SNPs from the first call before the pseudogenomes
    output:
        recal_table="recal_tables/{sample}_before_ref.table"
    conda:
        "envs/bsrc.yaml"
    log:
        "logs/gatk/baserecalibrator/{sample}.log"
    params:
        extra="",  # optional
        java_opts="-Xmx50g -XX:ConcGCThreads=4 -XX:ParallelGCThreads=4", # optional
    # optional specification of memory usage of the JVM that snakemake will respect with global
    # resource restrictions (https://snakemake.readthedocs.io/en/latest/snakefiles/rules.html#resources)
    # and which can be used to request RAM during cluster job submission as `{resources.mem_mb}`:
    # https://snakemake.readthedocs.io/en/latest/executing/cluster.html#job-properties
    resources:
        mem_mb=50000
    threads: 4
    script:
        "scripts/baserecalibrator.py"



rule gatk_applybqsr:
    input:
        bam="split/{sample}_ref.bam",
        ref="/home/ARO.local/aurorac/Projects/genomes/reference/Tilapia_header_GCF_001858045.2.fa",
        dict="/home/ARO.local/aurorac/Projects/genomes/reference/Tilapia_header_GCF_001858045.2.dict",
        recal_table="recal_tables/{sample}_before_ref.table"
    output:
        bam="recal/{sample}_ref.bam"# This bam files will be used for the DEGs analysis
    conda:
        "envs/abqsr.yaml"
    log:
        "logs/gatk/gatk_applybqsr/{sample}.log"
    params:
        extra="",  # optional
        java_opts="-Xmx50g -XX:ConcGCThreads=4 -XX:ParallelGCThreads=4", # optional
    # optional specification of memory usage of the JVM that snakemake will respect with global
    # resource restrictions (https://snakemake.readthedocs.io/en/latest/snakefiles/rules.html#resources)
    # and which can be used to request RAM during cluster job submission as `{resources.mem_mb}`:
    # https://snakemake.readthedocs.io/en/latest/executing/cluster.html#job-properties
    resources:
        mem_mb=50000
    threads: 4
    script:
        "scripts/applybqsr.py"



rule gatk_baserecalibrator_after:
    input:
        bam="recal/{sample}_ref.bam",
        ref="/home/ARO.local/aurorac/Projects/genomes/reference/Tilapia_header_GCF_001858045.2.fa",
        dict="/home/ARO.local/aurorac/Projects/genomes/reference/Tilapia_header_GCF_001858045.2.dict",
        known="complementary/snps_Nile.vcf" # optional known sites, in this case the SNPs from the first call before the pseudogenomes
    output:
        recal_table="recal_tables/{sample}_after_ref.table"
    conda:
        "envs/bsrc.yaml"
    log:
        "logs/gatk/baserecalibrator/{sample}.log"
    params:
        extra="",  # optional
        java_opts="-Xmx50g -XX:ConcGCThreads=4 -XX:ParallelGCThreads=4", # optional
    # resource restrictions (https://snakemake.readthedocs.io/en/latest/snakefiles/rules.html#resources)
    # optional specification of memory usage of the JVM that snakemake will respect with global
    # and which can be used to request RAM during cluster job submission as `{resources.mem_mb}`:
    # https://snakemake.readthedocs.io/en/latest/executing/cluster.html#job-properties
    resources:
        mem_mb=50000
    threads: 4
    script:
        "scripts/baserecalibrator.py"



# Validate bam
rule validate_bqsr:
    input:
        bam="recal/{sample}_ref.bam"
    output:
        "logs/picard/ValidateSamFile/Validate_bqsr_{sample}.log"
    conda:
        "envs/mark_duplicates.yaml"
    params:
        extra = "",
        mode = "SUMMARY",
        ref="/home/ARO.local/aurorac/Projects/genomes/reference/Tilapia_header_GCF_001858045.2.fa"
    resources:
        mem_mb=50000
    threads: 4
    log:
        "logs/picard/ValidateSamFile/Validate_bqsr_{sample}.log"
    shell:
        "validatesamfile.py"


# Make count matrix from each sample for differentially expressed genes
rule htseq:
    input:
        bam = "recal/{sample}_ref.bam",
        annotation ="/home/ARO.local/aurorac/Projects/genomes/reference/Tilapia_GCF_001858045.2_annotation.gff3"
    output:
        "gene_counts/{sample}.count_ref.txt"
    params:
        extra=" -f bam -r pos -s no -i Parent" # -n is the number of samples that process in parallel, may need to change in every experimental design.
    resources:
        mem_mb=50000
    threads: 4
    log:
        "logs/htseq/htseq_{sample}.log"
    script:
        "scripts/htseq.py"


# Haplotype caller
rule haplotype_caller:
    input:
        # single or list of bam files
        bam="recal/{sample}_ref.bam",
        ref="/home/ARO.local/aurorac/Projects/genomes/reference/Tilapia_header_GCF_001858045.2.fa"
    output:
        gvcf="calls/{sample}_ref.g.vcf"
    conda:
        "envs/haplotype.yaml"
    log:
        "logs/gatk/haplotypecaller/{sample}.log"
    params:
        extra="",  # optional
        java_opts="-Xmx50g -XX:ConcGCThreads=4 -XX:ParallelGCThreads=4", # optional
    # optional specification of memory usage of the JVM that snakemake will respect with global
    # resource restrictions (https://snakemake.readthedocs.io/en/latest/snakefiles/rules.html#resources)
    # and which can be used to request RAM during cluster job submission as `{resources.mem_mb}`:
    # https://snakemake.readthedocs.io/en/latest/executing/cluster.html#job-properties
    resources:
        mem_mb=50000
    threads: 4
    script:
        "scripts/haplotypecaller.py"



# Combine gvcf
rule combine_gvcfs:
    input:
        gvcfs=expand("calls/{sample}_ref.g.vcf", sample = samples),
        ref="/home/ARO.local/aurorac/Projects/genomes/reference/Tilapia_header_GCF_001858045.2.fa"
    output:
        gvcf=temp("calls/all_ref_g.vcf")
    conda:
        "envs/genotype.yaml"
    log:
        "logs/gatk/combinegvcfs/combinegvcfs.log"
    params:
        extra = "",
        java_opts="-Xmx400g -XX:ConcGCThreads=8 -XX:ParallelGCThreads=8",  # optional
    # optional specification of memory usage of the JVM that snakemake will respect with global
    # resource restrictions (https://snakemake.readthedocs.io/en/latest/snakefiles/rules.html#resources)
    # and which can be used to request RAM during cluster job submission as `{resources.mem_mb}`:
    # https://snakemake.readthedocs.io/en/latest/executing/cluster.html#job-properties
    resources:
        mem_mb=400000
    threads: 8
    script:
        "scripts/combinegvcfs.py"




# Genotype gvcf
rule genotype_gvcfs:
    input:
        gvcf="calls/all_ref_g.vcf",  # combined gvcf over multiple samples
        ref="/home/ARO.local/aurorac/Projects/genomes/reference/Tilapia_header_GCF_001858045.2.fa"
    output:
        vcf="calls/all_ref.vcf",
    conda:
        "envs/genotype.yaml"
    log:
        "logs/gatk/genotypegvcfs/genotypegvcfs.log"
    params:
        extra="",  # optional
        java_opts="-Xmx400g -XX:ConcGCThreads=8 -XX:ParallelGCThreads=8", # optional
    # optional specification of memory usage of the JVM that snakemake will respect with global
    # resource restrictions (https://snakemake.readthedocs.io/en/latest/snakefiles/rules.html#resources)
    # and which can be used to request RAM during cluster job submission as `{resources.mem_mb}`:
    # https://snakemake.readthedocs.io/en/latest/executing/cluster.html#job-properties
    resources:
        mem_mb=400000
    threads: 8
    script:
        "scripts/genotypegvcfs.py"



# Annotation
# Convert th vcf file to annovar format thus extracting the SNPs
rule convert_to_annovar:
    input:
        gvcf1="calls/all_ref.vcf"
    output:
        o1="annotated/annovar_Nile"
    resources:
        mem_mb=400000
    threads: 8
    shell:
        "perl scripts/convert2annovar.pl {input.gvcf1} -format vcf4 -allsample -withfreq -withfilter -context -out {output.o1}"



#Make tokens for the annotation script
rule token_annotation:
    input:
        "annotated/annovar_Nile"
    output:
        o1 = "annotated/annotated_Nile_snps"
    shell:
        "touch {output.o1} "



# Collect the annotations in the db
rule annotate:
    input:
        i1="annotated/annovar_Nile",
        i2="annotated/annotated_Nile_snps",
        buildver="ON/ON"
    output:
        o2="annotated/annotated_Nile_snps.variant_function"
    resources:
        mem_mb=400000
    threads: 8
    shell:
        "perl scripts/annotate_variation.pl -geneanno {input.i1} -buildver {input.buildver} ./ -outfile {input.i2}"





# Make tokens for the tables with the annotation:
rule token_table:
    input:
        "annotated/annotated_Nile_snps.variant_function"
    output:
        "annotated/annotated_all_snps"
    shell:
        "touch {output}"





# Annotate the file. This file can be used for the creation of the pseudogenomes
rule an_table:
    input:
        i1="annotated/annotated_Nile_snps",
        i2="annotated/annotated_Nile_snps.variant_function",
        i3="annotated/annotated_all_snps",
        gvcf1="calls/all_ref.vcf",
        path="ON/", #Path to the database (buildver)
        buildver="ON"
    output:
        o1=temp("annotated/annotated_all_snps.ON_multianno.vcf")
    resources:
        mem_mb=4000000
    threads: 8
    shell:
        "perl scripts/table_annovar.pl {input.gvcf1} {input.path} -buildver {input.buildver} -out {input.i3} -remove -protocol refGene -operation g -nastring . -vcfinput"



# Count table
# Select biallelic sites
rule gatk_select:
    input:
        vcf="annotated/annotated_all_snps.ON_multianno.vcf",
        ref="/home/ARO.local/aurorac/Projects/genomes/reference/Tilapia_header_GCF_001858045.2.fa"
    output:
        vcf="annotated/snvs_Nile.vcf"
    conda:
        "envs/genotype.yaml"
    log:
        "logs/gatk/select/snvs_Nile.log"
    params:
        extra="--restrict-alleles-to BIALLELIC",  # optional filter arguments, see GATK docs
        java_opts="-Xmx400g -XX:ConcGCThreads=8 -XX:ParallelGCThreads=8", # optional
    # optional specification of memory usage of the JVM that snakemake will respect with global
    # resource restrictions (https://snakemake.readthedocs.io/en/latest/snakefiles/rules.html#resources)
    # and which can be used to request RAM during cluster job submission as `{resources.mem_mb}`:
    # https://snakemake.readthedocs.io/en/latest/executing/cluster.html#job-properties
    resources:
        mem_mb=400000
    threads: 8
    script:
        "scripts/selectvariants.py"




# Make the table
rule gatk_variantstotable:
    input:
        vcf="annotated/snvs_Nile.vcf",
        ref="/home/ARO.local/aurorac/Projects/genomes/reference/Tilapia_header_GCF_001858045.2.fa"
    output:
        vcf="variants/AD_GT_counts_bi.table"
    conda:
        "envs/gatk.yaml"
    log:
        "logs/gatk/var2table/snvs_Nile.log"
    params:
        o1="annotated/annotated_all_snps.ON_multianno.vcf",
        extra="-SMA TRUE -F CHROM -F POS -F Gene.refGene -F Func.refGene -F ExonicFunc.refGene -F AF -GF AD -GF GT",  # optional filter arguments, see GATK docs
        java_opts=""#"-Xmx600g -XX:ConcGCThreads=8 -XX:ParallelGCThreads=8", # optional
    # optional specification of memory usage of the JVM that snakemake will respect with global
    # resource restrictions (https://snakemake.readthedocs.io/en/latest/snakefiles/rules.html#resources)
    # and which can be used to request RAM during cluster job submission as `{resources.mem_mb}`:
    # https://snakemake.readthedocs.io/en/latest/executing/cluster.html#job-properties
    resources:
        mem_mb=600000
    threads: 8
    script:
        "scripts/variantstotable.py"


## Pseudogenomes pipeline ##
############################

# From here it starts the part of the pipeline that creates the pseudogenomes
# and prepares the vcf files for the elimination of the mapping bias.

#Subset the vcf according to first and last sample
rule subset_ref_vcf:
    input:
        vcf="annotated/snvs_Nile.vcf",
    output:
        vcf="pseudogenomes/subset_vcf_file_wo_indels.recode.vcf.gz"
    resources:
        mem_mb=200000
    threads: 20
    log:
        "logs/pseudogenomes/subset_wo_indels.log"
    script:
        "scripts/subset_ref_vcf.py"



# Create the Pseudogenomes
rule pseudogenome1:
    input:
        vcf="pseudogenomes/subset_vcf_file_wo_indels.recode.vcf.gz",
        ref="/home/ARO.local/aurorac/Projects/genomes/reference/Tilapia_header_GCF_001858045.2.fa"
    output:
        pseudo = "pseudogenomes/14FW20-7/14FW20-7_GCF_001858045.2.fa"
    params:
        extra = "14FW20-7"
    resources:
        mem_mb=80000
    threads: 8
    log:
        "logs/pseudogenomes/PSG1.log"
    script:
        "scripts/pseudogenomes.py"


rule pseudogenome2:
    input:
        vcf="pseudogenomes/subset_vcf_file_wo_indels.recode.vcf.gz",
        ref="/home/ARO.local/aurorac/Projects/genomes/reference/Tilapia_header_GCF_001858045.2.fa"
    output:
        pseudo = "pseudogenomes/4SW20-7/4SW20-7_GCF_001858045.2.fa"
    params:
        extra = "4SW20-7"
    resources:
        mem_mb=80000
    threads: 8
    log:
        "logs/pseudogenomes/PSG2.log"
    script:
        "scripts/pseudogenomes.py"


# Make index for the pseudogenomes for further mapping
rule pseudo1_index:
    input:
        "pseudogenomes/14FW20-7/14FW20-7_GCF_001858045.2.fa"
    output:
        "pseudogenomes/14FW20-7/SAindex"
    params:
        threads= "16",
        dir = "pseudogenomes/14FW20-7/",
        annotation = "/home/ARO.local/aurorac/Projects/genomes/reference/Tilapia_GCF_001858045.2_annotation.gtf",
        read_length = "149"
    resources:
        mem_mb=80000
    threads: 8
    log:
        "logs/pseudogenomes/PSG1_index.log"
    script:
        "scripts/star_gi.py"


rule pseudo2_index:
    input:
        "pseudogenomes/4SW20-7/4SW20-7_GCF_001858045.2.fa"
    output:
        "pseudogenomes/4SW20-7/SAindex"
    params:
        threads= "16",
        dir = "pseudogenomes/4SW20-7/",
        annotation = "/home/ARO.local/aurorac/Projects/genomes/reference/Tilapia_GCF_001858045.2_annotation.gtf",
        read_length = "149"
    resources:
        mem_mb=80000
    threads: 8
    log:
        "logs/pseudogenomes/PSG2_index.log"
    script:
        "scripts/star_gi.py"


rule pseudo1_dict:
    input:
        "pseudogenomes/14FW20-7/14FW20-7_GCF_001858045.2.fa"
    output:
        "pseudogenomes/14FW20-7/14FW20-7_GCF_001858045.2.dict"
    resources:
        mem_mb=160000
    threads: 16
    log:
        "logs/pseudogenomes/PSG1_dictionary.log"
    script:
        "scripts/create_genome_dictionary.py"


rule pseudo2_dict:
    input:
        "pseudogenomes/4SW20-7/4SW20-7_GCF_001858045.2.fa"
    output:
        "pseudogenomes/4SW20-7/4SW20-7_GCF_001858045.2.dict"
    resources:
        mem_mb=160000
    threads: 16
    log:
        "logs/pseudogenomes/PSG2_dictionary.log"
    script:
        "scripts/create_genome_dictionary.py"


## Now you can start the pipeline again from the fastq file for the pseudogenomes:
rule star_pe_PSG1:
    input:
        # use a list for multiple fastq files for one sample
        # usually technical replicates across lanes/flowcells example: "reads/{sample}._R2.1.fastq", "reads/{sample}._R2.2.fastq"
        fastq1="reads_umi/{sample}.R1_UMI.fastq",
        # paired end reads needs to be ordered so each item in the two lists match
        fastq2="reads_umi/{sample}.R2_UMI.fastq",
        index="pseudogenomes/14FW20-7/14FW20-7_GCF_001858045.2.dict",
        index2="pseudogenomes/14FW20-7/SAindex"
    output:
        # see STAR manual for additional output files
        file=temp("{sample}_UMI_PSG1_Aligned.sortedByCoord.out.bam") #""star/{sample}/Aligned.sortedByCoord.out.sam""
    log:
        "logs/star/{sample}_UMI_PSG1.log"
    params:
        # path to STAR reference genome index
        # optional parameters
        extra="",
        filename="{sample}_UMI_PSG1_",
        threads="2"
    resources:
        mem_mb=20000
    threads: 2
    script:
        "scripts/star.py"

rule star_pe_PSG2:
    input:
        # use a list for multiple fastq files for one sample
        # usually technical replicates across lanes/flowcells example: "reads/{sample}._R2.1.fastq", "reads/{sample}._R2.2.fastq"
        fastq1="reads_umi/{sample}.R1_UMI.fastq",
        # paired end reads needs to be ordered so each item in the two lists match
        fastq2="reads_umi/{sample}.R2_UMI.fastq",
        index="pseudogenomes/4SW20-7/4SW20-7_GCF_001858045.2.dict",
        index1="pseudogenomes/4SW20-7/SAindex"
    output:
        # see STAR manual for additional output files
        file=temp("{sample}_UMI_PSG2_Aligned.sortedByCoord.out.bam") #""star/{sample}/Aligned.sortedByCoord.out.sam""
    log:
        "logs/star/{sample}_UMI_PSG2.log"
    params:
        # path to STAR reference genome index
        # optional parameters
        extra="",
        filename="{sample}_UMI_PSG2_",
        threads="2"
    resources:
        mem_mb=20000
    threads: 2
    script:
        "scripts/star.py"


# Replace or add read Groups:
rule replace_rg_PSG1:
    input:
        "{sample}_UMI_PSG1_Aligned.sortedByCoord.out.bam"
    output:
        temp("fixed-rg/{sample}_PSG1.bam")
    conda:
        "envs/picard.yaml"
    log:
        "logs/picard/replace_rg/{sample}.log"
    params:
        extra = "SORT_ORDER=coordinate RGID=NextSeq RGLB=idp RGPL=illumina RGPU={sample} RGSM={sample} CREATE_INDEX=True",  # optional specification of memory usage of the JVM that snakemake will respect with global
    # resource restrictions (https://snakemake.readthedocs.io/en/latest/snakefiles/rules.html#resources)
    # and which can be used to request RAM during cluster job submission as `{resources.mem_mb}`:
    # https://snakemake.readthedocs.io/en/latest/executing/cluster.html#job-properties
        java_opts = ""
    resources:
        mem_mb=20000
    threads: 2
    script:
        "scripts/readgroups.py"


# Mark and eliminate Duplicates:
rule mark_duplicates_PSG1:
    input:
        "fixed-rg/{sample}_PSG1.bam"
    output:
        bam=temp("marked_dedup/{sample}_PSG1.bam"),
        metrics="marked_dedup/{sample}_PSG1.metrics.txt"
    conda:
        "envs/picard.yaml"
    log:
        "logs/picard/dedup/{sample}_PSG1.log"
    params:
        extra = "REMOVE_DUPLICATES=true BARCODE_TAG=RX", #Duplicates can also be removed with UMI-tools
        java_opts = "",
    # optional specification of memory usage of the JVM that snakemake will respect with global
    # resource restrictions (https://snakemake.readthedocs.io/en/latest/snakefiles/rules.html#resources)
    # and which can be used to request RAM during cluster job submission as `{resources.mem_mb}`:
    # https://snakemake.readthedocs.io/en/latest/executing/cluster.html#job-properties
    resources:
        mem_mb=20000
    threads: 2
    script:
        "scripts/mark_duplicates.py"


# Replace or add read Groups:
rule replace_rg_PSG2:
    input:
        "{sample}_UMI_PSG2_Aligned.sortedByCoord.out.bam"
    output:
        temp("fixed-rg/{sample}_PSG2.bam")
    conda:
        "envs/picard.yaml"
    log:
        "logs/picard/replace_rg/{sample}.log"
    params:
        extra = "SORT_ORDER=coordinate RGID=NextSeq RGLB=idp RGPL=illumina RGPU={sample} RGSM={sample} CREATE_INDEX=True",  # optional specification of memory usage of the JVM that snakemake will respect with global
    # resource restrictions (https://snakemake.readthedocs.io/en/latest/snakefiles/rules.html#resources)
    # and which can be used to request RAM during cluster job submission as `{resources.mem_mb}`:
    # https://snakemake.readthedocs.io/en/latest/executing/cluster.html#job-properties
        java_opts = ""
    resources:
        mem_mb=20000
    threads: 2
    script:
        "scripts/readgroups.py"


# Mark and eliminate Duplicates:
rule mark_duplicates_PSG2:
    input:
        "fixed-rg/{sample}_PSG2.bam"
    output:
        bam=temp("marked_dedup/{sample}_PSG2.bam"),
        metrics="marked_dedup/{sample}_PSG2.metrics.txt"
    conda:
        "envs/picard.yaml"
    log:
        "logs/picard/dedup/{sample}_PSG2.log"
    params:
        extra = "REMOVE_DUPLICATES=true BARCODE_TAG=RX", #Duplicates can also be removed with UMI-tools
        java_opts = "",
    # optional specification of memory usage of the JVM that snakemake will respect with global
    # resource restrictions (https://snakemake.readthedocs.io/en/latest/snakefiles/rules.html#resources)
    # and which can be used to request RAM during cluster job submission as `{resources.mem_mb}`:
    # https://snakemake.readthedocs.io/en/latest/executing/cluster.html#job-properties
    resources:
        mem_mb=20000
    threads: 2
    script:
        "scripts/mark_duplicates.py"

# Filter?
#rule filter_consensus_reads:
#java -Xmx4g -jar fgbio.jar FilterConsensusReads \
#--input=marked_dedup/{sample}.bam \
#--output=filtered/{sample}.bam \
#--min-reads=3
#--min-base-quality=50 \
#--max-no-call-fraction=0.05



# Split Ncigars
rule splitncigarreads_PSG1:
    input:
        bam="marked_dedup/{sample}_PSG1.bam",
        ref="/home/ARO.local/aurorac/Projects/genomes/reference/Tilapia_header_GCF_001858045.2.fa"
    output:
        bam=temp("split/{sample}_PSG1.bam")
    conda:
        "envs/splitNcigars.yaml"
    log:
        "logs/gatk/splitNCIGARreads/{sample}_PSG1.log"
    resources:
        mem_mb=20000
    threads: 4
    params:
        extra="",  # optional
        java_opts="-Xmx20g -XX:ConcGCThreads=4 -XX:ParallelGCThreads=4",  # optional
    # optional specification of memory usage of the JVM that snakemake will respect with global
    # resource restrictions (https://snakemake.readthedocs.io/en/latest/snakefiles/rules.html#resources)
    # and which can be used to request RAM during cluster job submission as `{resources.mem_mb}`:
    # https://snakemake.readthedocs.io/en/latest/executing/cluster.html#job-properties
    #resources:
    #    mem_mb=5000
    script:
        "scripts/splitncigarreads.py"



# Split Ncigars
rule splitncigarreads_PSG2:
    input:
        bam="marked_dedup/{sample}_PSG2.bam",
        ref="/home/ARO.local/aurorac/Projects/genomes/reference/Tilapia_header_GCF_001858045.2.fa"
    output:
        bam=temp("split/{sample}_PSG2.bam")
    conda:
        "envs/splitNcigars.yaml"
    log:
        "logs/gatk/splitNCIGARreads/{sample}_PSG2.log"
    resources:
        mem_mb=20000
    threads: 4
    params:
        extra="",  # optional
        java_opts="-Xmx20g -XX:ConcGCThreads=4 -XX:ParallelGCThreads=4",  # optional
    # optional specification of memory usage of the JVM that snakemake will respect with global
    # resource restrictions (https://snakemake.readthedocs.io/en/latest/snakefiles/rules.html#resources)
    # and which can be used to request RAM during cluster job submission as `{resources.mem_mb}`:
    # https://snakemake.readthedocs.io/en/latest/executing/cluster.html#job-properties
    #resources:
    #    mem_mb=5000
    script:
        "scripts/splitncigarreads.py"



# Base Recalibration (only for next step when there is already a first SNP call us ref)
rule gatk_baserecalibrator_before_PSG1:
    input:
        bam="split/{sample}_PSG1.bam",
        ref="pseudogenomes/14FW20-7/14FW20-7_GCF_001858045.2.fa",
        dict="pseudogenomes/14FW20-7/14FW20-7_GCF_001858045.2.dict",
        known="complementary/snps_Nile.vcf" # optional known sites, in this case the SNPs from the first call before the pseudogenomes
    output:
        recal_table="recal_tables/{sample}_before_PSG1.table"
    conda:
        "envs/bsrc.yaml"
    log:
        "logs/gatk/baserecalibrator/{sample}_PSG1.log"
    resources:
        mem_mb=20000
    threads: 4
    params:
        extra="",  # optional
        java_opts="-Xmx20g -XX:ConcGCThreads=4 -XX:ParallelGCThreads=4", # optional
    # optional specification of memory usage of the JVM that snakemake will respect with global
    # resource restrictions (https://snakemake.readthedocs.io/en/latest/snakefiles/rules.html#resources)
    # and which can be used to request RAM during cluster job submission as `{resources.mem_mb}`:
    # https://snakemake.readthedocs.io/en/latest/executing/cluster.html#job-properties
    #resources:
    #    mem_mb=10240
    script:
        "scripts/baserecalibrator.py"


rule gatk_baserecalibrator_before_PSG2:
    input:
        bam="split/{sample}_PSG2.bam",
        ref="pseudogenomes/4SW20-7/4SW20-7_GCF_001858045.2.fa",
        dict="pseudogenomes/4SW20-7/4SW20-7_GCF_001858045.2.dict",
        known="complementary/snps_Nile.vcf" # optional known sites, in this case the SNPs from the first call before the pseudogenomes
    output:
        recal_table="recal_tables/{sample}_before_PSG2.table"
    conda:
        "envs/bsrc.yaml"
    log:
        "logs/gatk/baserecalibrator/{sample}_PSG2.log"
    resources:
        mem_mb=20000
    threads: 4
    params:
        extra="",  # optional
        java_opts="-Xmx20g -XX:ConcGCThreads=4 -XX:ParallelGCThreads=4", # optional
    # optional specification of memory usage of the JVM that snakemake will respect with global
    # resource restrictions (https://snakemake.readthedocs.io/en/latest/snakefiles/rules.html#resources)
    # and which can be used to request RAM during cluster job submission as `{resources.mem_mb}`:
    # https://snakemake.readthedocs.io/en/latest/executing/cluster.html#job-properties
    #resources:
    #    mem_mb=10240
    script:
        "scripts/baserecalibrator.py"


rule gatk_applybqsr_PSG1:
    input:
        bam="split/{sample}_PSG1.bam",
        ref="pseudogenomes/14FW20-7/14FW20-7_GCF_001858045.2.fa",
        dict="pseudogenomes/14FW20-7/14FW20-7_GCF_001858045.2.dict",
        recal_table="recal_tables/{sample}_before_PSG1.table"
    output:
        bam="recal/{sample}_PSG1.bam"# This bam files will be used for the DEGs analysis
    conda:
        "envs/abqsr.yaml"
    log:
        "logs/gatk/gatk_applybqsr/{sample}_PSG1.log"
    resources:
        mem_mb=20000
    threads: 4
    params:
        extra="",  # optional
        java_opts="-Xmx20g -XX:ConcGCThreads=4 -XX:ParallelGCThreads=4", # optional
    # optional specification of memory usage of the JVM that snakemake will respect with global
    # resource restrictions (https://snakemake.readthedocs.io/en/latest/snakefiles/rules.html#resources)
    # and which can be used to request RAM during cluster job submission as `{resources.mem_mb}`:
    # https://snakemake.readthedocs.io/en/latest/executing/cluster.html#job-properties
    #resources:
    #    mem_mb=10240
    script:
        "scripts/applybqsr.py"


rule gatk_applybqsr_PSG2:
    input:
        bam="split/{sample}_PSG2.bam",
        ref="pseudogenomes/4SW20-7/4SW20-7_GCF_001858045.2.fa",
        dict="pseudogenomes/4SW20-7/4SW20-7_GCF_001858045.2.dict",
        recal_table="recal_tables/{sample}_before_PSG2.table"
    output:
        bam="recal/{sample}_PSG2.bam"# This bam files will be used for the DEGs analysis
    conda:
        "envs/abqsr.yaml"
    log:
        "logs/gatk/gatk_applybqsr/{sample}_PSG2.log"
    resources:
        mem_mb=20000
    threads: 4
    params:
        extra="",  # optional
        java_opts="-Xmx20g -XX:ConcGCThreads=4 -XX:ParallelGCThreads=4", # optional
    # optional specification of memory usage of the JVM that snakemake will respect with global
    # resource restrictions (https://snakemake.readthedocs.io/en/latest/snakefiles/rules.html#resources)
    # and which can be used to request RAM during cluster job submission as `{resources.mem_mb}`:
    # https://snakemake.readthedocs.io/en/latest/executing/cluster.html#job-properties
    #resources:
    #    mem_mb=10240
    script:
        "scripts/applybqsr.py"


rule gatk_baserecalibrator_after_PSG1:
    input:
        bam="recal/{sample}_PSG1.bam",
        ref="pseudogenomes/14FW20-7/14FW20-7_GCF_001858045.2.fa",
        dict="pseudogenomes/14FW20-7/14FW20-7_GCF_001858045.2.dict",
        known="complementary/snps_Nile.vcf" # optional known sites, in this case the SNPs from the first call before the pseudogenomes
    output:
        recal_table="recal_tables/{sample}_after_PSG1.table"
    conda:
        "envs/bsrc.yaml"
    log:
        "logs/gatk/baserecalibrator/{sample}_PSG1.log"
    resources:
        mem_mb=20000
    threads: 4
    params:
        extra="",  # optional
        java_opts="-Xmx20g -XX:ConcGCThreads=4 -XX:ParallelGCThreads=4", # optional
    # resource restrictions (https://snakemake.readthedocs.io/en/latest/snakefiles/rules.html#resources)
    # optional specification of memory usage of the JVM that snakemake will respect with global
    # and which can be used to request RAM during cluster job submission as `{resources.mem_mb}`:
    # https://snakemake.readthedocs.io/en/latest/executing/cluster.html#job-properties
    #resources:
    #    mem_mb=10240
    script:
        "scripts/baserecalibrator.py"


rule gatk_baserecalibrator_after_PSG2:
    input:
        bam="recal/{sample}_PSG2.bam",
        ref="pseudogenomes/4SW20-7/4SW20-7_GCF_001858045.2.fa",
        dict="pseudogenomes/4SW20-7/4SW20-7_GCF_001858045.2.dict",
        known="complementary/snps_Nile.vcf" # optional known sites, in this case the SNPs from the first call before the pseudogenomes
    output:
        recal_table="recal_tables/{sample}_after_PSG2.table"
    conda:
        "envs/bsrc.yaml"
    log:
        "logs/gatk/baserecalibrator/{sample}_PSG2.log"
    resources:
        mem_mb=20000
    threads: 4
    params:
        extra="",  # optional
        java_opts="-Xmx20g -XX:ConcGCThreads=4 -XX:ParallelGCThreads=4", # optional
    # resource restrictions (https://snakemake.readthedocs.io/en/latest/snakefiles/rules.html#resources)
    # optional specification of memory usage of the JVM that snakemake will respect with global
    # and which can be used to request RAM during cluster job submission as `{resources.mem_mb}`:
    # https://snakemake.readthedocs.io/en/latest/executing/cluster.html#job-properties
    script:
        "scripts/baserecalibrator.py"


# Validate bam
rule validate_bqsr_PSG1:
    input:
        bam="recal/{sample}_PSG1.bam",
        recal_table="recal_tables/{sample}_after_PSG1.table"
    output:
        "logs/picard/ValidateSamFile_PSG1/Validate_bqsr_{sample}_PSG1.log"
    conda:
        "envs/mark_duplicates.yaml"
    params:
        extra = "",
        mode = "SUMMARY",
        ref="/home/ARO.local/aurorac/Projects/genomes/reference/Tilapia_header_GCF_001858045.2.fa"
    log:
        "logs/picard/ValidateSamFile_PSG1/Validate_bqsr_{sample}_PSG1.log"
    resources:
        mem_mb=20000
    threads: 4
    shell:
        "validatesamfile.py"

rule validate_bqsr_PSG2:
    input:
        bam="recal/{sample}_PSG2.bam",
        recal_table="recal_tables/{sample}_after_PSG2.table"
    output:
        "logs/picard/ValidateSamFile_PSG2/Validate_bqsr_{sample}_PSG2.log"
    conda:
        "envs/mark_duplicates.yaml"
    params:
        extra = "",
        mode = "SUMMARY",
        ref="/home/ARO.local/aurorac/Projects/genomes/reference/Tilapia_header_GCF_001858045.2.fa"
    log:
        "logs/picard/ValidateSamFile_PSG2/Validate_bqsr_{sample}_PSG2.log"
    resources:
        mem_mb=20000
    threads: 4
    shell:
        "validatesamfile.py"



# Haplotype caller
rule haplotype_caller_PSG1:
    input:
        # single or list of bam files
        bam="recal/{sample}_PSG1.bam",
        ref="pseudogenomes/14FW20-7/14FW20-7_GCF_001858045.2.fa",
    output:
        gvcf=temp("calls/{sample}_PSG1.g.vcf")
    conda:
        "envs/haplotype.yaml"
    log:
        "logs/gatk/haplotypecaller/{sample}_PSG1.log"
    params:
        extra="",  # optional
        java_opts="-Xmx20g -XX:ConcGCThreads=4 -XX:ParallelGCThreads=4", # optional
    # optional specification of memory usage of the JVM that snakemake will respect with global
    # resource restrictions (https://snakemake.readthedocs.io/en/latest/snakefiles/rules.html#resources)
    # and which can be used to request RAM during cluster job submission as `{resources.mem_mb}`:
    # https://snakemake.readthedocs.io/en/latest/executing/cluster.html#job-properties
    resources:
        mem_mb=20000
    threads: 4
    script:
        "scripts/haplotypecaller.py"


rule haplotype_caller_PSG2:
    input:
        # single or list of bam files
        bam="recal/{sample}_PSG2.bam",
        ref="pseudogenomes/4SW20-7/4SW20-7_GCF_001858045.2.fa"
    output:
        gvcf=temp("calls/{sample}_PSG2.g.vcf")
    conda:
        "envs/haplotype.yaml"
    log:
        "logs/gatk/haplotypecaller/{sample}_PSG2.log"
    params:
        extra="",  # optional
        java_opts="-Xmx20g -XX:ConcGCThreads=4 -XX:ParallelGCThreads=4", # optional
    # optional specification of memory usage of the JVM that snakemake will respect with global
    # resource restrictions (https://snakemake.readthedocs.io/en/latest/snakefiles/rules.html#resources)
    # and which can be used to request RAM during cluster job submission as `{resources.mem_mb}`:
    # https://snakemake.readthedocs.io/en/latest/executing/cluster.html#job-properties
    resources:
        mem_mb=20000
    threads: 4
    script:
        "scripts/haplotypecaller.py"



# Combine gvcf
rule combine_gvcfs_PSG1:
    input:
        gvcfs=expand("calls/{sample}_PSG1.g.vcf", sample = samples),
        ref="pseudogenomes/14FW20-7/14FW20-7_GCF_001858045.2.fa",
    output:
        gvcf=temp("calls/all_g_PSG1.vcf")
    conda:
        "envs/genotype.yaml"
    log:
        "logs/gatk/combinegvcfs/combinegvcfs_PSG1.log"
    params:
        extra = "",
        java_opts="-Xmx50g -XX:ConcGCThreads=4 -XX:ParallelGCThreads=4",  # optional
    # optional specification of memory usage of the JVM that snakemake will respect with global
    # resource restrictions (https://snakemake.readthedocs.io/en/latest/snakefiles/rules.html#resources)
    # and which can be used to request RAM during cluster job submission as `{resources.mem_mb}`:
    # https://snakemake.readthedocs.io/en/latest/executing/cluster.html#job-properties
    resources:
        mem_mb=160000
    threads: 16
    script:
        "scripts/combinegvcfs.py"


rule combine_gvcfs_PSG2:
    input:
        gvcfs=expand("calls/{sample}_PSG2.g.vcf", sample = samples),
        ref="pseudogenomes/4SW20-7/4SW20-7_GCF_001858045.2.fa"
    output:
        gvcf=temp("calls/all_g_PSG2.vcf")
    conda:
        "envs/genotype.yaml"
    log:
        "logs/gatk/combinegvcfs/combinegvcfs_PSG2.log"
    params:
        extra = "",
        java_opts="-Xmx160g -XX:ConcGCThreads=4 -XX:ParallelGCThreads=4",  # optional
    # optional specification of memory usage of the JVM that snakemake will respect with global
    # resource restrictions (https://snakemake.readthedocs.io/en/latest/snakefiles/rules.html#resources)
    # and which can be used to request RAM during cluster job submission as `{resources.mem_mb}`:
    # https://snakemake.readthedocs.io/en/latest/executing/cluster.html#job-properties
    resources:
        mem_mb=160000
    threads: 16
    script:
        "scripts/combinegvcfs.py"



# Genotype gvcf
rule genotype_gvcfs_PSG1:
    input:
        gvcf="calls/all_g_PSG1.vcf",  # combined gvcf over multiple samples
        ref="pseudogenomes/14FW20-7/14FW20-7_GCF_001858045.2.fa",
    output:
        vcf="calls/all_PSG1.vcf",
    conda:
        "envs/genotype.yaml"
    log:
        "logs/gatk/genotypegvcfs/genotypegvcfs_PSG1.log"
    params:
        extra="",  # optional
        java_opts="-Xmx160g -XX:ConcGCThreads=4 -XX:ParallelGCThreads=4", # optional
    # optional specification of memory usage of the JVM that snakemake will respect with global
    # resource restrictions (https://snakemake.readthedocs.io/en/latest/snakefiles/rules.html#resources)
    # and which can be used to request RAM during cluster job submission as `{resources.mem_mb}`:
    # https://snakemake.readthedocs.io/en/latest/executing/cluster.html#job-properties
    resources:
        mem_mb=160000
    threads: 16
    script:
        "scripts/genotypegvcfs.py"


rule genotype_gvcfs_PSG2:
    input:
        gvcf="calls/all_g_PSG2.vcf",  # combined gvcf over multiple samples
        ref="pseudogenomes/4SW20-7/4SW20-7_GCF_001858045.2.fa"
    output:
        vcf="calls/all_PSG2.vcf",
    conda:
        "envs/genotype.yaml"
    log:
        "logs/gatk/genotypegvcfs/genotypegvcfs_PSG2.log"
    params:
        extra="",  # optional
        java_opts="-Xmx160g -XX:ConcGCThreads=4 -XX:ParallelGCThreads=4", # optional
    # optional specification of memory usage of the JVM that snakemake will respect with global
    # resource restrictions (https://snakemake.readthedocs.io/en/latest/snakefiles/rules.html#resources)
    # and which can be used to request RAM during cluster job submission as `{resources.mem_mb}`:
    # https://snakemake.readthedocs.io/en/latest/executing/cluster.html#job-properties
    resources:
        mem_mb=160000
    threads: 16
    script:
        "scripts/genotypegvcfs.py"


# Annotation PSG1
# Convert th vcf file to annovar format thus extracting the SNPs
rule convert_to_annovar_PSG1:
    input:
        gvcf1="calls/all_PSG1.vcf"
    output:
        o1="annotated/annovar_Nile_PSG1"
    resources:
        mem_mb=160000
    threads: 16
    shell:
        "perl scripts/convert2annovar.pl {input.gvcf1} -format vcf4 -allsample -withfreq -withfilter -context -out {output.o1}"



#Make tokens for the annotation script
rule token_annotation_PSG1:
    input:
        "annotated/annovar_Nile_PSG1"
    output:
        o1 = "annotated/annotated_Nile_snps_PSG1"
    shell:
        "touch {output.o1} "



# Collect the annotations in the db
rule annotate_PSG1:
    input:
        i1="annotated/annovar_Nile_PSG1",
        i2="annotated/annotated_Nile_snps_PSG1",
        buildver="ON/ON"
    output:
        o2="annotated/annotated_Nile_snps_PSG1.variant_function"
    resources:
        mem_mb=160000
    threads: 16
    shell:
        "perl scripts/annotate_variation.pl -geneanno {input.i1} -buildver {input.buildver} ./ -outfile {input.i2}"





# Make tokens for the tables with the annotation:
rule token_table_PSG1:
    input:
        "annotated/annotated_Nile_snps_PSG1.variant_function"
    output:
        "annotated/annotated_all_snps_PSG1"
    shell:
        "touch {output}"



# Annotate the file. This file can be used for the creation of the pseudogenomes
rule an_table_PSG1:
    input:
        i1="annotated/annotated_Nile_snps_PSG1",
        i2="annotated/annotated_Nile_snps_PSG1.variant_function",
        i3="annotated/annotated_all_snps_PSG1",
        gvcf1="calls/all_PSG1.vcf",
        path="ON/", #Path to the database (buildver)
        buildver="ON"
    output:
        o1=temp("annotated/annotated_all_snps_PSG1.ON_multianno.vcf")
    resources:
        mem_mb=160000
    threads: 16
    shell:
        "perl scripts/table_annovar.pl {input.gvcf1} {input.path} -buildver {input.buildver} -out {input.i3} -remove -protocol refGene -operation g -nastring . -vcfinput"


# PSG2 annotaion
# Convert th vcf file to annovar format thus extracting the SNPs
rule convert_to_annovar_PSG2:
    input:
        gvcf1="calls/all_PSG2.vcf"
    output:
        o1="annotated/annovar_Nile_PSG2"
    resources:
        mem_mb=160000
    threads: 16
    shell:
        "perl scripts/convert2annovar.pl {input.gvcf1} -format vcf4 -allsample -withfreq -withfilter -context -out {output.o1}"



#Make tokens for the annotation script
rule token_annotation_PSG2:
    input:
        "annotated/annovar_Nile_PSG2"
    output:
        o1 = "annotated/annotated_Nile_snps_PSG2"
    shell:
        "touch {output.o1} "



# Collect the annotations in the db
rule annotate_PSG2:
    input:
        i1="annotated/annovar_Nile_PSG2",
        i2="annotated/annotated_Nile_snps_PSG2",
        buildver="ON/ON"
    output:
        o2="annotated/annotated_Nile_snps_PSG2.variant_function"
    resources:
        mem_mb=160000
    threads: 16
    shell:
        "perl scripts/annotate_variation.pl -geneanno {input.i1} -buildver {input.buildver} ./ -outfile {input.i2}"





# Make tokens for the tables with the annotation:
rule token_table_PSG2:
    input:
        "annotated/annotated_Nile_snps_PSG2.variant_function"
    output:
        "annotated/annotated_all_snps_PSG2"
    shell:
        "touch {output}"



# Annotate the file. This file can be used for the creation of the pseudogenomes
rule an_table_PSG2:
    input:
        i1="annotated/annotated_Nile_snps_PSG2",
        i2="annotated/annotated_Nile_snps_PSG2.variant_function",
        i3="annotated/annotated_all_snps_PSG2",
        gvcf1="calls/all_PSG2.vcf",
        path="ON/", #Path to the database (buildver)
        buildver="ON"
    output:
        o1=temp("annotated/annotated_all_snps_PSG2.ON_multianno.vcf")
    resources:
        mem_mb=160000
    threads: 16
    shell:
        "perl scripts/table_annovar.pl {input.gvcf1} {input.path} -buildver {input.buildver} -out {input.i3} -remove -protocol refGene -operation g -nastring . -vcfinput"



# Count table
# Select biallelic sites
rule gatk_select_PSG1:
    input:
        vcf="annotated/annotated_all_snps_PSG1.ON_multianno.vcf",
        ref="pseudogenomes/14FW20-7/14FW20-7_GCF_001858045.2.fa",
    output:
        vcf="annotated/snvs_Nile_PSG1.vcf"
    conda:
        "envs/genotype.yaml"
    log:
        "logs/gatk/select/snvs_Nile_PSG1.log"
    params:
        extra="--restrict-alleles-to BIALLELIC",  # optional filter arguments, see GATK docs
        java_opts="-Xmx160g -XX:ConcGCThreads=16 -XX:ParallelGCThreads=16", # optional
    # optional specification of memory usage of the JVM that snakemake will respect with global
    # resource restrictions (https://snakemake.readthedocs.io/en/latest/snakefiles/rules.html#resources)
    # and which can be used to request RAM during cluster job submission as `{resources.mem_mb}`:
    # https://snakemake.readthedocs.io/en/latest/executing/cluster.html#job-properties
    resources:
        mem_mb=160000
    threads: 16
    script:
        "scripts/selectvariants.py"

rule gatk_select_PSG2:
    input:
        vcf="annotated/annotated_all_snps_PSG2.ON_multianno.vcf",
        ref="pseudogenomes/4SW20-7/4SW20-7_GCF_001858045.2.fa"
    output:
        vcf="annotated/snvs_Nile_PSG2.vcf"
    conda:
        "envs/genotype.yaml"
    log:
        "logs/gatk/select/snvs_Nile_PSG2.log"
    params:
        extra="--restrict-alleles-to BIALLELIC",  # optional filter arguments, see GATK docs
        java_opts="-Xmx160g -XX:ConcGCThreads=16 -XX:ParallelGCThreads=16", # optional
    # optional specification of memory usage of the JVM that snakemake will respect with global
    # resource restrictions (https://snakemake.readthedocs.io/en/latest/snakefiles/rules.html#resources)
    # and which can be used to request RAM during cluster job submission as `{resources.mem_mb}`:
    # https://snakemake.readthedocs.io/en/latest/executing/cluster.html#job-properties
    resources:
        mem_mb=160000
    threads: 16
    script:
        "scripts/selectvariants.py"




# Make the table
rule gatk_variantstotable_PSG1:
    input:
        vcf="annotated/snvs_Nile_PSG1.vcf",
        ref="pseudogenomes/14FW20-7/14FW20-7_GCF_001858045.2.fa"
    output:
        vcf="variants/AD_GT_counts_bi_PSG1.table"
    conda:
        "envs/gatk.yaml"
    log:
        "logs/gatk/var2table/snvs_Nile_PSG1.log"
    params:
        o1="annotated/annotated_all_snps_PSG1.ON_multianno.vcf",
        extra="-SMA TRUE -F CHROM -F POS -F Gene.refGene -F Func.refGene -F ExonicFunc.refGene -F AF -GF AD -GF GT",  # optional filter arguments, see GATK docs
        java_opts="-Xmx160g -XX:ConcGCThreads=16 -XX:ParallelGCThreads=16", # optional
    # optional specification of memory usage of the JVM that snakemake will respect with global
    # resource restrictions (https://snakemake.readthedocs.io/en/latest/snakefiles/rules.html#resources)
    # and which can be used to request RAM during cluster job submission as `{resources.mem_mb}`:
    # https://snakemake.readthedocs.io/en/latest/executing/cluster.html#job-properties
    resources:
        mem_mb=160000
    threads: 16
    script:
        "scripts/variantstotable.py"



rule gatk_variantstotable_PSG2:
    input:
        vcf="annotated/snvs_Nile_PSG2.vcf",
        ref="pseudogenomes/4SW20-7/4SW20-7_GCF_001858045.2.fa"
    output:
        vcf="variants/AD_GT_counts_bi_PSG2.table"
    conda:
        "envs/gatk.yaml"
    log:
        "logs/gatk/var2table/snvs_Nile_PSG2.log"
    params:
        o1="annotated/annotated_all_snps_PSG2.ON_multianno.vcf",
        extra="-SMA TRUE -F CHROM -F POS -F Gene.refGene -F Func.refGene -F ExonicFunc.refGene -F AF -GF AD -GF GT",  # optional filter arguments, see GATK docs
        java_opts="-Xmx160g -XX:ConcGCThreads=16 -XX:ParallelGCThreads=16", # optional
    # optional specification of memory usage of the JVM that snakemake will respect with global
    # resource restrictions (https://snakemake.readthedocs.io/en/latest/snakefiles/rules.html#resources)
    # and which can be used to request RAM during cluster job submission as `{resources.mem_mb}`:
    # https://snakemake.readthedocs.io/en/latest/executing/cluster.html#job-properties
    resources:
        mem_mb=160000
    threads: 16
    script:
        "scripts/variantstotable.py"




## Workflow ##
##############

# Preapre the tables for the workflow format
rule table_step1_PSG1:
    input:
        table1="variants/AD_GT_counts_bi_PSG1.table"
    output:
        table1="variants/AD_GT_counts_bi_PSG1_step1.csv"
    shell:
        "sed 's/|/\//g' {input.table1} > {output.table1}"


rule table_step1_PSG2:
    input:
        table1="variants/AD_GT_counts_bi_PSG2.table"
    output:
        table1="variants/AD_GT_counts_bi_PSG2_step1.csv"
    shell:
        "sed 's/|/\//g' {input.table1} > {output.table1}"




rule table_step2_PSG1:
    input:
        table="variants/AD_GT_counts_bi_PSG1_step1.csv",
        tb_colnames = "complementary/tb1_colnames.csv",
    output:
        csv="workflow/AD_GT_counts_bi_PSG1.csv"
    resources:
        mem_mb=240000
    threads: 24
    conda:
        "envs/python.yaml"
    log:
        "logs/python/table_step2/AD_GT_counts_bi_PSG1.log"
    script:
        "scripts/table2df_step2.py"


rule table_step2_PSG2:
    input:
        table="variants/AD_GT_counts_bi_PSG2_step1.csv",
        tb_colnames = "complementary/tb2_colnames.csv",
    output:
        csv="workflow/AD_GT_counts_bi_PSG2.csv"
    resources:
        mem_mb=240000
    threads: 24
    conda:
        "envs/python.yaml"
    log:
        "logs/python/table_step2/AD_GT_counts_bi_PSG2.log"
    script:
        "scripts/table2df_step2.py"


# Workflow_ASE for merging vcf files (5 rules)

rule merge_dataframes:
    input:
        csv1="workflow/AD_GT_counts_bi_PSG1.csv",
        csv2="workflow/AD_GT_counts_bi_PSG2.csv"
    output:
        "workflow/merged_df.csv"
    params:
        extra="",
        java_opts=""
    resources:
        mem_mb=240000
    threads:24
    conda:
        "envs/python.yaml"
    log:
        "logs/workflow/merge_dataframes.log"
    script:
        "scripts/mergedataframes.py"


rule biallelic_sites:
    input:
        csv="workflow/merged_df.csv",
        #sn1="complementary/Sample_names_complete_project3_Nile.csv",
        sn1="complementary/Sample_names.csv",
        psc="complementary/Pseudogenome_codes.csv"
    output:
        "workflow/biallelic.csv"
    params:
        extra="",
        java_opts=""
    resources:
        mem_mb=240000
    threads:24
    conda:
        "envs/python.yaml"
    log:
        "logs/workflow/biallelic_sites.log"
    script:
        "scripts/biallelicsites.py"


rule average:
    input:
        csv="workflow/biallelic.csv",
        sn1="complementary/Sample_names.csv",
        psc="complementary/Pseudogenome_codes.csv"
    output:
        "workflow/average.csv"
    params:
        extra="",
        java_opts=""
    resources:
        mem_mb=240000
    threads: 24
    conda:
        "envs/python.yaml"
    log:
        "logs/workflow/average.log"
    script:
        "scripts/average.py"



rule AD10:
    input:
        csv="workflow/average.csv",
        sn1="complementary/Sample_names.csv"
    output:
        "workflow/ad10.csv"
    params:
        extra="",
        java_opts=""
    resources:
        mem_mb=240000
    threads: 24
    conda:
        "envs/python.yaml"
    log:
        "logs/workflow/ad10.log"
    script:
        "scripts/ad10.py"


rule genotype:
    input:
        csv="workflow/ad10.csv",
        sn1="complementary/Sample_names.csv"
    output:
        result1="results/Uniform_to_validate_SNPs.csv",
        result2="results/SNP_panel.csv", #This SNP panel contains the SNPs found expressed in the experimental popilation
        result3="results/Genotype_models.csv"
    params:
        extra="",
        java_opts=""
    resources:
        mem_mb=240000
    threads: 24
    conda:
        "envs/python.yaml"
    log:
        "logs/workflow/genotype.log"
    script:
        "scripts/genotype.py"

# Here you can include a step to separate monoallelic expression (MAE) in case you have the genotypes from the genome
# The monoallelic expression (MAE) can result from homozygots as well as imprinted genes. In order to distinguish each
# case, we need to know the genotype.

"""
rule MAE:
    input:
        result1="results/Uniform_to_validate_SNPs.csv",
        geno="complementary/known_genotypes.csv"
    output:
        result4="workflow/SNPs_ready.csv"
    log:
        "logs/python/ASE_data_wrangling/SNPs_ready.log"
    conda:
        "envs/python.yaml"
    resources:
        mem_mb=240000
    threads:24
    script:
        "scripts/ASE_data_wrangling.py"
"""

# Quality control of t 773he output data in the workflow
rule QC:
    input:
        i1="results/Uniform_to_validate_SNPs.csv",
        sn1="complementary/Sample_names.csv",
        i2="complementary/Experimental_groups.csv"
    output:
        results="results/QC.csv"
    conda:
        "envs/python.yaml"
    resources:
        mem_mb=240000
    threads: 24
    log:
        "logs/python/QC/QC.log"
    script:
        "scripts/QC.py"


# Stats on correct data
rule stats:
    input:
        i1="complementary/Experimental_design.csv",
        i2="complementary/Experimental_groups.csv",
        i3="complementary/Sample_names.csv",
        i4="results/QC.csv"
    output:
        o1="stats/Chi_test.csv",
        o2="stats/Fisher_test.csv",
        o3="results/SNPs_analysed.csv"
    log:
        "logs/python/stats/SNPs_analysed.log"
    conda:
        "envs/python.yaml"
    resources:
        mem_mb=240000
    threads: 24
    script:
        "scripts/Stats.py"


# Tables and figures for the publication
rule results:
    input:
        i1="complementary/Tilapia gene transcript protein GO KEGG uniprot.csv",
        i2="complementary/GO_dictionary.csv",
        i3="results/SNPs_analysed.csv",
        i4="complementary/Sample_names.csv",
        i5="complementary/Experimental_groups.csv",
        i6="complementary/Experimental_design.csv",
        i7="complementary/GeneID_to_Genebank.csv"
    output:
        o1="results/SNPs_function.csv",
        o2="resuts/SNP_dictionary.csv",
        o3="results/Treatment_SNPs_sig_all_for_Venn.csv",
        o4="results/Summary_of_polymorphisms.csv",
        o5="results/Intronic_SNPs_caryotype.csv"
    log:
        "logs/python/results/tables.log"
    conda:
        "envs/python.yaml"
    resources:
        mem_mb=240000
    threads: 24
    script:
        "scripts/Tables.py"

"""
# Complementary analysis with DEGs:
rule DE_analysis:
    input:
        i1="results/SNPs_analysed.csv",
        i2="complementary/Test_names.csv" # These names are a subset of the table Experimental_design.csv and corresponds to the tests common to ASE and DEGs
    output:
        "degs/SNPs_DEGs.csv"
    conda:
        "envs/python.yaml"
    shell:
        "python scripts/DE_genes.py"
"""
