# To define samples and tools
SAMPLES = ["X"]
TOOLS = ["hisat2"]

# This will be at the end: all BAM indexes
rule all:
    input:
        expand("{sample}_{tool}.bam.bai", sample=SAMPLES, tool=TOOLS)

# 1. Map reads to reference genome with HISAT2
rule hisat2:
    input:
        reads1="data/data_sample/{sample}_1.fastq.gz",
        reads2="data/data_sample/{sample}_2.fastq.gz",
        index="/pepperbase/T2T_hisat.1.ht2"
    output:
        "{sample}_{tool}.sam"
    params:
        index="/blasted/pepperbase/T2T_hisat"
    threads: 1
    shell:
        "hisat2 -p {threads} -x {params.index} -1 {input.reads1} -2 {input.reads2} -S {output}"
        Command for one sample: hisat2 -x pepperbase/T2T_hisat -1 data/data_sample/SRR8692578_1.fastq.gz -2 data/data_sample/SRR8692578_2.fastq.gz -S SRR8692578.sam

# 2. Convert SAM to sorted BAM
rule sort_sam:
    input:
        "{sample}_{tool}.sam"
    output:
        "{sample}_{tool}.bam"
    shell:
        "samtools view -Sb {input} | samtools sort -o {output}"
        #command: samtools view -Sb SRR8692568.sam | samtools sort -o SRR8692568.sorted.bam

# 3. Index the sorted BAM
rule index_bam:
    input:
        "{sample}_{tool}.bam"
    output:
        "{sample}_{tool}.bam.bai"
    shell:
        "samtools index {input} {output}"
        #command: samtools index SRR8692568.sorted.bam SRR8692568.bai