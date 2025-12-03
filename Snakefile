# Directory containing FASTQ files
FASTQ_DIR = "data/data_sample"

# Output directory for BAM files
OUTDIR = "data/HiSat2_snakefile"

# Hisat2 index basename (no .1.ht2 etc.)
HISAT2_INDEX = "/pepperbase/T2T_hisat"

# SAMPLE DISCOVERY
# Look for all *_1.fastq.gz files and extract sample names
import glob
SAMPLES = [
    f.replace("_1.fastq.gz", "").split("/")[-1]
    for f in glob.glob(f"{FASTQ_DIR}/*_1.fastq.gz")
]


# BAM indices
rule all:
    input:
        expand(f"{OUTDIR}/{{sample}}.bam.bai", sample=SAMPLES)

# 1. HISAT2 MAPPING to SAM
rule hisat2:
    input:
        reads1 = lambda wildcards: f"{FASTQ_DIR}/{wildcards.sample}_1.fastq.gz",
        reads2 = lambda wildcards: f"{FASTQ_DIR}/{wildcards.sample}_2.fastq.gz",
    output:
        sam = f"{OUTDIR}/{{sample}}.sam"
    threads: 8
    params:
        index = HISAT2_INDEX
    shell:
        """
        mkdir -p {OUTDIR}
        hisat2 -p {threads} -x {params.index} \
            -1 {input.reads1} -2 {input.reads2} \
            -S {output.sam}
        """

# 2. SAM to Sorted BAM
rule sort_sam:
    input:
        sam = f"{OUTDIR}/{{sample}}.sam"
    output:
        bam = f"{OUTDIR}/{{sample}}.bam"
    shell:
        """
        samtools view -Sb {input.sam} | samtools sort -o {output.bam}
        """


# 3. BAM to BAM index
rule index_bam:
    input:
        bam = f"{OUTDIR}/{{sample}}.bam"
    output:
        bai = f"{OUTDIR}/{{sample}}.bam.bai"
    shell:
        """
        samtools index {input.bam} {output.bai}
        """
