# Directory containing FASTQ files
FASTQ_DIR = "data/data_sample"

# Hisat2 index basename (no .1.ht2 etc.)
HISAT2_INDEX = "pepperbase/T2T_hisat"

# Output directories
SAM_DIR = "data/HiSat2_snakefile_sam"
BAM_DIR = "data/HiSat2_snakefile_bam"
BAI_DIR = "data/HiSat2_snakefile_bai"
STRINGTIE_DIR = "data/stringtie_GTF_snakefile"

# GTF/GFF file for StringTie
GTF = "pepperbase/capsicum_genome.gff"

# SAMPLE DISCOVERY: look for all *_1.fastq.gz files and extract sample names
import glob
SAMPLES = [
    f.replace("_1.fastq.gz", "").split("/")[-1]
    for f in glob.glob(f"{FASTQ_DIR}/*_1.fastq.gz")
]


rule all:
    input:
        expand(f"{BAI_DIR}/{{sample}}.bam.bai", sample=SAMPLES) +
        expand(f"{STRINGTIE_DIR}/{{sample}}.gtf", sample=SAMPLES)

# HISAT2 to produce SAM
rule hisat2:
    input:
        reads1 = lambda wc: f"{FASTQ_DIR}/{wc.sample}_1.fastq.gz",
        reads2 = lambda wc: f"{FASTQ_DIR}/{wc.sample}_2.fastq.gz",
    output:
        sam = f"{SAM_DIR}/{{sample}}.sam"
    threads: 32
    params:
        index = HISAT2_INDEX
    shell:
        """
        mkdir -p {SAM_DIR}
        hisat2 -p {threads} -x {params.index} \
            -1 {input.reads1} -2 {input.reads2} \
            -S {output.sam}
        """

# SAM to sorted BAM
rule sam_to_bam:
    input:
        sam = f"{SAM_DIR}/{{sample}}.sam"
    output:
        bam = f"{BAM_DIR}/{{sample}}.sorted.bam"
    shell:
        """
        mkdir -p {BAM_DIR}
        samtools view -Sb {input.sam} | samtools sort -o {output.bam}
        """

# BAM to BAM index (.bam.bai)
rule index_bam:
    input:
        bam = f"{BAM_DIR}/{{sample}}.sorted.bam"
    output:
        bai = f"{BAI_DIR}/{{sample}}.bam.bai"
    shell:
        """
        mkdir -p {BAI_DIR}
        samtools index {input.bam} {output.bai}
        """

# StringTie: assemble / quantify using sorted BAMs
rule stringtie:
    input:
        bam = f"{BAM_DIR}/{{sample}}.sorted.bam"
    output:
        gtf = f"{STRINGTIE_DIR}/{{sample}}.gtf"
    threads: 32
    params:
        gtf = GTF
    shell:
        """
        mkdir -p {STRINGTIE_DIR}
        stringtie -p {threads} {input.bam} -G {params.gtf} -o {output.gtf}
        """

