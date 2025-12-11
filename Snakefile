# bash for all shell commands
shell.executable("/bin/bash")

# Path to conda setup script inside prog/
BLASTED_DIR = "/lustre/BIF/nobackup/BIF30806/blasted"
CONDA_SH = f"{BLASTED_DIR}/prog/etc/profile.d/conda.sh"

# Load conda + activate the environment for every rule
shell.prefix(f"source {CONDA_SH} && conda activate testenv")

# R script PATHS
# Differential Expression Analysis (Normalization and DEA)
# Takes transcript_count_matrix.csv as input and produces DEA results
DEA_SCRIPT = f"{BLASTED_DIR}/scripts/DEA_LimmaVoom.R"
# Plotting (takes DEA results as input)
# Plotting: Heatmap for Cold Condition
HEATMAP_COLD_SCRIPT = f"{BLASTED_DIR}/scripts/HeatmapCold.R"
# Plotting: Heatmap for NaCl (Salt) Condition
HEATMAP_NACL_SCRIPT = f"{BLASTED_DIR}/scripts/HeatmapNaCl.R"


# Directory containing FASTQ files
FASTQ_DIR = "data/data_sample"

# Hisat2 index basename (no .1.ht2 etc.)
HISAT2_INDEX = "pepperbase/T2T_hisat"

# Output directories
SAM_DIR = "data/HiSat2_snakefile_sam"
BAM_DIR = "data/HiSat2_snakefile_bam"
BAI_DIR = "data/HiSat2_snakefile_bai"
STRINGTIE_DIR = "data/stringtie_GTF_snakefile"
COUNT_DIR = "data/counts_stringtie_snakefile"
PREPDE_OUT_DIR = "data/prepDE_out_snakefile"

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
        expand(f"{BAI_DIR}/{{sample}}.bam.bai", sample=SAMPLES),
        expand(f"{STRINGTIE_DIR}/{{sample}}.gtf", sample=SAMPLES),
        expand(f"{COUNT_DIR}/{{sample}}/{{sample}}_quant.gtf", sample=SAMPLES),
        f"{COUNT_DIR}/merged_stringtie.gtf",
        f"{PREPDE_OUT_DIR}/gene_count_matrix.csv",
        f"{PREPDE_OUT_DIR}/transcript_count_matrix.csv"
        

# HISAT2 to produce SAM
rule hisat2:
    input:
        reads1 = lambda wc: f"{FASTQ_DIR}/{wc.sample}_1.fastq.gz",
        reads2 = lambda wc: f"{FASTQ_DIR}/{wc.sample}_2.fastq.gz",
    output:
        sam = f"{SAM_DIR}/{{sample}}.sam"
    threads: 20
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
    threads: 20
    params:
        gtf = GTF
    shell:
        """
        mkdir -p {STRINGTIE_DIR}
        stringtie -p {threads} {input.bam} -G {params.gtf} -o {output.gtf}
        """

# .txt file with GTF paths
rule gtf_list:
    input:
        expand(f"{STRINGTIE_DIR}/{{sample}}.gtf", sample=SAMPLES)
    output:
        f"{STRINGTIE_DIR}/gtf_list.txt"
    run:
        with open(output[0], "w") as f:
            for sample in SAMPLES:
                f.write(f"{STRINGTIE_DIR}/{sample}.gtf\n")

# StringTie merge mode GTF files
rule stringtie_merge:
    input:
        gtf_list = f"{STRINGTIE_DIR}/gtf_list.txt"
    output:
        merged_gtf = f"{COUNT_DIR}/merged_stringtie.gtf"
    threads: 20
    params:
        gtf = GTF
    shell:
        """
        mkdir -p {COUNT_DIR}
        stringtie --merge -p {threads} -G {params.gtf} \
            -o {output.merged_gtf} {input.gtf_list}
        """

# StringTie quantification with merged GTF
rule stringtie_quant:
    input:
        bam = f"{BAM_DIR}/{{sample}}.sorted.bam",
        merged_gtf = f"{COUNT_DIR}/merged_stringtie.gtf"
    output:
        quant_gtf = f"{COUNT_DIR}/{{sample}}/{{sample}}_quant.gtf"
    threads: 20
    run:
        import os
        os.makedirs(os.path.dirname(output.quant_gtf), exist_ok=True)
        shell(
            "stringtie -e -B -p {threads} -G {input.merged_gtf} "
            "-o {output.quant_gtf} {input.bam}"
        )

# Create input list for prepDE (No header, Space-separated)
rule create_prepde_csv:
    input:
        expand(f"{COUNT_DIR}/{{sample}}/{{sample}}_quant.gtf", sample=SAMPLES)
    output:
        # Changed extension to .txt to reflect it's not a CSV
        txt = f"{COUNT_DIR}/prepDE_input.txt" 
    run:
        with open(output.txt, "w") as f:
            # REMOVED header writing line here
            for sample in SAMPLES:
                # Changed comma to SPACE below
                f.write(f"{sample} {COUNT_DIR}/{sample}/{sample}_quant.gtf\n")

# Run prepDE.py to produce count matrices
rule prepde_counts:
    input:
        # Update input to match the new .txt output above
        txt = f"{COUNT_DIR}/prepDE_input.txt" 
    output:
        gene = f"{PREPDE_OUT_DIR}/gene_count_matrix.csv",
        transcript = f"{PREPDE_OUT_DIR}/transcript_count_matrix.csv"
    shell:
        """
        mkdir -p {PREPDE_OUT_DIR}
        prepDE.py \
            -i {input.txt} \
            -g {output.gene} \
            -t {output.transcript}
        """
