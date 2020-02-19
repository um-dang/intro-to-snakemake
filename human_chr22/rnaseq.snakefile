rule all:
    input:
        expand("reference/chr22.{ext}", ext=['amb', 'ann', 'bwt', 'pac', 'sa']),
        expand("coverage/{sample}.bigwig", sample=config['samples'])


rule bwa_index:
    input:
        config['ref_fasta']
    output:
        "reference/chr22.amb",
        "reference/chr22.ann",
        "reference/chr22.bwt",
        "reference/chr22.pac",
        "reference/chr22.sa"
    singularity: "docker://biocontainers/bwa:v0.7.15_cv4"
    shell:
        "bwa index -p reference/chr22 {input}"

rule bwa_mem:
    input:
        sa = "reference/chr22.sa",
        r1 = "read_data/{sample}/{sample}_R1.fastq.gz",
        r2 = "read_data/{sample}/{sample}_R2.fastq.gz"
    output:
        "aligned/{sample}.sam"
    threads: 8
    params:
        ref_prefix = "reference/chr22"
    singularity: "docker://biocontainers/bwa:v0.7.15_cv4"
    shell:
        "bwa mem -t {threads} {params.ref_prefix} {input.r1} {input.r2} > {output}"

rule sam_to_bam_bai:
    input:
        "aligned/{sample}.sam"
    output:
        bam = "aligned/{sample}.bam",
        bai = "aligned/{sample}.bai"
    singularity: "docker://biocontainers/samtools:v1.9-4-deb_cv1"
    shell:
        "samtools sort -o {output.bam} {input} && "
        "samtools index {output.bam} {output.bai}"

rule deeptools_bamCoverage:
    input:
        bam = "aligned/{sample}.bam",
        bai = "aligned/{sample}.bai"
    output:
        "coverage/{sample}.bigwig"
    singularity: "docker://quay.io/biocontainers/deeptools:3.0.1--py36_1"
    shell:
        "bamCoverage -b {input.bam} -o {output}"
