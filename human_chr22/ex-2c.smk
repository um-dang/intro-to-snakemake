rule all:
    input:
        expand("coverage/{sample}.bigwig", sample=config['samples'])

rule fasta_gunzip:
    input:
        config['ref_fasta_gz']
    output:
        temp(config['ref_prefix'] + ".fa")
    shell:
        "gunzip -c {input} > {output}"

rule bwa_index:
    input:
        config['ref_prefix'] + ".fa"
    output:
        expand(config['ref_prefix'] + ".{ext}", ext=['amb', 'ann', 'bwt', 'pac', 'sa'])
    params:
        ref_prefix = config['ref_prefix']
    singularity: "docker://biocontainers/bwa:v0.7.15_cv4"
    shell:
        "bwa index -p {params.ref_prefix} {input}"

rule bwa_mem:
    input:
        sa = "reference/chr22.sa",
        r1 = "read_data/{sample}/{sample}_R1.fastq.gz",
        r2 = "read_data/{sample}/{sample}_R2.fastq.gz"
    output:
        "aligned/{sample}.sam"
    threads: 10
    params:
        ref_prefix = config['ref_prefix']
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
