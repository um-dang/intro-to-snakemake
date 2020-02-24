rule all:
    input:
        expand("aligned/{sample}.sam", sample=config['samples'])

rule fasta_gunzip:
    input:
        config['ref_fasta_gz']
    output:
        config['ref_prefix'] + ".fa"
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
        sa = config['ref_prefix'] + ".sa",
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
